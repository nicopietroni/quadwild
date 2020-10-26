#include "quad_from_patches.h"

#include "includes/orient_faces.h"
#include "includes/ilp.h"
#include "includes/patterns.h"
#include "includes/mapping.h"
#include "includes/convert.h"

#ifdef SAVEMESHESFORDEBUG
#include <igl/writeOBJ.h>
#include <wrap/io_trimesh/export_obj.h>
#endif


namespace qfp {

template<class PolyMesh, class TriangleMesh>
void quadrangulationFromPatches(
    TriangleMesh& trimesh,
    const std::vector<std::vector<size_t>>& trimeshPartitions,
    const std::vector<std::vector<size_t>>& trimeshCorners,
    const std::vector<double>& edgeFactor,
    const Parameters& parameters,
    PolyMesh& quadmesh,
    std::vector<std::vector<size_t>>& quadmeshPartitions,
    std::vector<std::vector<size_t>>& quadmeshCorners,
    std::vector<int>& ilpResult)
{
    ChartData chartData;

    assert(trimeshPartitions.size() == trimeshCorners.size() && edgeFactor.size() == trimeshPartitions.size());

    //Get chart data
    chartData = getChartData(
            trimesh,
            trimeshPartitions,
            trimeshCorners);

    //Solve ILP to find best side size
    ilpResult = findSubdivisions(
            chartData,
            edgeFactor,
            parameters.ilpMethod,
            parameters.alpha,
            parameters.isometry,
            parameters.regularityForQuadrilaterals,
            parameters.regularityForNonQuadrilaterals,
            parameters.regularityNonQuadrilateralWeight,
            parameters.hardParityConstraint,
            parameters.timeLimit,
            parameters.gapLimit,
            parameters.minimumGap);

    //Quadrangulate
    quadrangulate(
            trimesh,
            chartData,
            ilpResult,
            parameters.chartSmoothingIterations,
            parameters.quadrangulationSmoothingIterations,
            quadmesh,
            quadmeshPartitions,
            quadmeshCorners);
}

template<class TriangleMesh>
ChartData getChartData(
        TriangleMesh& trimesh,
        const std::vector<std::vector<size_t>>& trimeshPartitions,
        const std::vector<std::vector<size_t>>& trimeshCorners)
{
    std::vector<int> faceLabel(trimesh.face.size(), -1);
    for (size_t pId = 0; pId < trimeshPartitions.size(); pId++) {
        for (const size_t& fId : trimeshPartitions[pId]) {
            assert(faceLabel[fId] == -1);
            faceLabel[fId] = static_cast<int>(pId);
        }
    }

    ChartData chartData = getPatchDecompositionChartData(trimesh, faceLabel, trimeshCorners);
    return chartData;
}

inline std::vector<int> findSubdivisions(
        const ChartData& chartData,
        const std::vector<double>& edgeFactor,
        const ILPMethod& method,
        const double alpha,
        const bool isometry,
        const bool regularityForQuadrilaterals,
        const bool regularityForNonQuadrilaterals,
        const double regularityNonQuadrilateralWeight,
        const bool hardParityConstraint,
        const double timeLimit,
        const double gapLimit,
        const double minimumGap)
{
    if (chartData.charts.size() <= 0)
        return std::vector<int>();

    double gap;
    ILPStatus status;

//    //Fix the chart on border
//    for (ChartSubSide& subside : chartData.subSides) {
//        subside.isFixed = subside.isOnBorder;
//    }

    //Solve ILP to find the best patches
    std::vector<int> ilpResult = solveILP(chartData, edgeFactor, method, alpha, isometry, regularityForQuadrilaterals, regularityForNonQuadrilaterals, regularityNonQuadrilateralWeight, hardParityConstraint, timeLimit, gapLimit, gap, status);

    if (status == ILPStatus::SOLUTIONFOUND && gap < minimumGap) {
        std::cout << "Solution found! Gap: " << gap << std::endl;
    }
    else {
        if (status == ILPStatus::INFEASIBLE) {
            std::cout << "Error! Model was infeasible or time limit exceeded!" << std::endl;
        }
        else if (status == ILPStatus::SOLUTIONWRONG && !hardParityConstraint) {
            std::cout << "Solution wrong! It have been used soft constraints for parity, so trying with hard constraints: " << std::endl;
            return findSubdivisions(chartData, edgeFactor, method, alpha, isometry, regularityForQuadrilaterals, regularityForNonQuadrilaterals, regularityNonQuadrilateralWeight, true, timeLimit, gapLimit, minimumGap);
        }
        else {
            if (method == ILPMethod::LEASTSQUARES) {
                std::cout << "Minimum gap has been not reached. Trying with ABS (linear optimization method)." << gap << std::endl;
                return findSubdivisions(chartData, edgeFactor, ILPMethod::ABS, alpha, isometry, regularityForQuadrilaterals, regularityForNonQuadrilaterals, regularityNonQuadrilateralWeight, hardParityConstraint, timeLimit, gapLimit, minimumGap);
            }
            else if (regularityForNonQuadrilaterals) {
                std::cout << "Minimum gap has been not reached. Trying without regularity for non-quadrilaterals." << gap << std::endl;
                return findSubdivisions(chartData, edgeFactor, ILPMethod::ABS, alpha, isometry, regularityForQuadrilaterals, false, regularityNonQuadrilateralWeight, hardParityConstraint, timeLimit, gapLimit, minimumGap);
            }
            else if (regularityForQuadrilaterals) {
                std::cout << "Minimum gap has been not reached. Trying without any regularity terms." << gap << std::endl;
                return findSubdivisions(chartData, edgeFactor, ILPMethod::ABS, alpha, true, false, false, regularityNonQuadrilateralWeight, hardParityConstraint, timeLimit, gapLimit, minimumGap);
            }
        }
    }

    return ilpResult;
}


template<class TriangleMesh, class PolyMesh>
void quadrangulate(
        TriangleMesh& trimesh,
        const ChartData& chartData,
        const std::vector<int>& ilpResult,
        const int chartSmoothingIterations,
        const int quadrangulationSmoothingIterations,
        PolyMesh& quadmesh,
        std::vector<std::vector<size_t>>& quadmeshPartitions,
        std::vector<std::vector<size_t>>& quadmeshCorners)
{
    if (trimesh.face.size() <= 0)
        return;
    if (ilpResult.size() == 0)
        return;

    std::vector<std::vector<size_t>> vertexSubsideMap(chartData.subSides.size());
    std::vector<int> cornerVertices(trimesh.vert.size(), -1);
    std::vector<size_t> fixedVertices;

    quadmeshPartitions.resize(chartData.charts.size());
    quadmeshCorners.resize(chartData.charts.size());   

    //Fill fixed vertices (subsides corners)
    for (const ChartSubSide& subside : chartData.subSides) {
        size_t vStart = subside.vertices[0];
        size_t vEnd = subside.vertices[subside.vertices.size() - 1];

        if (cornerVertices[vStart] == -1) {
            cornerVertices[vStart] = quadmesh.vert.size();
            vcg::tri::Allocator<PolyMesh>::AddVertex(
                        quadmesh,
                        trimesh.vert.at(vStart).P());

            fixedVertices.push_back(cornerVertices[vStart]);
        }

        if (cornerVertices[vEnd] == -1) {
            cornerVertices[vEnd] = quadmesh.vert.size();
            vcg::tri::Allocator<PolyMesh>::AddVertex(
                        quadmesh,
                        trimesh.vert.at(vEnd).P());

            fixedVertices.push_back(cornerVertices[vEnd]);
        }
    }

#ifdef SAVEMESHESFORDEBUG
        vcg::tri::io::ExporterOBJ<PolyMesh>::Save(quadmesh, "res/corner_vertices.obj", vcg::tri::io::Mask::IOM_NONE);
#endif

    //Fill subside map for fixed borders
    std::set<size_t> finalMeshBorders;
    for (size_t i = 0; i < chartData.subSides.size(); i++) {
        const ChartSubSide& subside = chartData.subSides[i];
        if (subside.isFixed) {
            for (size_t k = 0; k < subside.vertices.size(); k++) {
                const size_t& vId = subside.vertices.at(k);

                size_t newVertexId;

                if (cornerVertices[vId] == -1) {
                    assert(k > 0 && k < subside.vertices.size() - 1);

                    newVertexId = quadmesh.vert.size();
                    vcg::tri::Allocator<PolyMesh>::AddVertex(
                                quadmesh,
                                trimesh.vert.at(vId).P());

                    fixedVertices.push_back(newVertexId);
                }
                else {
                    newVertexId = cornerVertices[vId];
                    assert(newVertexId >= 0);
                }

                finalMeshBorders.insert(newVertexId);

                vertexSubsideMap[i].push_back(newVertexId);
            }
        }
    }

#ifdef SAVEMESHESFORDEBUG
        vcg::tri::io::ExporterOBJ<PolyMesh>::Save(quadmesh, "res/fixed_vertices.obj", vcg::tri::io::Mask::IOM_NONE);
#endif


    //For each chart
    for (size_t cId = 0; cId < chartData.charts.size(); cId++) {
        const Chart& chart = chartData.charts[cId];

        if (chart.faces.size() == 0)
            continue;

        const std::vector<ChartSide>& chartSides = chart.chartSides;
        if (chartSides.size() < 3 || chartSides.size() > 6) {
            std::cout << "Chart " << cId << " with corners less than 3 or greater than 6!" << std::endl;
            continue;
        }

        bool ilpSolvedForAll = true;
        for (size_t sId : chart.chartSubSides) {
            if (ilpResult[sId] < 0)
                ilpSolvedForAll = false;
        }

        if (!ilpSolvedForAll) {
            std::cout << "Chart " << cId << " not computed. ILP was not solved." << std::endl;
            continue;
        }

        //Input mesh
        Eigen::MatrixXd chartV;
        Eigen::MatrixXi chartF;
        vcg::tri::UpdateFlags<TriangleMesh>::FaceClearS(trimesh);
        vcg::tri::UpdateFlags<TriangleMesh>::VertexClearS(trimesh);
        for (const size_t& fId : chart.faces) {
            trimesh.face[fId].SetS();
            for (int k = 0; k < trimesh.face[fId].VN(); k++) {
                trimesh.face[fId].V(k)->SetS();
            }
        }
        std::vector<int> vMap, fMap;
        VCGToEigen(trimesh, chartV, chartF, vMap, fMap, true, 3);

        //Input subdivisions
        Eigen::VectorXi l(chartSides.size());

        std::vector<std::vector<double>> chartSideLength(chartSides.size());
        std::vector<std::vector<std::vector<size_t>>> chartSideVertices(chartSides.size());
        std::vector<std::vector<size_t>> chartSideSubdivision(chartSides.size());

        for (size_t i = 0; i < chartSides.size(); i++) {
            const ChartSide& chartSide = chartSides[i];

            chartSideLength[i].resize(chartSide.subsides.size());
            chartSideVertices[i].resize(chartSide.subsides.size());
            chartSideSubdivision[i].resize(chartSide.subsides.size());

            size_t targetSideSubdivision = 0;
            for (size_t j = 0; j < chartSide.subsides.size(); j++) {
                const size_t& subSideId = chartSides[i].subsides[j];
                const ChartSubSide& subSide = chartData.subSides[subSideId];

                targetSideSubdivision += ilpResult[subSideId];

                chartSideLength[i][j] = subSide.length;
                chartSideVertices[i][j] = subSide.vertices;
                chartSideSubdivision[i][j] = ilpResult[subSideId];

                if (chartSide.reversedSubside[j]) {
                    std::reverse(chartSideVertices[i][j].begin(), chartSideVertices[i][j].end());
                }

                for (size_t k = 0; k < chartSideVertices[i][j].size(); k++) {
                    size_t vId = chartSideVertices[i][j][k];
                    assert(vMap[vId] >= 0);
                    chartSideVertices[i][j][k] = vMap[vId];
                }

                if (ilpResult[subSideId] < 0) {
                    std::cout << "Warning: ILP not valid" << std::endl;
                    return;
                }
            }

            l(static_cast<int>(i)) = targetSideSubdivision;
        }

        //Pattern quadrangulation
        Eigen::MatrixXd patchV;
        Eigen::MatrixXi patchF;
        std::vector<size_t> patchBorders;
        std::vector<size_t> patchCorners;
        qfp::computePattern(l, patchV, patchF, patchBorders, patchCorners);


#ifdef SAVEMESHESFORDEBUG
        igl::writeOBJ(std::string("res/") + std::to_string(cId) + std::string("_patch.obj"), patchV, patchF);
#endif

        std::vector<std::vector<size_t>> patchEigenSides = qfp::getPatchSides(patchV, patchF, patchBorders, patchCorners, l);

#ifdef SAVEMESHESFORDEBUG
        igl::writeOBJ(std::string("res/") + std::to_string(cId) + std::string("_patch_adjusted.obj"), patchV, patchF);
#endif

        assert(chartSides.size() == patchCorners.size());
        assert(chartSides.size() == patchEigenSides.size());

#ifdef SAVEMESHESFORDEBUG
        igl::writeOBJ(std::string("res/") + std::to_string(cId) + std::string("_chart.obj"), chartV, chartF);
#endif

        //Compute quadrangulation
        Eigen::MatrixXd uvMapV;
        Eigen::MatrixXi uvMapF;
        Eigen::MatrixXd quadrangulationV;
        Eigen::MatrixXi quadrangulationF;
        qfp::computeQuadrangulation(chartV, chartF, patchV, patchF, chartSideVertices, chartSideLength, chartSideSubdivision, patchEigenSides, uvMapV, uvMapF, quadrangulationV, quadrangulationF);

#ifdef SAVEMESHESFORDEBUG
        Eigen::MatrixXd uvMesh(uvMapV.rows(), 3);
        for (int i = 0; i < uvMapV.rows(); i++) {
            uvMesh(i, 0) = uvMapV(i, 0);
            uvMesh(i, 1) = uvMapV(i, 1);
            uvMesh(i, 2) = 0;
        }

        std::string uvFile = std::string("res/") + std::to_string(cId) + std::string("_uv.obj");
        igl::writeOBJ(uvFile, uvMesh, uvMapF);
#endif
        assert(chartV.rows() == uvMapV.rows());

        //Get polymesh
        PolyMesh quadrangulatedChartMesh;
        eigenToVCG(quadrangulationV, quadrangulationF, quadrangulatedChartMesh, 4);

#ifdef SAVEMESHESFORDEBUG
        igl::writeOBJ(std::string("res/") + std::to_string(cId) + std::string("_quadrangulation.obj"), quadrangulationV, quadrangulationF);
#endif

        //Smoothing
        if (chartSmoothingIterations > 0) {
            vcg::tri::UpdateSelection<PolyMesh>::VertexAll(quadrangulatedChartMesh);
            for (size_t vId : patchBorders) {
                quadrangulatedChartMesh.vert[vId].ClearS();
            }
            vcg::PolygonalAlgorithm<PolyMesh>::LaplacianReproject(quadrangulatedChartMesh, chartSmoothingIterations, 0.5, true);
        }

        std::vector<int> currentVertexMap(quadrangulatedChartMesh.vert.size(), -1);

        //Map subsides on the vertices of the current mesh (create if necessary)
        for (size_t i = 0; i < chartSides.size(); i++) {
            const ChartSide& side = chartSides[i];
            const std::vector<size_t>& patchSideVertices = patchEigenSides[i];

            size_t currentPatchSideVertexId = 0;

            for (size_t j = 0; j < side.subsides.size(); j++) {
                const size_t& subSideId = side.subsides[j];
                const bool& subsideReversed = side.reversedSubside[j];
                const ChartSubSide& subside = chartData.subSides[subSideId];

                //Create new vertices of the subsides
                if (vertexSubsideMap[subSideId].empty()) {
                    std::vector<size_t> subsideOrderedVertices = subside.vertices;
                    if (subsideReversed) {
                        std::reverse(subsideOrderedVertices.begin(), subsideOrderedVertices.end());
                    }

                    //Get fixed corners of the subside
                    size_t vStart = subsideOrderedVertices[0];
                    size_t vEnd = subsideOrderedVertices[subsideOrderedVertices.size() - 1];

                    size_t patchVStart = currentPatchSideVertexId;
                    size_t patchVEnd = (currentPatchSideVertexId + ilpResult[subSideId]) % patchSideVertices.size();

                    assert(patchVStart < patchSideVertices.size());
                    assert(patchVEnd < patchSideVertices.size());
                    assert(cornerVertices[vStart] >= 0 && cornerVertices[vEnd] >= 0);
                    assert(currentVertexMap[patchSideVertices[patchVStart]] == cornerVertices[vStart] || currentVertexMap[patchSideVertices[patchVStart]] == -1);
                    assert(currentVertexMap[patchSideVertices[patchVEnd]] == cornerVertices[vEnd] || currentVertexMap[patchSideVertices[patchVEnd]] == -1);

                    currentVertexMap[patchSideVertices[patchVStart]] = cornerVertices[vStart];
                    currentVertexMap[patchSideVertices[patchVEnd]] = cornerVertices[vEnd];

                    for (int k = 0; k <= ilpResult[subSideId]; k++) {
                        assert(currentPatchSideVertexId < patchSideVertices.size());
                        size_t patchSideVId = patchSideVertices[currentPatchSideVertexId];

                        if (currentVertexMap[patchSideVId] == -1) {
                            assert(k > 0 && k < ilpResult[subSideId]);

                            //Add new vertex
                            size_t newVertexId = quadmesh.vert.size();
                            const typename PolyMesh::CoordType& coord = quadrangulatedChartMesh.vert[patchSideVId].P();
                            vcg::tri::Allocator<PolyMesh>::AddVertex(quadmesh, coord);

                            vertexSubsideMap[subSideId].push_back(newVertexId);

                            fixedVertices.push_back(newVertexId);
                            currentVertexMap[patchSideVId] = newVertexId;
                        }
                        else {
                            //Use the existing vertex
                            int existingVertexId = currentVertexMap[patchSideVId];
                            assert(existingVertexId >= 0);
                            vertexSubsideMap[subSideId].push_back(existingVertexId);
                        }

                        currentPatchSideVertexId++;
                    }
                }
                //Set the existing vertices
                else {
                    assert(vertexSubsideMap[subSideId].size() == ilpResult[subSideId] + 1);

                    for (int k = 0; k <= ilpResult[subSideId]; k++) {                        
                        assert(currentPatchSideVertexId < patchSideVertices.size());
                        int patchSideVId = patchSideVertices[currentPatchSideVertexId];

                        size_t subSideVertexIndex = subsideReversed ? ilpResult[subSideId] - k : k;

                        currentVertexMap[patchSideVId] = vertexSubsideMap[subSideId][subSideVertexIndex];

                        size_t existingVertexId = currentVertexMap[patchSideVId];

                        //If it is not a corner
                        if (!subside.isFixed && k > 0 && k < ilpResult[subSideId]) {
                            //Average
                            const typename PolyMesh::CoordType& coord = quadrangulatedChartMesh.vert[patchSideVId].P();
                            quadmesh.vert.at(existingVertexId).P() = (coord + quadmesh.vert.at(existingVertexId).P()) / 2.0;
                        }

                        currentPatchSideVertexId++;
                    }
                }

                currentPatchSideVertexId--;
            }

            assert(currentPatchSideVertexId+1 == patchSideVertices.size());
        }

        //Internal vertices
        for (size_t i = 0; i < quadrangulatedChartMesh.vert.size(); i++) {
            if (currentVertexMap[i] == -1) {
                size_t newId = quadmesh.vert.size();

                const typename PolyMesh::CoordType& coord = quadrangulatedChartMesh.vert[i].P();
                vcg::tri::Allocator<PolyMesh>::AddVertex(quadmesh, coord);

                currentVertexMap[i] = newId;
            }
        }

        //Create faces and fill partitions
        for (size_t i = 0; i < quadrangulatedChartMesh.face.size(); i++) {
            assert(quadrangulatedChartMesh.face[i].VN() == 4);

            size_t newFaceId = quadmesh.face.size();

            vcg::tri::Allocator<PolyMesh>::AddFaces(quadmesh, 1);

            quadmesh.face[newFaceId].Alloc(quadrangulatedChartMesh.face[i].VN());
            for (int j = 0; j < quadrangulatedChartMesh.face[i].VN(); j++) {
                int vId = currentVertexMap[vcg::tri::Index(quadrangulatedChartMesh, quadrangulatedChartMesh.face[i].V(j))];
                assert(vId >= 0);

                quadmesh.face[newFaceId].V(j) = &quadmesh.vert[vId];
            }

            quadmeshPartitions[chart.label].push_back(newFaceId);
        }

        //Fill corners vertices
        for (size_t i = 0; i < chartSides.size(); i++) {
            const ChartSide& side = chartSides[i];
            size_t vStart = side.vertices[0];

            quadmeshCorners[chart.label].push_back(cornerVertices.at(vStart));
        }

#ifdef SAVEMESHESFORDEBUG
        vcg::tri::io::ExporterOBJ<PolyMesh>::Save(quadmesh, ("res/" + std::to_string(cId) + "_tmp_quadrangulation.obj").c_str(), vcg::tri::io::Mask::IOM_NONE);
#endif
    }

#ifdef SAVEMESHESFORDEBUG
    vcg::tri::io::ExporterOBJ<PolyMesh>::Save(quadmesh, "res/total_quadrangulation.obj", vcg::tri::io::Mask::IOM_NONE);
#endif

    int numDuplicateVertices = vcg::tri::Clean<PolyMesh>::RemoveDuplicateVertex(quadmesh);
    if (numDuplicateVertices > 0) {
        std::cout << "Removed " << numDuplicateVertices << " duplicate vertices." << std::endl;
    }
    int numUnreferencedVertices = vcg::tri::Clean<PolyMesh>::RemoveUnreferencedVertex(quadmesh);
    if (numUnreferencedVertices > 0) {
        std::cout << "Removed " << numUnreferencedVertices << " unreferenced vertices." << std::endl;
    }

#ifdef SAVEMESHESFORDEBUG
    vcg::tri::io::ExporterOBJ<PolyMesh>::Save(quadmesh, "res/quadrangulation_cleaned.obj", vcg::tri::io::Mask::IOM_NONE);
#endif

    vcg::tri::UpdateTopology<PolyMesh>::FaceFace(quadmesh);
    vcg::PolygonalAlgorithm<PolyMesh>::UpdateFaceNormalByFitting(quadmesh);
    OrientFaces<PolyMesh>::AutoOrientFaces(quadmesh);
    vcg::PolygonalAlgorithm<PolyMesh>::UpdateFaceNormalByFitting(quadmesh);
    vcg::tri::UpdateNormal<PolyMesh>::PerVertexNormalized(quadmesh);

    vcg::tri::UpdateNormal<TriangleMesh>::PerFaceNormalized(trimesh);
    vcg::tri::UpdateNormal<TriangleMesh>::PerVertexNormalized(trimesh);
    vcg::tri::UpdateBounding<TriangleMesh>::Box(trimesh);

    vcg::GridStaticPtr<typename TriangleMesh::FaceType,typename TriangleMesh::FaceType::ScalarType> Grid;
    Grid.Set(trimesh.face.begin(),trimesh.face.end());

    //Reproject
    vcg::tri::UpdateBounding<PolyMesh>::Box(quadmesh);
    typename TriangleMesh::ScalarType maxD=quadmesh.bbox.Diag();
    typename TriangleMesh::ScalarType minD=0;

    for (size_t i=0;i<quadmesh.vert.size();i++)
    {
        typename TriangleMesh::CoordType closestPT;
        typename TriangleMesh::FaceType *f=
                vcg::tri::GetClosestFaceBase<TriangleMesh>(
                    trimesh,
                    Grid,
                    quadmesh.vert[i].P(),
                    maxD,minD,
                    closestPT);

        quadmesh.vert[i].P()=closestPT;
    }

    if (quadrangulationSmoothingIterations > 0) {
        vcg::tri::UpdateSelection<PolyMesh>::VertexAll(quadmesh);
        for (const size_t& borderVertexId : fixedVertices) {
            quadmesh.vert[borderVertexId].ClearS();
        }
        vcg::PolygonalAlgorithm<PolyMesh>::template LaplacianReproject<TriangleMesh>(quadmesh, trimesh, quadrangulationSmoothingIterations, 0.7, 0.7, true);
    }

    vcg::PolygonalAlgorithm<PolyMesh>::UpdateFaceNormalByFitting(quadmesh);
    vcg::tri::UpdateNormal<PolyMesh>::PerVertexNormalized(quadmesh);
    vcg::tri::UpdateBounding<PolyMesh>::Box(quadmesh);

#ifdef SAVEMESHESFORDEBUG
    vcg::tri::io::ExporterOBJ<PolyMesh>::Save(quadmesh, "res/quadrangulation_reprojected.obj", vcg::tri::io::Mask::IOM_NONE);
#endif
}

}
