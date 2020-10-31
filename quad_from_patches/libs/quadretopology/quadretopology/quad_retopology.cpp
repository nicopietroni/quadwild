#include "quad_retopology.h"

#include "includes/qr_convert.h"
#include "includes/qr_field_tracer.h"
#include "includes/qr_patch_tracer.h"
#include "includes/qr_utils.h"
#include "includes/qr_patterns.h"
#include "includes/qr_mapping.h"
#include "includes/qr_patch_assembler.h"
#include <map>

#include <vcg/complex/algorithms/polygonal_algorithms.h>

#ifdef SAVE_MESHES_FOR_DEBUG
#include <igl/writeOBJ.h>
#endif

namespace QuadRetopology {


template<class TriangleMeshType, class PolyMeshType>
std::vector<int> patchDecomposition(
        TriangleMeshType& newSurface,
        PolyMeshType& preservedSurface,
        std::vector<std::vector<size_t>>& partitions,
        std::vector<std::vector<size_t>>& corners,
        const Parameters& parameters)
{
    return QuadRetopology::patchDecomposition(
        newSurface, 
        preservedSurface, 
        partitions, 
        corners,  
        parameters.initialRemeshing,
        parameters.initialRemeshingEdgeFactor,
        parameters.reproject,
        parameters.splitConcaves,
        parameters.finalSmoothing);
            
}
template<class TriangleMeshType, class PolyMeshType>
std::vector<int> patchDecomposition(
        TriangleMeshType& newSurface,
        PolyMeshType& preservedSurface,
        std::vector<std::vector<size_t>>& partitions,
        std::vector<std::vector<size_t>>& corners,
        const bool initialRemeshing,
        const double edgeFactor,
        const bool reproject,
        const bool splitConcaves,
        const bool finalSmoothing)
{
    if (newSurface.face.size() <= 0)
        return std::vector<int>();

    std::vector<std::vector<std::vector<std::pair<size_t,size_t>>>> sides;
    internal::PatchAssembler<TriangleMeshType, PolyMeshType> patchAssembler(newSurface, preservedSurface);
    typename internal::PatchAssembler<TriangleMeshType, PolyMeshType>::Parameters parameters;
    parameters.InitialRemesh = initialRemeshing;
    parameters.EdgeSizeFactor = edgeFactor;
    parameters.FinalSmooth = finalSmoothing;
    parameters.SplitAllConcave = splitConcaves;
    parameters.Reproject = reproject;
    patchAssembler.BatchProcess(partitions, corners, sides, parameters);

    std::vector<int> newSurfaceLabel(newSurface.face.size(), -1);
    for (size_t pId = 0; pId < partitions.size(); pId++) {
        for (const size_t& fId : partitions[pId]) {
            assert(newSurfaceLabel[fId] == -1);
            newSurfaceLabel[fId] = static_cast<int>(pId);
        }
    }

    return newSurfaceLabel;
}

template<class TriangleMesh>
ChartData computeChartData(
        TriangleMesh& mesh,
        const std::vector<std::vector<size_t>>& meshPartitions,
        const std::vector<std::vector<size_t>>& meshCorners)
{
    std::vector<int> faceLabel(mesh.face.size(), -1);
    for (size_t pId = 0; pId < meshPartitions.size(); pId++) {
        for (const size_t& fId : meshPartitions[pId]) {
            assert(faceLabel[fId] == -1);
            faceLabel[fId] = static_cast<int>(pId);
        }
    }

    ChartData chartData = computeChartData(mesh, faceLabel, meshCorners);
    return chartData;
}

//It works just on triangle meshes
template<class TriangleMeshType>
ChartData computeChartData(
        TriangleMeshType& mesh,
        const std::vector<int>& faceLabel,
        const std::vector<std::vector<size_t>>& corners)
{
    typedef std::map<std::pair<size_t, size_t>, std::pair<int, int>> EdgeLabelMap;
    typedef std::map<std::pair<size_t, size_t>, int> EdgeSubSideMap;

    ChartData chartData;

    if (mesh.face.size() == 0)
        return chartData;

    vcg::tri::UpdateTopology<TriangleMeshType>::FaceFace(mesh);

    //Region growing algorithm for getting charts
    internal::findChartFacesAndBorderFaces(mesh, faceLabel, chartData);

    //TODO SPLIT IN FUNCTIONS
    EdgeSubSideMap edgeSubSideMap;
    for (const int& pId : chartData.labels) {
        Chart& chart = chartData.charts[pId];

        assert(chart.label == pId);

        if (chart.faces.size() == 0)
            continue;

        std::unordered_set<size_t> cornerSet(corners[pId].begin(), corners[pId].end());

#ifndef NDEBUG
        if (cornerSet.size() < 3 || cornerSet.size() > 6) {
            std::cout << "Warning 3: Given as input for " << pId << ": " << chart.chartSides.size() << " sides." << std::endl;
        }
#endif

        EdgeLabelMap edgeLabelMap;
        std::vector<std::vector<size_t>> vertexNextMap(mesh.vert.size());

        std::set<size_t> remainingVertices;

        //Fill edge map and next vertex map
        for (const size_t& fId : chart.borderFaces) {
            typename TriangleMeshType::FaceType* currentFacePointer = &mesh.face[fId];
            vcg::face::Pos<typename TriangleMeshType::FaceType> pos(currentFacePointer, 0);

            for (int k = 0; k < currentFacePointer->VN(); k++) {
                pos.FlipF();
                size_t adjFace = vcg::tri::Index(mesh, pos.F());
                int adjLabel = faceLabel[adjFace];

                bool isBorderEdge = false;
                int adjChartLabel = -2;

                if (currentFacePointer == pos.F()) {
                    adjChartLabel = -1;
                    isBorderEdge = true;
                }
                else if (adjLabel != chart.label) {
                    adjChartLabel = adjLabel;
                    isBorderEdge = true;
                }
                pos.FlipF();

                //For each border edge
                if (isBorderEdge) {
                    assert(adjChartLabel > -2);

                    typename TriangleMeshType::VertexType* vStart = pos.V();
                    pos.FlipV();
                    typename TriangleMeshType::VertexType* vEnd = pos.V();
                    pos.FlipV();

                    size_t vStartId = vcg::tri::Index(mesh, vStart);
                    size_t vEndId = vcg::tri::Index(mesh, vEnd);

                    std::pair<size_t, size_t> edge(vStartId, vEndId);
                    if (edge.first > edge.second) {
                        std::swap(edge.first, edge.second);
                    }

                    edgeLabelMap.insert(std::make_pair(edge, std::make_pair(chart.label, adjChartLabel)));
                    vertexNextMap[vStartId].push_back(vEndId);

                    remainingVertices.insert(vStartId);
                    remainingVertices.insert(vEndId);
                }

                pos.FlipV();
                pos.FlipE();
            }
        }

        do {
            //Find first label
            size_t vStartId;
            size_t vCurrentId;
            size_t vNextId;

            //Corner detection variables
            typename TriangleMeshType::CoordType lastEdgeVec;
            bool isCorner = false;

            vCurrentId = *remainingVertices.begin();

            std::vector<size_t> nextConfiguration = internal::findVertexChainPath(vCurrentId, vertexNextMap);
            vNextId = vertexNextMap[vCurrentId][nextConfiguration[vCurrentId]];

            //Get last edge vector
            lastEdgeVec = mesh.vert[vNextId].P() - mesh.vert[vCurrentId].P();
            lastEdgeVec.Normalize();

//            std::pair<size_t, size_t> startEdge(vCurrentId, vNextId);
//            if (startEdge.first > startEdge.second) {
//                std::swap(startEdge.first, startEdge.second);
//            }

            int currentLabel;

            //Iterate in the borders to get the first corner
            vStartId = vCurrentId;
            size_t firstCornerIterations = 0;
            do {
                //Next border edge
                vCurrentId = vertexNextMap[vCurrentId][nextConfiguration[vCurrentId]];
                vNextId = vertexNextMap[vCurrentId][nextConfiguration[vCurrentId]];

                typename TriangleMeshType::CoordType currentEdgeVec = mesh.vert[vNextId].P() - mesh.vert[vCurrentId].P();
                currentEdgeVec.Normalize();

                //Check if it is a corner
                isCorner = cornerSet.find(vCurrentId) != cornerSet.end();

                lastEdgeVec = currentEdgeVec;

                firstCornerIterations++;
            } while (!isCorner && vCurrentId != vStartId && firstCornerIterations < MAXITERATIONS);

#ifndef NDEBUG
            if (firstCornerIterations >= MAXITERATIONS) {
                std::cout << "Error: error iterating! Cannot find the first corner or get back to the start vertex." << std::endl;
            }
#endif

#ifndef NDEBUG
            if (vCurrentId == vStartId) {
                std::cout << "Warning 1: input mesh is not well-defined: no corners!" << std::endl;
            }
#endif
            vStartId = vCurrentId;
            vNextId = vertexNextMap[vCurrentId][nextConfiguration[vCurrentId]];

            ChartSide currentSide;
            currentSide.length = 0;
            currentSide.size = 0;


            int adjChartLabel;
            do {
                size_t subSideId = chartData.subsides.size();
                ChartSubside currentSubSide;

                //Get edge
                std::pair<size_t, size_t> edge(vCurrentId, vNextId);
                if (edge.first > edge.second) {
                    std::swap(edge.first, edge.second);
                }
                //Get current label on the other side
                const std::pair<int,int>& currentEdgeLabels = edgeLabelMap.at(edge);
                assert(currentEdgeLabels.first == chart.label || currentEdgeLabels.second == chart.label);
                adjChartLabel = currentEdgeLabels.first == chart.label ? currentEdgeLabels.second : currentEdgeLabels.first;

                std::unordered_set<size_t> cornerSetAdj;
                if (adjChartLabel >= 0)
                    cornerSetAdj.insert(corners[adjChartLabel].begin(), corners[adjChartLabel].end());

                double length = 0;

                bool newSubSide = false;

                bool firstIteration = true;
                isCorner = false;
                bool isAdjCorner = false;
                size_t vSubSideStartId = vCurrentId;
                size_t iterations = 0;
                do {
                    typename TriangleMeshType::CoordType currentEdgeVec = mesh.vert[vNextId].P() - mesh.vert[vCurrentId].P();
                    currentEdgeVec.Normalize();

                    std::pair<size_t, size_t> edge(vCurrentId, vNextId);
                    if (edge.first > edge.second) {
                        std::swap(edge.first, edge.second);
                    }

                    //Check if it is a corner
                    if (!firstIteration) {
                        isCorner = cornerSet.find(vCurrentId) != cornerSet.end();
                        isAdjCorner = cornerSetAdj.find(vCurrentId) != cornerSetAdj.end();
                    }

                    //Get current label on the other subside
                    const std::pair<int,int>& currentEdgeLabels = edgeLabelMap.at(edge);
                    assert(currentEdgeLabels.first == chart.label || currentEdgeLabels.second == chart.label);
                    currentLabel = currentEdgeLabels.first == chart.label ? currentEdgeLabels.second : currentEdgeLabels.first;

                    if (!isCorner && !isAdjCorner && currentLabel == adjChartLabel) {
                        EdgeSubSideMap::iterator findIt = edgeSubSideMap.find(edge);

                        //If the subside has already been processed
                        if (findIt == edgeSubSideMap.end()) {
                            currentSubSide.vertices.push_back(vCurrentId);

                            length += (mesh.vert[vNextId].P() - mesh.vert[vCurrentId].P()).Norm();

                            edgeSubSideMap.insert(std::make_pair(edge, subSideId));

                            newSubSide = true;
                        }
                        else if (firstIteration) {
                            subSideId = findIt->second;
                        }
                        firstIteration = false;

                        remainingVertices.erase(vCurrentId);

                        //Next border edge
                        vCurrentId = vertexNextMap[vCurrentId][nextConfiguration[vCurrentId]];
                        vNextId = vertexNextMap[vCurrentId][nextConfiguration[vCurrentId]];

                        lastEdgeVec = currentEdgeVec;
                    }
                } while (!isCorner && !isAdjCorner && currentLabel == adjChartLabel && vCurrentId != vSubSideStartId && iterations < MAXITERATIONS);
#ifndef NDEBUG
                if (iterations >= MAXITERATIONS) {
                    std::cout << "Error: error iterating! Cannot find a corner or get back to the start vertex." << std::endl;
                }
#endif
#ifndef NDEBUG
                if (vCurrentId == vSubSideStartId) {
                    std::cout << "Warning 2: input mesh is not well-defined: single border chart with no corners!" << std::endl;
                }
#endif

                //True if the subside is reversed (from the last to the first vertex)
                bool reversed;

                if (newSubSide) {
                    //Add last vertex
                    currentSubSide.vertices.push_back(vCurrentId);

                    //Create new side
                    int ChartSubsideId = chart.chartSubsides.size();
                    chart.chartSubsides.push_back(subSideId);

                    currentSubSide.incidentCharts[0] = chart.label;
                    currentSubSide.incidentCharts[1] = adjChartLabel;

                    currentSubSide.length = length;

                    currentSubSide.incidentChartSideId[0] = ChartSubsideId;

                    if (adjChartLabel >= 0) {
                        currentSubSide.isOnBorder = false;

                        chart.adjacentCharts.push_back(adjChartLabel);
                    }
                    else {
                        currentSubSide.isOnBorder = true;
                        currentSubSide.incidentChartSideId[1] = -1;
                    }

                    assert(currentSubSide.vertices.size() >= 2);
                    currentSubSide.size = currentSubSide.vertices.size() - 1;

                    chartData.subsides.push_back(currentSubSide);


                    //Pop last vertex
                    if (currentSide.vertices.size() > 0) {
                        assert(currentSide.vertices.back() == currentSubSide.vertices.front());
                        currentSide.vertices.pop_back();
                    }
                    currentSide.vertices.insert(
                                currentSide.vertices.end(),
                                currentSubSide.vertices.begin(),
                                currentSubSide.vertices.end());

                    reversed = false;
                }
                else {
                    assert(currentSubSide.vertices.size() == 0);

                    //Add the side to other chart
                    if (adjChartLabel >= 0) {
                        int chartSideId = chart.chartSubsides.size();
                        chart.chartSubsides.push_back(subSideId);

                        assert(chartData.subsides[subSideId].incidentCharts[1] == chart.label);
                        chartData.subsides[subSideId].incidentChartSideId[1] = chartSideId;

                        chart.adjacentCharts.push_back(adjChartLabel);

                        //Pop last vertex
                        if (currentSide.vertices.size() > 0) {
                            assert(currentSide.vertices.back() == chartData.subsides[subSideId].vertices.back());
                            currentSide.vertices.pop_back();
                        }

                        //Add side vertices
                        currentSide.vertices.insert(
                                    currentSide.vertices.end(),
                                    chartData.subsides[subSideId].vertices.rbegin(),
                                    chartData.subsides[subSideId].vertices.rend());
                    }

                    reversed = true;
                }

                currentSide.subsides.push_back(subSideId);
                currentSide.reversedSubside.push_back(reversed);

                currentSide.length += chartData.subsides[subSideId].length;
                currentSide.size += chartData.subsides[subSideId].size;

                if (isCorner) {
                    chart.chartSides.push_back(currentSide);
                    currentSide = ChartSide();
                }

            } while (vCurrentId != vStartId);

#ifndef NDEBUG
            if (!isCorner) {
                std::cout << "Warning 4: Chart has no final corner!" << std::endl;
            }
#endif

        } while (!remainingVertices.empty());

#ifndef NDEBUG
        if (chart.chartSides.size() < 3 || chart.chartSides.size() > 6) {
            std::cout << "Warning 3: Chart " << pId << " has " << chart.chartSides.size() << " sides." << std::endl;
        }
#endif
    }

    return chartData;
}

inline std::vector<double> computeChartEdgeLength(
        const ChartData& chartData,
        const std::vector<size_t> fixedSubsides,
        const size_t& iterations,
        const double& weight)
{
    vector<bool> isFixed(chartData.subsides.size(), false);
    for (int sId : fixedSubsides) {
        isFixed[sId] = true;
    }

    std::vector<double> avgLengths(chartData.charts.size() , -1);

    //Fill charts with a border
    for (size_t i = 0; i < chartData.charts.size(); i++) {
        const Chart& chart = chartData.charts[i];
        if (chart.faces.size() > 0) {
            double currentQuadLength = 0;
            int numSides = 0;

            for (size_t sId : chart.chartSubsides) {
                const ChartSubside& subside = chartData.subsides[sId];
                if (isFixed[sId]) {
                    currentQuadLength += subside.length / subside.size;
                    numSides++;
                }
            }

            if (numSides > 0) {
                currentQuadLength /= numSides;
                avgLengths[i] = currentQuadLength;
            }
        }
    }

    //Fill charts with no borders
    bool done;
    do {
        done = true;
        for (size_t i = 0; i < chartData.charts.size(); i++) {
            const Chart& chart = chartData.charts[i];
            if (chart.faces.size() > 0 && avgLengths[i] < 0) {
                double currentLength = 0;
                size_t numAdjacentCharts = 0;

                for (size_t adjId : chart.adjacentCharts) {
                    if (avgLengths[adjId] > 0) {
                        currentLength += avgLengths[adjId];
                        numAdjacentCharts++;
                        done = false;
                    }
                }

                if (currentLength > 0) {
                    currentLength /= numAdjacentCharts;
                    avgLengths[i] = currentLength;
                }
            }
        }
    } while (!done);


    //Smoothing
    for (size_t k = 0; k < iterations; k++) {
        std::vector<double> lastAvgLengths = avgLengths;

        for (size_t i = 0; i < chartData.charts.size(); i++) {
            const Chart& chart = chartData.charts[i];
            if (chart.faces.size() > 0) {
                assert(lastAvgLengths[i] > 0);

                double adjValue = 0.0;
                size_t numAdjacentCharts = 0;

                for (size_t adjId : chart.adjacentCharts) {
                    if (avgLengths[adjId] > 0) {
                        adjValue += lastAvgLengths[adjId];
                        numAdjacentCharts++;
                    }
                }

                if (adjValue > 0.0) {
                    adjValue = weight * lastAvgLengths[i] + (1.0 - weight) * adjValue;
                }
            }
        }
    }

    return avgLengths;
}

inline std::vector<int> findSubdivisions(
        const ChartData& chartData,
        const std::vector<size_t> fixedSubsides,
        const std::vector<double>& chartEdgeLength,
        const Parameters& parameters,
        double& gap)
{
    return findSubdivisions(
        chartData,
        fixedSubsides,
        chartEdgeLength,
        parameters.ilpMethod,
        parameters.alpha,
        parameters.isometry,
        parameters.regularityForQuadrilaterals,
        parameters.regularityForNonQuadrilaterals,
        parameters.regularityNonQuadrilateralWeight,
        parameters.feasibilityFix,
        parameters.hardParityConstraint,
        parameters.timeLimit,
        parameters.gapLimit,
        parameters.minimumGap,
        gap);
}

inline std::vector<int> findSubdivisions(
        const ChartData& chartData,
        const std::vector<size_t> fixedSubsides,
        const std::vector<double>& chartEdgeLength,
        const ILPMethod& method,
        const double alpha,
        const bool isometry,
        const bool regularityForQuadrilaterals,
        const bool regularityForNonQuadrilaterals,
        const double regularityNonQuadrilateralWeight,
        const bool feasibilityFix,
        const bool hardParityConstraint,
        const double timeLimit,
        const double gapLimit,
        const double minimumGap,
        double& gap)
{
    if (chartData.charts.size() <= 0)
        return std::vector<int>();

    ILPStatus status;

    //Solve ILP to find the best patches
    std::vector<int> ilpResult = internal::solveILP(
        chartData,
        fixedSubsides,
        chartEdgeLength,
        method,
        alpha,
        isometry,
        regularityForQuadrilaterals,
        regularityForNonQuadrilaterals,
        regularityNonQuadrilateralWeight,
        feasibilityFix,
        hardParityConstraint,
        timeLimit,
        gapLimit,
        gap,
        status);

    if (status == ILPStatus::SOLUTIONFOUND && gap < minimumGap) {
        std::cout << "Solution found! Gap: " << gap << std::endl;
    }
    else {
        if (status == ILPStatus::INFEASIBLE) {
            if (method == ILPMethod::LEASTSQUARES) {
                std::cout << "Model was infeasible or time limit exceeded. Trying with ABS to reduce time." << gap << std::endl;

                return findSubdivisions(
                    chartData,
                    fixedSubsides,
                    chartEdgeLength,
                    ILPMethod::ABS,
                    alpha,
                    isometry,
                    regularityForQuadrilaterals,
                    regularityForNonQuadrilaterals,
                    regularityNonQuadrilateralWeight,
                    feasibilityFix,
                    hardParityConstraint,
                    timeLimit,
                    gapLimit,
                    minimumGap,
                    gap);
            }
            else {
                std::cout << "Error! Model was infeasible or time limit exceeded!" << std::endl;
            }
        }
        else if (status == ILPStatus::SOLUTIONWRONG && !hardParityConstraint) {
            std::cout << "Solution wrong! It have been used soft constraints for parity, so trying with hard constraints." << std::endl;

            return findSubdivisions(
                chartData,
                fixedSubsides,
                chartEdgeLength,
                method,
                alpha,
                isometry,
                regularityForQuadrilaterals,
                regularityForNonQuadrilaterals,
                regularityNonQuadrilateralWeight,
                feasibilityFix,
                true,
                timeLimit,
                gapLimit,
                minimumGap,
                gap);
        }
        else {
            if (method == ILPMethod::LEASTSQUARES) {
                std::cout << "Minimum gap has been not reached. Trying with ABS (linear optimization method)." << gap << std::endl;

                return findSubdivisions(
                    chartData,
                    fixedSubsides,
                    chartEdgeLength,
                    ILPMethod::ABS,
                    alpha,
                    isometry,
                    regularityForQuadrilaterals,
                    regularityForNonQuadrilaterals,
                    regularityNonQuadrilateralWeight,
                    feasibilityFix,
                    true,
                    timeLimit,
                    gapLimit,
                    minimumGap,
                    gap);
            }
            else if (regularityForNonQuadrilaterals) {
                std::cout << "Minimum gap has been not reached. Trying without regularity for non-quadrilaterals." << gap << std::endl;

                return findSubdivisions(
                    chartData,
                    fixedSubsides,
                    chartEdgeLength,
                    method,
                    alpha,
                    isometry,
                    regularityForQuadrilaterals,
                    false,
                    regularityNonQuadrilateralWeight,
                    feasibilityFix,
                    true,
                    timeLimit,
                    gapLimit,
                    minimumGap,
                    gap);
            }
            else if (regularityForQuadrilaterals) {
                std::cout << "Minimum gap has been not reached. Trying without any regularity terms." << gap << std::endl;

                return findSubdivisions(
                    chartData,
                    fixedSubsides,
                    chartEdgeLength,
                    method,
                    alpha,
                    true,
                    false,
                    false,
                    regularityNonQuadrilateralWeight,
                    feasibilityFix,
                    true,
                    timeLimit,
                    gapLimit,
                    minimumGap,
                    gap);
            }
        }
    }

    return ilpResult;
}

template<class TriangleMeshType, class PolyMeshType>
void quadrangulate(
        TriangleMeshType& newSurface,
        const ChartData& chartData,
        const std::vector<size_t> fixedSubsides,
        const std::vector<int>& ilpResult,
        const Parameters& parameters,
        PolyMeshType& quadrangulation,
        std::vector<int>& quadrangulationFaceLabel,
        std::vector<std::vector<size_t>>& quadrangulationPartitions,
        std::vector<std::vector<size_t>>& quadrangulationCorners)
{
    return QuadRetopology::quadrangulate(
            newSurface,
            chartData,
            fixedSubsides,
            ilpResult,
            parameters.chartSmoothingIterations,
            parameters.quadrangulationFixedSmoothingIterations,
            parameters.quadrangulationNonFixedSmoothingIterations,
            parameters.doubletRemoval,
            quadrangulation,
            quadrangulationFaceLabel,
            quadrangulationPartitions,
            quadrangulationCorners);
}

template<class TriangleMeshType, class PolyMeshType>
void quadrangulate(
        TriangleMeshType& newSurface,
        const ChartData& chartData,
        const std::vector<size_t> fixedSubsides,
        const std::vector<int>& ilpResult,
        const int chartSmoothingIterations,
        const int quadrangulationFixedSmoothingIterations,
        const int quadrangulationNonFixedSmoothingIterations,
        const bool doubletRemoval,
        PolyMeshType& quadrangulation,
        std::vector<int>& quadrangulationFaceLabel,
        std::vector<std::vector<size_t>>& quadrangulationPartitions,
        std::vector<std::vector<size_t>>& quadrangulationCorners)
{
    if (newSurface.face.size() <= 0)
        return;
    if (ilpResult.size() == 0)
        return;

    std::vector<std::vector<size_t>> subsideVertexMap(chartData.subsides.size());
    std::vector<int> cornerVertices(newSurface.vert.size(), -1);

    quadrangulationPartitions.resize(chartData.charts.size());
    quadrangulationCorners.resize(chartData.charts.size());

    vector<bool> isFixed(chartData.subsides.size(), false);
    for (int sId : fixedSubsides) {
        isFixed[sId] = true;
    }

    //Fill fixed vertices (subsides corners)
    for (const ChartSubside& subside : chartData.subsides) {
        size_t vStart = subside.vertices[0];
        size_t vEnd = subside.vertices[subside.vertices.size() - 1];

        if (cornerVertices[vStart] == -1) {
            cornerVertices[vStart] = quadrangulation.vert.size();
            vcg::tri::Allocator<PolyMeshType>::AddVertex(
                        quadrangulation,
                        newSurface.vert[vStart].P());
        }

        if (cornerVertices[vEnd] == -1) {
            cornerVertices[vEnd] = quadrangulation.vert.size();
            vcg::tri::Allocator<PolyMeshType>::AddVertex(
                        quadrangulation,
                        newSurface.vert[vEnd].P());
        }
    }

    //Fill subside map for fixed borders
    std::set<size_t> fixedVerticesSet;
    for (size_t subsideId = 0; subsideId < chartData.subsides.size(); subsideId++) {
        const ChartSubside& subside = chartData.subsides[subsideId];
        if (isFixed[subsideId]) {
            for (size_t k = 0; k < subside.vertices.size(); k++) {
                const size_t& vId = subside.vertices[k];

                size_t newVertexId;

                if (cornerVertices[vId] == -1) {
                    assert(k > 0 && k < subside.vertices.size() - 1);

                    newVertexId = quadrangulation.vert.size();
                    vcg::tri::Allocator<PolyMeshType>::AddVertex(
                                quadrangulation,
                                newSurface.vert[vId].P());
                }
                else {
                    newVertexId = cornerVertices[vId];
                    assert(newVertexId >= 0);
                }

                fixedVerticesSet.insert(newVertexId);
                subsideVertexMap[subsideId].push_back(newVertexId);
            }

            if (ilpResult[subsideId] > subside.size) {
                int vToSplit = -1;
                double maxLength = 0.0;
                for (size_t k = 0; k < subside.vertices.size() - 1; k++) {
                    const size_t& vId1 = subside.vertices[k];
                    const size_t& vId2 = subside.vertices[k + 1];
                    double length = (newSurface.vert[vId2].P() - newSurface.vert[vId1].P()).Norm();
                    if (length >= maxLength) {
                        vToSplit = k;
                    }
                }

                if (vToSplit >= 0) {
                    const size_t& vId1 = subside.vertices[vToSplit];
                    const size_t& vId2 = subside.vertices[vToSplit + 1];
                    size_t splitVertexId = quadrangulation.vert.size();
                    vcg::tri::Allocator<PolyMeshType>::AddVertex(
                                quadrangulation,
                                (newSurface.vert[vId1].P() + newSurface.vert[vId2].P()) / 2.0);
                    vcg::tri::Allocator<PolyMeshType>::AddFace(
                                quadrangulation, subsideVertexMap[subsideId][vToSplit], subsideVertexMap[subsideId][vToSplit + 1], splitVertexId);

                    subsideVertexMap[subsideId].insert(subsideVertexMap[subsideId].begin() + vToSplit + 1, splitVertexId);

                    std::cout << "Triangle added in subside " << subsideId << ": +1!" << std::endl;
                }
                else {
                    std::cout << "ERROR: impossible to augment the subside " << subsideId << ": +1! Target vertex not found." << std::endl;
                }
            }
            else if (ilpResult[subsideId] < subside.size) {
                if (subside.size >= 2) {
                    int vToSkip = -1;
                    double minLength = std::numeric_limits<double>::max();
                    for (size_t k = 0; k < subside.vertices.size() - 2; k++) {
                        const size_t& vId1 = subside.vertices[k];
                        const size_t& vId2 = subside.vertices[k + 1];
                        const size_t& vId3 = subside.vertices[k + 2];
                        double length = (newSurface.vert[vId2].P() - newSurface.vert[vId1].P()).Norm() + (newSurface.vert[vId3].P() - newSurface.vert[vId2].P()).Norm();
                        if (length <= minLength) {
                            vToSkip = k;
                        }
                    }

                    if (vToSkip >= 0) {
                        vcg::tri::Allocator<PolyMeshType>::AddFace(
                                    quadrangulation, subsideVertexMap[subsideId][vToSkip], subsideVertexMap[subsideId][vToSkip + 1], subsideVertexMap[subsideId][vToSkip + 2]);

                        subsideVertexMap[subsideId].erase(subsideVertexMap[subsideId].begin() + vToSkip + 1);

                        std::cout << "Triangle added in subside " << subsideId << ": -1!" << std::endl;
                    }
                    else {
                        std::cout << "ERROR: impossible to reduce the subside " << subsideId << ": -1! Target vertex not found." << std::endl;
                    }
                }
                else {
                    std::cout << "ERROR: impossible to reduce the subside " << subsideId << ": -1! Subside is less than 2." << std::endl;
                }
            }
        }
    }


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
        for (size_t sId : chart.chartSubsides) {
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
        vcg::tri::UpdateFlags<TriangleMeshType>::FaceClearS(newSurface);
        vcg::tri::UpdateFlags<TriangleMeshType>::VertexClearS(newSurface);
        for (const size_t& fId : chart.faces) {
            newSurface.face[fId].SetS();
            for (int k = 0; k < newSurface.face[fId].VN(); k++) {
                newSurface.face[fId].V(k)->SetS();
            }
        }
        std::vector<int> vMap, fMap;
        QuadRetopology::internal::VCGToEigen(newSurface, chartV, chartF, vMap, fMap, true, 3);

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
                const ChartSubside& subSide = chartData.subsides[subSideId];

                if (ilpResult[subSideId] < 0) {
                    std::cout << "Error: ILP not valid" << std::endl;
                    return;
                }

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

            }

            l(static_cast<int>(i)) = targetSideSubdivision;
        }

        //Pattern quadrangulation
        Eigen::MatrixXd patchV;
        Eigen::MatrixXi patchF;
        std::vector<size_t> patchBorders;
        std::vector<size_t> patchCorners;
        PolyMeshType patchMesh;
        std::vector<std::vector<size_t>> patchSides;
        QuadRetopology::internal::computePattern(l, patchV, patchF, patchMesh, patchBorders, patchCorners, patchSides);

#ifdef SAVE_MESHES_FOR_DEBUG
        igl::writeOBJ(std::string("results/") + std::to_string(cId) + std::string("_patch.obj"), patchV, patchF);
#endif

        assert(chartSides.size() == patchCorners.size());
        assert(chartSides.size() == patchSides.size());

#ifdef SAVE_MESHES_FOR_DEBUG
        igl::writeOBJ(std::string("results/") + std::to_string(cId) + std::string("_chart.obj"), chartV, chartF);
#endif

        //Compute quadrangulation
        Eigen::MatrixXd uvMapV;
        Eigen::MatrixXi uvMapF;
        Eigen::MatrixXd quadrangulationV;
        Eigen::MatrixXi quadrangulationF;
        QuadRetopology::internal::computeQuadrangulation(chartV, chartF, patchV, patchF, chartSideVertices, chartSideLength, chartSideSubdivision, patchSides, uvMapV, uvMapF, quadrangulationV, quadrangulationF);

#ifdef SAVE_MESHES_FOR_DEBUG
        Eigen::MatrixXd uvMesh(uvMapV.rows(), 3);
        for (int i = 0; i < uvMapV.rows(); i++) {
            uvMesh(i, 0) = uvMapV(i, 0);
            uvMesh(i, 1) = uvMapV(i, 1);
            uvMesh(i, 2) = 0;
        }

        std::string uvFile = std::string("results/") + std::to_string(cId) + std::string("_uv.obj");
        igl::writeOBJ(uvFile, uvMesh, uvMapF);
#endif
        assert(chartV.rows() == uvMapV.rows());

        //Get polymesh
        PolyMeshType quadrangulatedChartMesh;
        QuadRetopology::internal::eigenToVCG(quadrangulationV, quadrangulationF, quadrangulatedChartMesh, 4);

#ifdef SAVE_MESHES_FOR_DEBUG
        igl::writeOBJ(std::string("results/") + std::to_string(cId) + std::string("_quadrangulation.obj"), quadrangulationV, quadrangulationF);
#endif

        //Smoothing
        if (chartSmoothingIterations > 0) {
            vcg::tri::UpdateSelection<PolyMeshType>::VertexAll(quadrangulatedChartMesh);
            for (size_t vId : patchBorders) {
                quadrangulatedChartMesh.vert[vId].ClearS();
            }
            vcg::PolygonalAlgorithm<PolyMeshType>::LaplacianReproject(quadrangulatedChartMesh, chartSmoothingIterations, 0.5, true);
        }

        std::vector<int> currentVertexMap(quadrangulatedChartMesh.vert.size(), -1);

        //Map subsides on the vertices of the current mesh (create if necessary)
        for (size_t i = 0; i < chartSides.size(); i++) {
            const ChartSide& side = chartSides[i];
            const std::vector<size_t>& patchSide = patchSides[i];

            size_t currentPatchSideVertex = 0;

            for (size_t j = 0; j < side.subsides.size(); j++) {
                const size_t& subsideId = side.subsides[j];
                const bool& reversed = side.reversedSubside[j];
                const ChartSubside& subside = chartData.subsides[subsideId];

                //Create new vertices of the subsides
                if (subsideVertexMap[subsideId].empty()) {
                    assert(!isFixed[subsideId]);

                    //Get fixed corners of the subside
                    size_t vStart = subside.vertices[0];
                    size_t vEnd = subside.vertices[subside.vertices.size() - 1];
                    assert(cornerVertices[vStart] >= 0 && cornerVertices[vEnd] >= 0);

                    currentVertexMap[patchSide[currentPatchSideVertex]] = cornerVertices[vStart];
                    currentVertexMap[patchSide[currentPatchSideVertex + ilpResult[subsideId]]] = cornerVertices[vEnd];

                    for (int k = 0; k <= ilpResult[subsideId]; k++) {
                        size_t patchSideVId = patchSide[currentPatchSideVertex];

                        if (currentVertexMap[patchSideVId] == -1) {
                            assert(k > 0 && k < ilpResult[subsideId]);

                            //Add new vertex
                            size_t newVertexId = quadrangulation.vert.size();

                            const typename PolyMeshType::CoordType& coord = quadrangulatedChartMesh.vert[patchSideVId].P();
                            vcg::tri::Allocator<PolyMeshType>::AddVertex(quadrangulation, coord);

                            currentVertexMap[patchSideVId] = newVertexId;

                            subsideVertexMap[subsideId].push_back(newVertexId);
                        }
                        else {
                            //Use the existing vertex
                            int existingVertexId = currentVertexMap[patchSideVId];
                            assert(existingVertexId >= 0);
                            subsideVertexMap[subsideId].push_back(existingVertexId);
                        }

                        currentPatchSideVertex++;
                    }

                    if (reversed) {
                        std::reverse(subsideVertexMap[subsideId].begin(), subsideVertexMap[subsideId].end());
                    }
                }
                //Set the existing vertices
                else {
                    assert(subsideVertexMap[subsideId].size() == ilpResult[subsideId] + 1);

                    for (int k = 0; k <= ilpResult[subsideId]; k++) {
                        int patchSideVId = patchSide[currentPatchSideVertex];

                        size_t subSideVertexIndex = reversed ? ilpResult[subsideId] - k : k;

                        currentVertexMap[patchSideVId] = subsideVertexMap[subsideId][subSideVertexIndex];

                        size_t existingVertexId = currentVertexMap[patchSideVId];

                        //If it is not a corner or if it is not on border
                        if (!isFixed[subsideId] && k > 0 && k < ilpResult[subsideId]) {
                            //Average
                            const typename PolyMeshType::CoordType& coord = quadrangulatedChartMesh.vert[patchSideVId].P();
                            quadrangulation.vert[existingVertexId].P() =
                                    (coord + quadrangulation.vert[existingVertexId].P())/2;
                        }

                        currentPatchSideVertex++;
                    }
                }

                currentPatchSideVertex--;
            }

            assert(currentPatchSideVertex+1 == patchSide.size());
        }

        //Internal vertices
        for (size_t i = 0; i < quadrangulatedChartMesh.vert.size(); i++) {
            if (currentVertexMap[i] == -1) {
                size_t newId = quadrangulation.vert.size();

                const typename PolyMeshType::CoordType& coord = quadrangulatedChartMesh.vert[i].P();
                vcg::tri::Allocator<PolyMeshType>::AddVertex(quadrangulation, coord);

                currentVertexMap[i] = newId;
            }
        }

        //Set faces
        for (size_t i = 0; i < quadrangulatedChartMesh.face.size(); i++) {
            assert(quadrangulatedChartMesh.face[i].VN() == 4);

            size_t newFaceId = quadrangulation.face.size();

            vcg::tri::Allocator<PolyMeshType>::AddFaces(quadrangulation, 1);

            quadrangulation.face[newFaceId].Alloc(quadrangulatedChartMesh.face[i].VN());
            for (int j = 0; j < quadrangulatedChartMesh.face[i].VN(); j++) {
                int vId = currentVertexMap[vcg::tri::Index(quadrangulatedChartMesh, quadrangulatedChartMesh.face[i].V(j))];
                assert(vId >= 0);

                quadrangulation.face[newFaceId].V(j) = &quadrangulation.vert[vId];
            }

            quadrangulationFaceLabel.push_back(chart.label);
            quadrangulationPartitions[chart.label].push_back(newFaceId);
        }

        //Fill corners vertices
        for (size_t i = 0; i < chartSides.size(); i++) {
            const ChartSide& side = chartSides[i];
            size_t vStart = side.vertices[0];

            quadrangulationCorners[chart.label].push_back(cornerVertices.at(vStart));
        }
    }

#ifdef SAVE_MESHES_FOR_DEBUG
    vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(quadrangulation, "results/quadrangulation_original.obj", vcg::tri::io::Mask::IOM_FACECOLOR);
#endif

    int numDuplicateVertices = QuadRetopology::internal::removeDuplicateVertices(quadrangulation, false);
    if (numDuplicateVertices > 0) {
        std::cout << "Warning: removed " << numDuplicateVertices << " duplicate vertices in quadrangulation." << std::endl;
    }
    int numDegenerateFaces = QuadRetopology::internal::removeDegenerateFaces(quadrangulation, false, true);
    if (numDegenerateFaces > 0) {
        std::cout << "Warning: removed " << numDegenerateFaces << " degenerate faces in quadrangulation." << std::endl;
    }
    int numUnreferencedVertices = QuadRetopology::internal::removeUnreferencedVertices(quadrangulation, false);
    if (numUnreferencedVertices > 0) {
        std::cout << "Warning: removed " << numUnreferencedVertices << " unreferenced vertices in quadrangulation." << std::endl;
    }

#ifdef SAVE_MESHES_FOR_DEBUG
    vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(quadrangulation, "results/quadrangulation_no_duplicates.obj", vcg::tri::io::Mask::IOM_FACECOLOR);
#endif

    vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(quadrangulation);
    vcg::PolygonalAlgorithm<PolyMeshType>::UpdateFaceNormalByFitting(quadrangulation);
    QuadRetopology::internal::OrientFaces<PolyMeshType>::AutoOrientFaces(quadrangulation);
    vcg::PolygonalAlgorithm<PolyMeshType>::UpdateFaceNormalByFitting(quadrangulation);
    vcg::tri::UpdateNormal<PolyMeshType>::PerVertexNormalized(quadrangulation);

    if (doubletRemoval) {
        int numDoublets = QuadRetopology::internal::removeDoublets(quadrangulation, false);
        if (numDoublets > 0) {
            std::cout << "Removed " << numDoublets << " doublets in quadrangulation." << std::endl;

            int numUnreferencedVerticesAfterDoublets = QuadRetopology::internal::removeUnreferencedVertices(quadrangulation, false);
            if (numUnreferencedVerticesAfterDoublets > 0) {
                std::cout << "Removed " << numUnreferencedVerticesAfterDoublets << " unreferenced vertices after doublet removal in quadrangulation." << std::endl;
            }
        }
    }

    //Compact and remap all info
    vcg::tri::UpdateQuality<PolyMeshType>::VertexConstant(quadrangulation, -1);
    for (size_t i=0;i<quadrangulation.vert.size();i++) {
        if (quadrangulation.vert[i].IsD())
            continue;

        quadrangulation.vert[i].Q() = i;
    }
    vcg::tri::UpdateQuality<PolyMeshType>::FaceConstant(quadrangulation, -1);
    for (size_t i=0;i<quadrangulation.face.size();i++) {
        if (quadrangulation.face[i].IsD())
            continue;

        quadrangulation.face[i].Q() = i;
    }

    PolyMeshType tmpMesh;
    vcg::tri::Append<PolyMeshType, PolyMeshType>::Mesh(tmpMesh, quadrangulation);

    std::vector<int> tmpToQuadrangulationVertex(tmpMesh.vert.size(), -1);
    std::vector<int> quadrangulationToTmpVertex(quadrangulation.vert.size(), -1);
    for (size_t i = 0; i < tmpMesh.vert.size(); i++) {
        if (tmpMesh.vert[i].IsD())
            continue;
        int id = tmpMesh.vert[i].Q();
        assert(id >= 0);
        tmpToQuadrangulationVertex[i] = id;
        quadrangulationToTmpVertex[id] = i;
    }
    std::vector<int> tmpToQuadrangulationFace(tmpMesh.face.size(), -1);
    std::vector<int> quadrangulationToTmpFace(quadrangulation.face.size(), -1);
    for (size_t i = 0; i < tmpMesh.face.size(); i++) {
        if (tmpMesh.face[i].IsD())
            continue;
        int id = tmpMesh.face[i].Q();
        assert(id >= 0);
        tmpToQuadrangulationFace[i] = id;
        quadrangulationToTmpFace[id] = i;
    }

    std::vector<int> newFaceLabel(tmpMesh.face.size(), -1);
    for (size_t i = 0; i < tmpMesh.face.size(); i++) {
        if (!tmpMesh.face[i].IsD()) {
            newFaceLabel[i] = quadrangulationFaceLabel[tmpToQuadrangulationFace[i]];
        }
    }
    quadrangulationFaceLabel = newFaceLabel;

    for (std::vector<size_t>& corners : quadrangulationCorners) {
        std::vector<size_t> newCorners;
        for (const size_t& vId : corners) {
            if (!quadrangulation.vert[vId].IsD()) {
                newCorners.push_back(quadrangulationToTmpVertex[vId]);
            }
        }
        corners = newCorners;
    }

    for (std::vector<size_t>& partition : quadrangulationPartitions) {
        std::vector<size_t> newPartition;
        for (const size_t& fId : partition) {
            if (!quadrangulation.face[fId].IsD()) {
                newPartition.push_back(quadrangulationToTmpFace[fId]);
            }
        }
        partition = newPartition;
    }


    std::vector<size_t> quadrangulationVerticesBetweenPatch;
    for (const size_t& vId : fixedVerticesSet) {
        if (!quadrangulation.vert[vId].IsD()) {
            quadrangulationVerticesBetweenPatch.push_back(quadrangulationToTmpVertex[vId]);
        }
    }

    //Create
    quadrangulation.Clear();
    vcg::tri::Append<PolyMeshType, PolyMeshType>::Mesh(quadrangulation, tmpMesh);
    vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(quadrangulation);
    vcg::tri::UpdateFlags<PolyMeshType>::FaceBorderFromFF(quadrangulation);
    vcg::tri::UpdateFlags<PolyMeshType>::VertexBorderFromFaceAdj(quadrangulation);

#ifdef SAVE_MESHES_FOR_DEBUG
    vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(quadrangulation, "results/quadrangulation_no_doublets.obj", vcg::tri::io::Mask::IOM_FACECOLOR);
#endif

    vcg::tri::UpdateNormal<TriangleMeshType>::PerFaceNormalized(newSurface);
    vcg::tri::UpdateNormal<TriangleMeshType>::PerVertexNormalized(newSurface);
    vcg::tri::UpdateBounding<TriangleMeshType>::Box(newSurface);

    vcg::GridStaticPtr<typename TriangleMeshType::FaceType,typename TriangleMeshType::FaceType::ScalarType> Grid;
    Grid.Set(newSurface.face.begin(),newSurface.face.end());

    //Reproject
    vcg::tri::UpdateBounding<PolyMeshType>::Box(quadrangulation);
    typename TriangleMeshType::ScalarType maxD=quadrangulation.bbox.Diag();
    typename TriangleMeshType::ScalarType minD=0;

    for (size_t i=0;i<quadrangulation.vert.size();i++) {
        if (quadrangulation.vert[i].IsD())
            continue;

        typename TriangleMeshType::CoordType closestPT;
        typename TriangleMeshType::FaceType *f=
                vcg::tri::GetClosestFaceBase<TriangleMeshType>(
                    newSurface,
                    Grid,
                    quadrangulation.vert[i].P(),
                    maxD,minD,
                    closestPT);

        quadrangulation.vert[i].P()=closestPT;
    }

#ifdef SAVE_MESHES_FOR_DEBUG
    vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(quadrangulation, "results/quadrangulation_reprojected.obj", vcg::tri::io::Mask::IOM_FACECOLOR);
#endif

    if (quadrangulationFixedSmoothingIterations > 0) {
        vcg::tri::UpdateSelection<PolyMeshType>::VertexAll(quadrangulation);
        for (const size_t& fixedVertexId : quadrangulationVerticesBetweenPatch) {
            if (quadrangulation.vert[fixedVertexId].IsD())
                continue;

            quadrangulation.vert[fixedVertexId].ClearS();
        }

        vcg::PolygonalAlgorithm<PolyMeshType>::template LaplacianReproject<TriangleMeshType>(quadrangulation, newSurface, quadrangulationFixedSmoothingIterations, 0.7, 0.7, true);
    }

#ifdef SAVE_MESHES_FOR_DEBUG
    vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(quadrangulation, "results/quadrangulation_smoothedfixed.obj", vcg::tri::io::Mask::IOM_FACECOLOR);
#endif

    if (quadrangulationNonFixedSmoothingIterations > 0) {
        vcg::tri::UpdateSelection<PolyMeshType>::VertexAll(quadrangulation);
        for (size_t i=0;i<quadrangulation.vert.size();i++) {
            if (quadrangulation.vert[i].IsD())
                continue;
            if (quadrangulation.vert[i].IsB()) {
                quadrangulation.vert[i].ClearS();
            }
        }

        vcg::PolygonalAlgorithm<PolyMeshType>::template LaplacianReproject<TriangleMeshType>(quadrangulation, newSurface, quadrangulationNonFixedSmoothingIterations, 0.7, 0.7, true);
    }

    vcg::PolygonalAlgorithm<PolyMeshType>::UpdateFaceNormalByFitting(quadrangulation);
    vcg::tri::UpdateNormal<PolyMeshType>::PerVertexNormalized(quadrangulation);
    vcg::tri::UpdateBounding<PolyMeshType>::Box(quadrangulation);

#ifdef SAVE_MESHES_FOR_DEBUG
    vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(quadrangulation, "results/quadrangulation_smoothednonfixed.obj", vcg::tri::io::Mask::IOM_FACECOLOR);
#endif
}


template<class PolyMeshType, class TriangleMeshType>
void computeResult(
        PolyMeshType& preservedMesh,
        PolyMeshType& quadrangulation,
        PolyMeshType& result,
        TriangleMeshType& targetBoolean,
        const bool doubletRemoval,
        const int resultSmoothingIterations,
        const double resultSmoothingNRing,
        const int resultSmoothingLaplacianIterations,
        const double resultSmoothingLaplacianNRing,
        std::vector<int>& resultPreservedVertexMap,
        std::vector<int>& resultPreservedFaceMap)
{
    vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(quadrangulation);
    vcg::tri::UpdateFlags<PolyMeshType>::FaceBorderFromFF(quadrangulation);
    vcg::tri::UpdateFlags<PolyMeshType>::VertexBorderFromFaceBorder(quadrangulation);

    for (size_t i = 0; i < preservedMesh.vert.size(); i++) {
        if (!preservedMesh.vert[i].IsD()) {
            preservedMesh.vert[i].Q() = i;
        }
    }
    for (size_t i = 0; i < quadrangulation.vert.size(); i++) {
        if (!quadrangulation.vert[i].IsD()) {
            quadrangulation.vert[i].Q() = -1;
        }
    }

    for (size_t i = 0; i < quadrangulation.face.size(); i++) {
        if (!quadrangulation.face[i].IsD()) {
            quadrangulation.face[i].Q() = -1;
        }
    }

    //Select quadrangulation borders
    vcg::tri::UpdateSelection<PolyMeshType>::FaceClear(quadrangulation);
    vcg::tri::UpdateSelection<PolyMeshType>::VertexClear(quadrangulation);
    for (size_t i = 0; i < quadrangulation.vert.size(); i++) {
        if (quadrangulation.vert[i].IsD())
            continue;

        if (quadrangulation.vert[i].IsB()) {
            quadrangulation.vert[i].SetS();
        }
    }

    //Create result
    PolyMeshType tmpMesh;

    vcg::tri::Append<PolyMeshType, PolyMeshType>::Mesh(tmpMesh, preservedMesh);
    vcg::tri::UpdateSelection<PolyMeshType>::VertexClear(tmpMesh);

    vcg::tri::Append<PolyMeshType, PolyMeshType>::Mesh(tmpMesh, quadrangulation);

#ifdef SAVE_MESHES_FOR_DEBUG
    vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(tmpMesh, "results/result_nomerged.obj", vcg::tri::io::Mask::IOM_FACECOLOR);
#endif

//    vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(tmpMesh);
//    vcg::tri::UpdateFlags<PolyMeshType>::FaceBorderFromFF(tmpMesh);
//    vcg::tri::UpdateFlags<PolyMeshType>::VertexBorderFromFaceBorder(tmpMesh);

    //Here quadrangulation borders would be selected

//    int numClustered = QuadRetopology::internal::clusterVertices(tmpMesh, true, 0.00001);
//    if (numClustered > 0) {
//        std::cout << "Clustered " << numClustered << " vertices." << std::endl;
//    }
    int numDuplicateVertices = QuadRetopology::internal::removeDuplicateVertices(tmpMesh, true);
    if (numDuplicateVertices > 0) {
        std::cout << "Merged " << numDuplicateVertices << " duplicate vertices." << std::endl;
    }
    int numDegenerateFaces = QuadRetopology::internal::removeDegenerateFaces(tmpMesh, false, true);
    if (numDegenerateFaces > 0) {
        std::cout << "Removed " << numDegenerateFaces << " degenerate faces after duplicate vertex removal." << std::endl;
    }
    int numUnreferencedVertices = QuadRetopology::internal::removeUnreferencedVertices(tmpMesh, true);
    if (numUnreferencedVertices > 0) {
        std::cout << "Removed " << numUnreferencedVertices << " unreferenced vertices after duplicate vertex removal." << std::endl;
    }

    if (doubletRemoval) {
        int numDoublets = QuadRetopology::internal::removeDoublets(tmpMesh, true);
        if (numDoublets > 0) {
            std::cout << "Removed " << numDoublets << " doublets." << std::endl;

            int numUnreferencedVerticesAfterDoublets = QuadRetopology::internal::removeUnreferencedVertices(tmpMesh, false);
            if (numUnreferencedVerticesAfterDoublets > 0) {
                std::cout << "Removed " << numUnreferencedVerticesAfterDoublets << " unreferenced vertices after doublet removal." << std::endl;
            }
        }
    }

#ifdef SAVE_MESHES_FOR_DEBUG
    vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(tmpMesh, "results/result_noduplicates.obj", vcg::tri::io::Mask::IOM_FACECOLOR);
#endif

    result.Clear();
    vcg::tri::Append<PolyMeshType, PolyMeshType>::Mesh(result, tmpMesh);

    vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(result);
    vcg::PolygonalAlgorithm<PolyMeshType>::UpdateFaceNormalByFitting(result);

//    internal::OrientFaces<PolyMeshType>::AutoOrientFaces(result);
//    vcg::PolygonalAlgorithm<PolyMeshType>::UpdateFaceNormalByFitting(result);

    vcg::tri::UpdateNormal<PolyMeshType>::PerVertexNormalized(result);

    //Fill maps
    std::vector<size_t> smoothingVertices; //Vertices to be smoothed
    resultPreservedVertexMap.resize(result.vert.size(), -1);
    for (size_t i = 0; i < result.vert.size(); i++) {
        if (result.vert[i].IsD())
            continue;

        if (result.vert[i].Q() >= 0) {
            resultPreservedVertexMap[i] = result.vert[i].Q();
        }
        else {
            smoothingVertices.push_back(i);
        }
    }
    resultPreservedFaceMap.resize(result.face.size(), -1);
    for (size_t i = 0; i < result.face.size(); i++) {
        if (result.face[i].IsD())
            continue;

        if (result.face[i].Q() >= 0) {
            resultPreservedFaceMap[i] = result.face[i].Q();
        }
    }

#ifdef SAVE_MESHES_FOR_DEBUG
    vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(result, "results/result_before_reprojection.obj", vcg::tri::io::Mask::IOM_FACECOLOR);
#endif

    vcg::tri::UpdateNormal<TriangleMeshType>::PerFaceNormalized(targetBoolean);
    vcg::tri::UpdateNormal<TriangleMeshType>::PerVertexNormalized(targetBoolean);
    vcg::tri::UpdateBounding<TriangleMeshType>::Box(targetBoolean);

    if (result.face.size() == 0)
        return;

    vcg::tri::UpdateSelection<PolyMeshType>::VertexClear(result);
    vcg::tri::UpdateSelection<PolyMeshType>::FaceClear(result);
    for (const size_t& vId : smoothingVertices) {
        result.vert[vId].SetS();
    }

    vcg::GridStaticPtr<typename TriangleMeshType::FaceType,typename TriangleMeshType::FaceType::ScalarType> Grid;
    Grid.Set(targetBoolean.face.begin(),targetBoolean.face.end());

    //Reproject
    vcg::tri::UpdateBounding<PolyMeshType>::Box(result);
    typename TriangleMeshType::ScalarType maxD=result.bbox.Diag();
    typename TriangleMeshType::ScalarType minD=0;

    for (const size_t& vId : smoothingVertices) {
        typename TriangleMeshType::CoordType closestPT;
        typename TriangleMeshType::FaceType *f=
                vcg::tri::GetClosestFaceBase<TriangleMeshType>(
                    targetBoolean,
                    Grid,
                    result.vert[vId].P(),
                    maxD,minD,
                    closestPT);

        result.vert[vId].P()=closestPT;
    }

    for (int it = 0; it < resultSmoothingIterations; it++) {
        typename PolyMeshType::ScalarType maxDistance = internal::averageEdgeLength(result) * resultSmoothingNRing;

        internal::LaplacianGeodesic(result, 1, maxDistance, 0.7);

        //Reproject
        vcg::tri::UpdateBounding<PolyMeshType>::Box(result);
        typename TriangleMeshType::ScalarType maxD=result.bbox.Diag();
        typename TriangleMeshType::ScalarType minD=0;

        for (const size_t& vId : smoothingVertices) {
            typename TriangleMeshType::CoordType closestPT;
            typename TriangleMeshType::FaceType *f=
                    vcg::tri::GetClosestFaceBase<TriangleMeshType>(
                        targetBoolean,
                        Grid,
                        result.vert[vId].P(),
                        maxD,minD,
                        closestPT);

            result.vert[vId].P()=closestPT;
        }
    }

    typename PolyMeshType::ScalarType maxDistance = internal::averageEdgeLength(result) * resultSmoothingLaplacianNRing;

    internal::LaplacianGeodesic(result, resultSmoothingLaplacianIterations, maxDistance, 0.8);

    vcg::PolygonalAlgorithm<PolyMeshType>::UpdateFaceNormalByFitting(result);
    vcg::tri::UpdateNormal<PolyMeshType>::PerVertexNormalized(result);
    vcg::tri::UpdateBounding<PolyMeshType>::Box(result);

#ifdef SAVE_MESHES_FOR_DEBUG
    vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(result, "results/result_after_reprojection.obj", vcg::tri::io::Mask::IOM_FACECOLOR);
#endif
}

template<class PolyMeshType, class TriangleMeshType>
void computeResult(
        PolyMeshType& preservedSurface,
        PolyMeshType& quadrangulatedNewSurface,
        PolyMeshType& result,
        TriangleMeshType& targetBoolean,
        const Parameters& parameters,
        std::vector<int>& resultPreservedVertexMap,
        std::vector<int>& resultPreservedFaceMap)
{
    return computeResult(
        preservedSurface,
        quadrangulatedNewSurface,
        result,
        targetBoolean,
        parameters.doubletRemoval,
        parameters.resultSmoothingIterations,
        parameters.resultSmoothingNRing,
        parameters.resultSmoothingLaplacianIterations,
        parameters.resultSmoothingLaplacianNRing,
        resultPreservedVertexMap,
        resultPreservedFaceMap);
}

}
