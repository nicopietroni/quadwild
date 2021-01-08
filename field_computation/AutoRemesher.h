#ifndef AUTOREMESHER_H
#define AUTOREMESHER_H

#include <vcg/complex/allocate.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/crease_cut.h>
#include <vcg/complex/algorithms/isotropic_remeshing.h>

#include <vcg/complex/append.h>

#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/algorithms/closest.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse.h>

#include <memory>

template <class Mesh>
class AutoRemesher {

    typedef typename Mesh::ScalarType ScalarType;
    typedef typename Mesh::CoordType  CoordType;

    typedef typename Mesh::VertexType    VertexType;
    typedef typename Mesh::VertexPointer VertexPointer;

    typedef typename Mesh::FaceType    FaceType;
    typedef typename Mesh::FacePointer FacePointer;

    typedef vcg::GridStaticPtr<FaceType, ScalarType> StaticGrid;

    static ScalarType computeAR(Mesh & m, const double perc = 0.05)
    {
        vcg::tri::ForEachFace(m, [] (FaceType & f) {
            f.Q() = vcg::QualityRadii(f.cP(0), f.cP(1), f.cP(2));
        });

        vcg::Histogram<ScalarType> hist;
        vcg::tri::Stat<Mesh>::ComputePerFaceQualityHistogram(m, hist);

        //		    return hist.MinV();
        return hist.Percentile(perc);
    }

    static void collapseSurvivingMicroEdges(Mesh & m, const ScalarType qualityThr = 0.001, const ScalarType edgeRatio = 0.025, const int maxIter = 2)
    {
        typedef vcg::tri::BasicVertexPair<VertexType> VertexPair;
        typedef vcg::tri::EdgeCollapser<Mesh, VertexPair> Collapser;
        typedef typename vcg::face::Pos<FaceType> PosType;

        int count = 0; int iter = 0;
        do
        {
            count = 0;
            vcg::tri::UpdateTopology<Mesh>::VertexFace(m);

            for(auto fi=m.face.begin(); fi!=m.face.end(); ++fi)
                if(!(*fi).IsD())
                {
                    if(vcg::QualityRadii(fi->cP(0), fi->cP(1), fi->cP(2)) <= qualityThr)
                    {
                        ScalarType minEdgeLength = std::numeric_limits<ScalarType>::max();
                        ScalarType maxEdgeLength = 0;

                        int minEdge = 0, maxEdge = 0;
                        for(auto i=0; i<3; ++i)
                        {
                            const ScalarType len = vcg::Distance(fi->cP0(i), fi->cP1(i)) ;
                            if (len < minEdgeLength)
                            {
                                minEdge = i;
                                minEdgeLength = len;
                            }
                            if (len > maxEdgeLength)
                            {
                                maxEdge = i;
                                maxEdgeLength = len;
                            }
                        }

                        if (minEdgeLength <= maxEdgeLength * edgeRatio)
                        {
                            PosType pi(&*fi, minEdge);
                            VertexPair  bp = VertexPair(fi->V0(minEdge), fi->V1(minEdge));
                            CoordType mp = (fi->cP0(minEdge) + fi->cP1(minEdge))/2.f;

                            if(Collapser::LinkConditions(bp))
                            {
                                Collapser::Do(m, bp, mp, true);
                            }
                        }
                    }
                }
        } while (count > 0 && ++iter < maxIter);
    }


    static void removeColinearFaces(Mesh & m, const ScalarType colinearThr = 0.001)
    {
        vcg::tri::UpdateTopology<Mesh>::FaceFace(m);

        Mesh projectMesh;
        vcg::tri::Append<Mesh, Mesh>::MeshCopy(projectMesh, m);
        vcg::tri::UpdateBounding<Mesh>::Box(projectMesh);
        vcg::tri::UpdateNormal<Mesh>::PerVertexNormalizedPerFace(projectMesh);
        StaticGrid grid;
        grid.Set(projectMesh.face.begin(), projectMesh.face.end());


        int count = 0;
        int iter = 0;
        do
        {
            vcg::tri::UpdateTopology<Mesh>::FaceFace(m);
            vcg::tri::UnMarkAll(m);

            count = 0;
            for (size_t i = 0; i < size_t(m.FN()); ++i)
            {
                FaceType & f = m.face[i];

                ScalarType quality = vcg::QualityRadii(f.cP(0), f.cP(1), f.cP(2));

                if (quality <= colinearThr)
                {
                    //find longest edge
                    double edges[3];
                    edges[0] = vcg::Distance(f.cP(0), f.cP(1));
                    edges[1] = vcg::Distance(f.cP(1), f.cP(2));
                    edges[2] = vcg::Distance(f.cP(2), f.cP(0));

                    ScalarType smallestEdge = std::min(edges[0], std::min(edges[1], edges[2]));
                    int longestIdx = std::find(edges, edges+3, std::max(std::max(edges[0], edges[1]), edges[2])) - (edges);

                    if (vcg::tri::IsMarked(m, f.V2(longestIdx)))
                        continue;


                    auto f1 = f.cFFp(longestIdx);
                    vcg::tri::Mark(m,f.V2(longestIdx));
                    if (!vcg::face::IsBorder(f, longestIdx) && vcg::face::IsManifold(f, longestIdx) && vcg::face::checkFlipEdgeNotManifold<FaceType>(f, longestIdx))  {

                        // Check if EdgeFlipping improves quality
                        FacePointer g = f.FFp(longestIdx); int k = f.FFi(longestIdx);
                        vcg::Triangle3<ScalarType> t1(f.P(longestIdx), f.P1(longestIdx), f.P2(longestIdx)), t2(g->P(k), g->P1(k), g->P2(k)),
                                t3(f.P(longestIdx), g->P2(k), f.P2(longestIdx)), t4(g->P(k), f.P2(longestIdx), g->P2(k));

                        auto n1 = vcg::TriangleNormal(t1);
                        auto n2 = vcg::TriangleNormal(t2);
                        auto n3 = vcg::TriangleNormal(t3);
                        auto n4 = vcg::TriangleNormal(t4);

                        auto biggestSmallest = vcg::DoubleArea(t1) > vcg::DoubleArea(t2) ? std::make_pair(t1, t2) : std::make_pair(t2, t1);
                        auto areaRatio = vcg::DoubleArea(biggestSmallest.first) / vcg::DoubleArea(biggestSmallest.second);

                        bool normalCheck = true;
                        //                        if (n1.Norm() > 0.001 && n2.Norm() > 0.001)
                        {
                            auto referenceNormal = vcg::NormalizedTriangleNormal(biggestSmallest.first);

                            normalCheck &= vcg::NormalizedTriangleNormal(t3) * referenceNormal >= 0.95;
                            normalCheck &= vcg::NormalizedTriangleNormal(t4) * referenceNormal >= 0.95;
                        }

                        bool areaCheck = false;
                        if (areaRatio > 1000)
                        {
                            areaCheck |= vcg::DoubleArea(t3) / vcg::DoubleArea(biggestSmallest.second) > 1000 && vcg::DoubleArea(t4) / vcg::DoubleArea(biggestSmallest.second) > 1000;
                        }

                        if ((normalCheck) && (areaCheck || std::min( QualityFace(t1), QualityFace(t2) ) <= std::min( QualityFace(t3), QualityFace(t4))))
                        {
                            ScalarType dist;
                            CoordType closest;
                            auto fp0 = vcg::tri::GetClosestFaceBase(projectMesh, grid, vcg::Barycenter(t3), smallestEdge/4., dist, closest);
                            if (fp0 == NULL)
                                continue;

                            auto fp1 = vcg::tri::GetClosestFaceBase(projectMesh, grid, vcg::Barycenter(t4), smallestEdge/4., dist, closest);
                            if (fp1 == NULL)
                                continue;

                            vcg::face::FlipEdgeNotManifold<FaceType>(f, longestIdx);
                            ++count;
                        }
                    }
                }
            }
        } while (count && ++iter < 75);
    }

    static void specialTreatment(Mesh & m, Mesh & project, typename vcg::tri::IsotropicRemeshing<Mesh>::Params & para)
    {
        std::cerr << "[REMESH] Mesh requires special treatment" << std::endl;
        std::cerr << "[REMESH] Special Treatment. Faces: " << m.FN() << " Initial quality: " << computeAR(m) << std::endl;

        vcg::tri::UpdateFlags<Mesh>::Clear(m);
        CleanMesh(m);

        vcg::tri::UpdateSelection<Mesh>::FaceClear(m);
        vcg::tri::UpdateTopology<Mesh>::FaceFace(m);
        vcg::tri::UpdateNormal<Mesh>::PerFaceNormalized(m);
        vcg::tri::ForEachFace(m, [] (FaceType & f) {
            float QQ = 0;
            if (f.Q() <= 0.15)
                f.SetS();
            else
            {
                for (int j = 0; j < 3; j++)
                {
                    auto adjf = f.FFp(j);
                    float nangle = vcg::AngleN((*adjf).cN(), f.cN());
                    if (nangle > 130)
                    {
                        f.SetS();
                        break;
                    }
                }
            }
        });

        vcg::tri::UpdateSelection<Mesh>::FaceDilate(m);

        para.iter = 10;
        para.selectedOnly = true;

        para.aspectRatioThr = 0.05;
        para.cleanFlag = false;

        ScalarType edgeL = para.maxLength;

        para.SetTargetLen(edgeL * 0.75);
        para.maxSurfDist = m.bbox.Diag() / 100.;

        para.userSelectedCreases = false;
        para.SetFeatureAngleDeg(50);

        vcg::tri::Smooth<Mesh>::VertexCoordPlanarLaplacian(m, 20, vcg::math::ToRad(0.5), true);
        vcg::tri::IsotropicRemeshing<Mesh>::Do(m, para);

        vcg::tri::UpdateFlags<Mesh>::Clear(m);

        para.iter= 10;
        para.SetTargetLen(edgeL * 1.5);
        para.cleanFlag = true;
        para.userSelectedCreases = false;
        para.selectedOnly = false;
        para.maxSurfDist = m.bbox.Diag()/500.;
        para.SetFeatureAngleDeg(30);

        vcg::tri::IsotropicRemeshing<Mesh>::Do(m, para);

        std::cerr << "[REMESH] Special Treatment. Faces: " << m.FN() << " Post-Treatment quality: " << computeAR(m) << std::endl;
    }
public:

    typedef struct Params {

        int iterations = 15;

        ScalarType targetAspect = 0.35;
        int targetDeltaFN = 5000;

        int initialApproximateFN = 20000;

        ScalarType creaseAngle = 25.;

        bool userSelectedCreases = true;
        bool surfDistCheck = true;

        int erodeDilate = 0;

    } Params;

    static size_t openNonManifoldEdges(Mesh & m, const ScalarType moveThreshold)
    {

        vcg::tri::UpdateTopology<Mesh>::FaceFace(m);
        vcg::tri::UpdateTopology<Mesh>::VertexFace(m);
        vcg::tri::UpdateFlags<Mesh>::VertexClearV(m);

        std::cout << "Opening non-manifold edges...";
        std::cout << "mesh starts with " << vcg::tri::Clean<Mesh>::CountNonManifoldEdgeFF(m) << std::endl;

        typedef typename vcg::face::Pos<FaceType> PosType;

        typedef std::vector<std::vector<std::pair<size_t, size_t> > > VertexToFaceGroups;

        std::unordered_map<size_t,VertexToFaceGroups> map;

        vcg::tri::ForEachFacePos(m, [&](PosType & pos) {

            if (!pos.V()->IsV() && !pos.IsManifold())
            {
                pos.V()->SetV();

                std::vector<FacePointer> faceVec;
                std::vector<int> vIndices;
                vcg::face::VFStarVF(pos.V(), faceVec, vIndices);

                std::unordered_set<size_t> inserted;

                VertexToFaceGroups faceGroups;

                for (size_t i = 0; i < faceVec.size(); ++i)
                {
                    const FacePointer fp = faceVec[i];

                    size_t fidx = vcg::tri::Index(m, fp);
                    if (inserted.count(fidx) != 0)
                        continue;

                    std::vector<std::pair<size_t, size_t> > manifoldGroup;

                    PosType cyclePos(fp, vIndices[i]);
                    PosType beginPos = cyclePos;
                    //get to a non manifold edge...
                    do
                    {
                        manifoldGroup.push_back(std::make_pair(vcg::tri::Index(m, cyclePos.F()), cyclePos.VInd()));
                        inserted.insert(vcg::tri::Index(m, cyclePos.F()));
                        cyclePos.F()->Q() = faceGroups.size() + 1;
                        cyclePos.FlipE();
                        cyclePos.NextF();
                    } while (cyclePos.IsManifold() && !cyclePos.IsBorder() && cyclePos != beginPos);

                    if (cyclePos != beginPos)
                    {
                        cyclePos = beginPos;
                        cyclePos.NextF();

                        while (cyclePos.IsManifold() && !cyclePos.IsBorder())
                        {
                            manifoldGroup.push_back(std::make_pair(vcg::tri::Index(m, cyclePos.F()), cyclePos.VInd()));
                            inserted.insert(vcg::tri::Index(m, cyclePos.F()));
                            cyclePos.F()->Q() = faceGroups.size() + 1;
                            cyclePos.FlipE();
                            cyclePos.NextF();
                        }
                    }

                    faceGroups.push_back(manifoldGroup);
                }

                map[vcg::tri::Index(m, pos.V())] = faceGroups;
            }
        });

        for (auto group : map)
        {
            const size_t vert = group.first;
            const VertexToFaceGroups & faceGroups = group.second;
            if (faceGroups.size() > 1)
            {
                auto vp = vcg::tri::Allocator<Mesh>::AddVertices(m, faceGroups.size()-1);

                for (size_t i = 1; i < faceGroups.size(); ++i)
                {
                    vp->P() = m.vert[vert].cP();

                    CoordType delta(0, 0, 0);
                    for (std::pair<size_t,size_t> faceVertIndex : faceGroups[i])
                    {
                        m.face[faceVertIndex.first].V(faceVertIndex.second) = &*vp;
                        delta += vcg::Barycenter(m.face[faceVertIndex.first]) - vp->cP();
                    }
                    delta /= faceGroups[i].size();
                    vp->P() += delta * moveThreshold;
                    ++vp;
                }
            }
        }

        return map.size();
    }

    static void SplitNonManifold (Mesh & m)
    {
        int splitV = 0, openings = 0;
        do
        {
            vcg::tri::UpdateTopology<Mesh>::FaceFace(m);
            splitV = vcg::tri::Clean<Mesh>::SplitNonManifoldVertex(m, 0.0001);
            openings = openNonManifoldEdges(m, 0.0001);
            std::cout << "Opened " << openings << " non manifold edges" << std::endl;
        } while (splitV > 0 || openings > 0);

        vcg::tri::Clean<Mesh>::RemoveUnreferencedVertex(m);
        vcg::tri::Allocator<Mesh>::CompactEveryVector(m);
    }

    static std::shared_ptr<Mesh> CleanMesh(Mesh & m, const bool & splitNonManifold = false)
    {
        std::shared_ptr<Mesh> ret = std::make_shared<Mesh>();
        vcg::tri::Append<Mesh, Mesh>::MeshCopy(*ret, m);
        vcg::tri::UpdateTopology<Mesh>::FaceFace(*ret);

        int i = 0; int removed = 0;
        do
        {
            removed = vcg::tri::Clean<Mesh>::SplitNonManifoldVertex(*ret, 0.0001);
            std::cout << "removed " << removed << " non manifold vertices..." << std::endl;
        }
        while (removed > 0 && i < 50);

        std::cout << "Removed " << vcg::tri::Clean<Mesh>::RemoveDuplicateFace(*ret) << " duplicate faces..." << std::endl;
        std::cout << "Removed " << vcg::tri::Clean<Mesh>::RemoveUnreferencedVertex(*ret) << " unreferenced vertices..." << std::endl;
        vcg::tri::Allocator<Mesh>::CompactEveryVector(*ret);

        vcg::tri::UpdateTopology<Mesh>::FaceFace(*ret);

        std::cout << "[AutoRemeshCleaner] Removing colinear faces by flip..." << std::endl;
        removeColinearFaces(*ret);

        std::cout << "Removed " << vcg::tri::Clean<Mesh>::RemoveZeroAreaFace(*ret) << " zero area faces..." << std::endl;
        vcg::tri::Clean<Mesh>::RemoveUnreferencedVertex(*ret);
        vcg::tri::Allocator<Mesh>::CompactEveryVector(*ret);


        if (splitNonManifold)
            SplitNonManifold(*ret);


        vcg::tri::UpdateTopology<Mesh>::FaceFace(*ret);

        std::cout << "[AutoRemeshCleaner] Attempting to coherently orient faces..." << std::endl;
        bool oriented = false, orientable = false;
        vcg::tri::Clean<Mesh>::OrientCoherentlyMesh(*ret, oriented, orientable);

        std::cout << "[AutoRemeshCleaner] Orientable:" << orientable << std::endl;
        std::cout << "[AutoRemeshCleaner] Oriented:"   << oriented << std::endl;

        return ret;
    }

    //for big meshes disabling par.surfDistCheck provides big perf improvements, sacrificing result accuracy
    static std::shared_ptr<Mesh> Remesh (Mesh & m, Params & par)
    {
        std::shared_ptr<Mesh> ret = std::make_shared<Mesh>();

        vcg::tri::UpdateBounding<Mesh>::Box(m);
        vcg::tri::UpdateTopology<Mesh>::FaceFace(m);

        typename vcg::tri::IsotropicRemeshing<Mesh>::Params para;
        para.iter = par.iterations;
        para.SetFeatureAngleDeg(par.creaseAngle);
        para.splitFlag    = true;
        para.swapFlag     = true;
        para.collapseFlag = true;
        para.smoothFlag   = true;
        para.projectFlag  = true;
        para.selectedOnly = false;
        para.adapt=false;

        para.aspectRatioThr = 0.3;
        para.cleanFlag = false;

        para.maxSurfDist = m.bbox.Diag() / 2500.;
        para.surfDistCheck = m.FN() < 400000 ? par.surfDistCheck : false;
        para.userSelectedCreases = true;//par.userSelectedCreases;

        ScalarType prevFN = m.FN();
        ScalarType deltaFN = m.FN();

        ScalarType aspect = 0;
        ScalarType edgeLow  = 0;
        ScalarType edgeL = std::sqrt(2.309 * vcg::tri::Stat<Mesh>::ComputeMeshArea(m) / par.initialApproximateFN);//m.bbox.Diag() * 0.025;//std::sqrt(vcg::tri::Stat<Mesh>::ComputeMeshArea(m) * 2 / par.initialApproximateFN);

        ScalarType edgeHigh = edgeL * 2.;

        para.SetTargetLen(edgeL);

        vcg::tri::Append<Mesh, Mesh>::MeshCopy(*ret, m);

        ret->UpdateDataStructures();
        ret->InitSharpFeatures(par.creaseAngle);
        ret->ErodeDilate(par.erodeDilate);

        vcg::tri::IsotropicRemeshing<Mesh>::Do(*ret, para);
        std::cout << "Iter: "<<  0 << " faces: " << ret->FN() << " quality: " <<  computeAR(*ret) << std::endl;

        ret->UpdateDataStructures();
        ret->InitSharpFeatures(par.creaseAngle);
        ret->ErodeDilate(par.erodeDilate);

        para.SetTargetLen(edgeL * 0.85);
        para.adapt = true;
        para.maxSurfDist = m.bbox.Diag() / 2500.;

        vcg::tri::IsotropicRemeshing<Mesh>::Do(*ret, para);
        auto quality = computeAR(*ret);
        std::cout << "Iter: "<<  1 << " faces: " << ret->FN() << " quality: " << quality << std::endl;


        if (quality < 0.1)
        {
            ret->UpdateDataStructures();
            ret->InitSharpFeatures(par.creaseAngle);
            ret->ErodeDilate(par.erodeDilate);

            int count = 0;
            vcg::tri::ForEachFace(*ret, [&] (FaceType & f) {

                if (f.cQ() < 0.01)
                {
                    ++count;
                    for (int i = 0; i < 3; ++i)
                    {
                        f.ClearFaceEdgeS(i);
                    }
                }
            });

            std::cout << count << " faces were relieved of crease edges due to poor quality" << std::endl;
            para.SetTargetLen(edgeL * 0.6);
            para.adapt = true;
            para.maxSurfDist = m.bbox.Diag() / 200.;
            vcg::tri::IsotropicRemeshing<Mesh>::Do(*ret, para);
            auto quality = computeAR(*ret);
            std::cout << "Iter: "<<  3 << " faces: " << ret->FN() << " quality: " << quality << std::endl;
        }


        /*
        bool forcedExit = false;
        int countIterations = 0;
        do {
            ++countIterations;
            vcg::tri::Append<Mesh, Mesh>::MeshCopy(*ret, m);

            para.SetTargetLen(edgeL);
#ifdef DEBUG
            std::cout << edgeLow << "--" << edgeL << "--" << edgeHigh << std::endl;
#endif

            vcg::tri::IsotropicRemeshing<Mesh>::Do(*ret, para);

            aspect = computeAR(*ret);
            ScalarType newFN = ret->FN();
        std::cout << "Iter: "<<  countIterations << " faces: " << newFN << std::endl;
            deltaFN = std::abs(ret->FN() - prevFN);
            prevFN = newFN;

        if (aspect >= par.targetAspect && newFN >= 30000 && newFN <= 100000)
        break;

            if (aspect >= par.targetAspect)
            {
                edgeLow = edgeL;
            } else {
                edgeHigh = edgeL;
            }

            edgeL = (edgeHigh + edgeLow) / 2.;

            if (ret->FN() > 200000 && aspect < par.targetAspect * 0.5)
            {
                specialTreatment(*ret, m, para);

                break;
                //vcg::tri::Append<Mesh, Mesh>::MeshCopy(m, *ret);
            }

#ifdef DEBUG
            std::cout << "Delta FN " << deltaFN << std::endl;
            std::cout << "Mesh FN  " << ret->FN() << std::endl;
            std::cout << "Min AR   " << computeAR(*ret) << std::endl;
#endif
        }
        while (aspect < par.targetAspect || deltaFN > par.targetDeltaFN);
*/
        std::cerr << "[REMESH] RemeshedFaces:" << ret->FN() << std::endl;
        std::cerr << "[REMESH] RemeshedAspect:" << computeAR(*ret) << std::endl;

        const ScalarType thr = 0.01;
        collapseSurvivingMicroEdges(*ret, thr);

        vcg::tri::UpdateSelection<Mesh>::FaceClear(*ret);
        vcg::tri::UpdateTopology<Mesh>::FaceFace(*ret);

        vcg::tri::ForEachFace(*ret, [&] (FaceType & f) {
            if (!f.IsD() && vcg::QualityRadii(f.cP(0), f.cP(1), f.cP(2)) <= 0.01)
                f.SetS();
        });

        int zeroArea = 0;
        do
        {
            zeroArea = vcg::tri::Clean<Mesh>::RemoveZeroAreaFace(*ret);
            std::cout << "removed " << zeroArea << " zero area faces " << std::endl;
        }
        while(zeroArea != 0);

        vcg::tri::UpdateSelection<Mesh>::FaceDilate(*ret);
        //		vcg::tri::UpdateSelection<Mesh>::FaceDilate(*ret);
        //		vcg::tri::Smooth<Mesh>::VertexCoordLaplacian(*ret, 15, true);

        vcg::tri::Allocator<Mesh>::CompactEveryVector(*ret);
        std::cout << "[REMESH] remeshing ends.." << std::endl;
        return ret;
    }
};

#endif // AUTOREMESHER_H
