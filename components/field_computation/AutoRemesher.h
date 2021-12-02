/***************************************************************************/
/* Copyright(C) 2021


The authors of

Reliable Feature-Line Driven Quad-Remeshing
Siggraph 2021


 All rights reserved.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
****************************************************************************/

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
            //std::cout<<"Q:"<<f.Q()<<std::endl;
        });

        vcg::Histogram<ScalarType> hist;
        vcg::tri::Stat<Mesh>::ComputePerFaceQualityHistogram(m, hist);

        //		    return hist.MinV();
        return hist.Percentile(perc);
    }
public:

    static bool collapseSurvivingMicroEdges(Mesh & m, const ScalarType qualityThr = 0.001, const ScalarType edgeRatio = 0.025, const int maxIter = 2)
    {
        typedef vcg::tri::BasicVertexPair<VertexType> VertexPair;
        typedef vcg::tri::EdgeCollapser<Mesh, VertexPair> Collapser;
        typedef typename vcg::face::Pos<FaceType> PosType;

        bool done=false;
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

                        //                        if (minEdgeLength <= maxEdgeLength * edgeRatio)
                        //                        {
                        PosType pi(&*fi, minEdge);

                        //                        //select the vertices
                        //                        (*fi).V(0)->SetS();
                        //                        (*fi).V(1)->SetS();
                        //                        (*fi).V(2)->SetS();
                        //                        (*fi).FFp(minEdge)->V(0)->SetS();
                        //                        (*fi).FFp(minEdge)->V(1)->SetS();
                        //                        (*fi).FFp(minEdge)->V(2)->SetS();

                        VertexPair  bp = VertexPair(fi->V0(minEdge), fi->V1(minEdge));
                        CoordType mp = (fi->cP0(minEdge) + fi->cP1(minEdge))/2.f;

                        if(Collapser::LinkConditions(bp))
                        {
                            Collapser::Do(m, bp, mp, true);
                            done=true;
                        }
                        //                        }
                    }
                }
        } while (count > 0 && ++iter < maxIter);
        return done;
    }
private:
    static void MakeEdgeSelConsistent(Mesh & m)
    {
        for (size_t i=0;i<m.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (m.face[i].IsFaceEdgeS(j))
                {
                    if (vcg::face::IsBorder(m.face[i],j))continue;
                    FaceType *f1=m.face[i].FFp(j);
                    int EdgeI=m.face[i].FFi(j);
                    f1->SetFaceEdgeS(EdgeI);
                }
            }
    }

    static void SelectAllBoundaryV(Mesh & m)
    {
        for (size_t i=0;i<m.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (vcg::face::IsBorder(m.face[i],j)){
                    m.face[i].V0(j)->SetS();
                    m.face[i].V1(j)->SetS();
                    continue;
                }

                if (m.face[i].IsFaceEdgeS(j))
                {
                    m.face[i].V0(j)->SetS();
                    m.face[i].V1(j)->SetS();
                    FaceType *f1=m.face[i].FFp(j);
                    int EdgeI=m.face[i].FFi(j);
                    f1->SetFaceEdgeS(EdgeI);
                }
            }
    }

public:

    typedef struct Params {
        int iterations = 15;
        ScalarType targetAspect = 0.35;
        int targetDeltaFN = 5000;
        int initialApproximateFN = 20000;
        ScalarType creaseAngle = 25.;
        //bool userSelectedCreases = true;
        bool surfDistCheck = true;
        int erodeDilate = 0;
        ScalarType minAdaptiveMult = 0.3;
        ScalarType maxAdaptiveMult = 3;
        ScalarType minAspectRatioThr = 0.05;
        ScalarType targetEdgeLen = 0;
    } Params;

    static size_t openNonManifoldEdges(Mesh & m, const ScalarType moveThreshold,
                                       bool debugMesg=false)
    {

        vcg::tri::UpdateTopology<Mesh>::FaceFace(m);
        vcg::tri::UpdateTopology<Mesh>::VertexFace(m);
        vcg::tri::UpdateFlags<Mesh>::VertexClearV(m);

        if (debugMesg)
        {
            std::cout << "Opening non-manifold edges...";
            std::cout << "mesh starts with " << vcg::tri::Clean<Mesh>::CountNonManifoldEdgeFF(m) << std::endl;
        }

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


//    static ScalarType ExpectedEdgeL(const Mesh & m,
//                                    size_t TargetSph=2000,
//                                    size_t MinFaces=15000)
//    {
//        ScalarType Vol=m.Volume();
//        ScalarType A=m.Area();
//        ScalarType FaceA=A/TargetSph;
//        //radius and volume of a sphere
//        ScalarType KScale=(Vol/A)*(3/(m.bbox.Diag()/3.4));
//        //ScalarType KScale=(A*(m.bbox.Diag()/3.4))/(3*Vol);
//        ScalarType IdealA=FaceA*KScale;
//        ScalarType IdealL0=sqrt(IdealA*2.309);
//        ScalarType IdealL1=sqrt((A*2.309)/MinFaces);
//        std::cout<<"KScale"<<KScale<<std::endl;
//        //exit(0);
//        return std::max(IdealL0,IdealL1);
//    }

    static ScalarType ExpectedEdgeL(const Mesh & m,
                                    size_t TargetSph=2000,
                                    size_t MinFaces=10000)
    {
        ScalarType Vol=m.Volume();
        ScalarType A=m.Area();
        ScalarType FaceA=A/TargetSph;
        //radius and volume of a sphere
        ScalarType Sphericity=(pow(M_PI,1.0/3.0)*pow((6.0*Vol),2.0/3.0))/A;
        ScalarType KScale=pow(Sphericity,2);
        //ScalarType KScale=(A*(m.bbox.Diag()/3.4))/(3*Vol);
        //ScalarType KScale=(pow(A,1.5))/(3*Vol);
        ScalarType IdealA=FaceA*KScale;
        ScalarType IdealL0=sqrt(IdealA*2.309);
        ScalarType IdealL1=sqrt((A*2.309)/MinFaces);
        std::cout<<"KScale "<<KScale<<std::endl;
        //exit(0);
        return std::min(IdealL0,IdealL1);
        //return IdealL0;
    }

    //remove sharp features that have been removed after remesher or clean that were still marked
    static void UpdateCoherentSharp(Mesh & m, Params & par)
    {
        if (par.creaseAngle<=0)return;
        m.UpdateDataStructures();
        //std::set<std::pair<CoordType,CoordType> > Features;
        for (size_t i=0;i<m.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!m.face[i].IsFaceEdgeS(j))continue;

                ScalarType angle = DihedralAngleRad(m.face[i],j);
                if(fabs(angle)<vcg::math::ToRad(par.creaseAngle))
                {
                    if (vcg::face::IsBorder(m.face[i],j))continue;
                    m.face[i].ClearFaceEdgeS(j);
                    m.face[i].FFp(j)->ClearFaceEdgeS(m.face[i].FFi(j));
                }
            }

        m.InitFeatureCoordsTable();
        //MeshPrepocess<Mesh>::InitSharpFeatures(mesh,BPar.sharp_feature_thr,BPar.feature_erode_dilate);

    }

    //for big meshes disabling par.surfDistCheck provides big perf improvements, sacrificing result accuracy
    //static std::shared_ptr<Mesh> Remesh (Mesh & m, Params & par)
    static void RemeshAdapt(Mesh & m, Params & par)
    {

        vcg::tri::UpdateBounding<Mesh>::Box(m);
        vcg::tri::UpdateTopology<Mesh>::FaceFace(m);

        typename vcg::tri::IsotropicRemeshing<Mesh>::Params para;
        para.iter = par.iterations;
        //para.SetFeatureAngleDeg(par.creaseAngle);

        para.splitFlag    = true;
        para.swapFlag     = true;
        para.collapseFlag = true;
        para.smoothFlag   = true;
        para.projectFlag  = true;
        para.selectedOnly = false;
        para.adapt=false;
        para.aspectRatioThr = 0.3;
        para.cleanFlag = true;

        para.minAdaptiveMult = par.minAdaptiveMult;
        para.maxAdaptiveMult = par.maxAdaptiveMult;

        para.maxSurfDist = m.bbox.Diag() / 2500.;
        para.surfDistCheck = m.FN() < 400000 ? par.surfDistCheck : false;
        para.userSelectedCreases = true;



        ScalarType edgeL = ExpectedEdgeL(m);//,par.initialApproximateFN);//std::sqrt(2.309 * vcg::tri::Stat<Mesh>::ComputeMeshArea(m) / par.initialApproximateFN);//m.bbox.Diag() * 0.025;//std::sqrt(vcg::tri::Stat<Mesh>::ComputeMeshArea(m) * 2 / par.initialApproximateFN);

        if (par.targetEdgeLen == 0)
            par.targetEdgeLen = edgeL;

        para.SetTargetLen(par.targetEdgeLen);


        std::cout << "Before Remeshing - faces: " << m.FN() << " quality: " <<  computeAR(m) << std::endl;
        vcg::tri::IsotropicRemeshing<Mesh>::Do(m, para);
        std::cout << "After Iter 0 - faces: " << m.FN() << " quality: " <<  computeAR(m) << std::endl;


        const ScalarType thr = 0.01;

        //vcg::tri::UpdateSelection<Mesh>::VertexClear(m);
        collapseSurvivingMicroEdges(m,thr);

        UpdateCoherentSharp(m,par);

        para.adapt = true;
        para.smoothFlag   = true;
        para.maxSurfDist = m.bbox.Diag() / 2500.;

        vcg::tri::IsotropicRemeshing<Mesh>::Do(m, para);

        m.UpdateDataStructures();

        std::cout << "After Iter 1 - faces: " << m.FN() << " quality: " <<  computeAR(m) << std::endl;

//        MakeEdgeSelConsistent(m);
//        SelectAllBoundaryV(m);
//        typename Local_Param_Smooth<Mesh>::UVSmoothParam UVP;
//        UVP.FixSel=true;

//        Local_Param_Smooth<Mesh>::Smooth(m,UVP);
//        m.UpdateDataStructures();
//        std::cout << "After Iter 2 - faces: " << m.FN() << " quality: " <<  computeAR(m) << std::endl;

    }


    //    static int SelectToRemesh(Mesh & m,ScalarType minR=0.2,size_t dilate=3)
    //    {
    //        vcg::tri::UpdateSelection<Mesh>::FaceClear(m);
    //        for (size_t i=0;i<m.face.size();i++)
    //        {
    //            m.face[i].Q() = vcg::QualityRadii(m.face[i].cP(0),
    //                                              m.face[i].cP(1),
    //                                              m.face[i].cP(2));
    //            if (m.face[i].Q()<minR)
    //                m.face[i].SetS();
    //       }
    ////        //then check the borders
    ////        vcg::tri::UpdateSelection<Mesh>::VertexClear(m);
    ////        for (size_t i=0;i<m.face.size();i++)
    ////            for (size_t j=0;j<3;j++)
    ////            {
    ////                if (!m.face[i].IsFaceEdgeS(j))continue;
    ////                m.face[i].V0(j)->SetS();
    ////                m.face[i].V1(j)->SetS();
    ////            }
    ////        for (size_t i=0;i<m.face.size();i++)
    ////            for (size_t j=0;j<3;j++)
    ////            {
    ////                if (m.face[i].IsFaceEdgeS(j))continue;
    ////                if (!m.face[i].V0(j)->IsS())continue;
    ////                if (!m.face[i].V1(j)->IsS())continue;
    ////                m.face[i].SetS();
    ////            }

    //        vcg::tri::UpdateSelection<Mesh>::VertexClear(m);
    //        for (size_t s=0;s<dilate;s++)
    //        {
    //            vcg::tri::UpdateSelection<Mesh>::VertexFromFaceLoose(m);
    //            vcg::tri::UpdateSelection<Mesh>::FaceFromVertexLoose(m);
    //        }

    //        int numS=0;
    //        for (size_t i=0;i<m.face.size();i++)
    //            if (m.face[i].IsS())numS++;
    //        return (numS);
    //    }

    static int NumBadTris(Mesh & m,ScalarType minR=0.2,size_t dilate=3)
    {
        int numS=0;
        for (size_t i=0;i<m.face.size();i++)
        {
            m.face[i].Q() = vcg::QualityRadii(m.face[i].cP(0),
                                              m.face[i].cP(1),
                                              m.face[i].cP(2));
            if (m.face[i].Q()<minR)
                numS++;
        }
        return (numS);
    }

    //    static void Remesh2(Mesh & m, Params & par)
    //    {

    //        vcg::tri::UpdateBounding<Mesh>::Box(m);
    //        vcg::tri::UpdateTopology<Mesh>::FaceFace(m);

    //        typename vcg::tri::IsotropicRemeshing<Mesh>::Params para;
    //        para.iter = par.iterations;
    //        //para.SetFeatureAngleDeg(par.creaseAngle);

    //        para.splitFlag    = true;
    //        para.swapFlag     = true;
    //        para.collapseFlag = true;
    //        para.smoothFlag   = false;
    //        para.projectFlag  = false;
    //        para.selectedOnly = false;
    //        para.adapt=false;
    //        para.aspectRatioThr = 0.3;
    //        para.cleanFlag = true;
    //        para.minAdaptiveMult = par.minAdaptiveMult;
    //        para.maxAdaptiveMult = par.maxAdaptiveMult;

    //        para.maxSurfDist = m.bbox.Diag() / 2500.;
    //        para.surfDistCheck = m.FN() < 400000 ? par.surfDistCheck : false;
    //        para.userSelectedCreases = true;

    //        ScalarType edgeL = ExpectedEdgeL(m);

    //        //if (par.targetEdgeLen == 0)
    //        par.targetEdgeLen = edgeL;

    //        para.SetTargetLen(par.targetEdgeLen);

    //        std::cout << "Before Remeshing - faces: " << m.FN() << " quality: " <<  computeAR(m) << std::endl;
    //        vcg::tri::IsotropicRemeshing<Mesh>::Do(m, para);
    //        std::cout << "After Iter 0 - faces: " << m.FN() << " quality: " <<  computeAR(m) << std::endl;


    //        MakeEdgeSelConsistent(m);
    //        SelectAllBoundaryV(m);
    //        typename Local_Param_Smooth<Mesh>::UVSmoothParam UVP;
    //        UVP.FixSel=true;

    //        Local_Param_Smooth<Mesh>::Smooth(m,UVP);


    //        //const ScalarType thr = 0.01;
    //        collapseSurvivingMicroEdges(m, 0.01);

    //        m.UpdateDataStructures();

    //        int MaxS=3;
    //        int currS=0;
    //        int NumSel0=0;
    //        int NumSel1=0;
    //        do{
    //            //             MakeEdgeSelConsistent(m);
    //            //             SelectAllBoundaryV(m);
    //            //             typename Local_Param_Smooth<Mesh>::UVSmoothParam UVP;
    //            //             UVP.FixSel=true;

    //            //Local_Param_Smooth<Mesh>::Smooth(m,UVP);

    //            NumSel0=NumBadTris(m);
    //            std::cout << "Not Nice 0: " <<  NumSel0 << " Faces"<<std::endl;


    //            if (NumSel0>0)
    //            {
    //                par.targetEdgeLen*=0.75;
    //                para.SetTargetLen(par.targetEdgeLen);
    //                para.maxSurfDist = m.bbox.Diag() / 2500.;
    //                //para.selectedOnly=true;
    //                para.adapt=true;
    //                para.smoothFlag= true;
    //                vcg::tri::IsotropicRemeshing<Mesh>::Do(m, para);
    //                collapseSurvivingMicroEdges(m, 0.01);
    //                m.UpdateDataStructures();
    //            }

    //            NumSel1=NumBadTris(m);
    //            std::cout << "Not Nice 1: " <<  NumSel1 << " Faces"<<std::endl;
    //            currS++;

    //        }while ((NumSel1<NumSel0) && (currS<MaxS));
    //        //        std::cout << "Performed: " <<  currS << " Steps"<<std::endl;
    //        std::cout << "After Iter Adapt - faces: " << m.FN() << " quality: " <<  computeAR(m) << std::endl;

    //    }

    //    //for big meshes disabling par.surfDistCheck provides big perf improvements, sacrificing result accuracy
    //    static void Remesh (Mesh & m, Params & par)
    //    {

    //        m.UpdateDataStructures();

    ////        vcg::tri::UpdateBounding<Mesh>::Box(m);
    ////        vcg::tri::UpdateTopology<Mesh>::FaceFace(m);

    //        typename vcg::tri::IsotropicRemeshing<Mesh>::Params para;
    //        para.iter = par.iterations;
    //        para.SetFeatureAngleDeg(par.creaseAngle);
    //        para.splitFlag    = true;
    //        para.swapFlag     = true;
    //        para.collapseFlag = true;
    //        para.smoothFlag   = false;
    //        para.projectFlag  = true;
    //        para.selectedOnly = false;
    //        para.adapt=false;
    //        para.aspectRatioThr = 0.3;
    //        para.cleanFlag = false;

    //        para.minAdaptiveMult = par.minAdaptiveMult;
    //        para.maxAdaptiveMult = par.maxAdaptiveMult;

    //        para.maxSurfDist = m.bbox.Diag() / 2500.;
    //        para.surfDistCheck = m.FN() < 400000 ? par.surfDistCheck : false;
    //        para.userSelectedCreases = true;//par.userSelectedCreases;

    ////        ScalarType prevFN = m.FN();
    ////        ScalarType deltaFN = m.FN();

    ////        ScalarType aspect = 0;
    ////        ScalarType edgeLow  = 0;
    //        ScalarType edgeL = std::sqrt(2.309 * vcg::tri::Stat<Mesh>::ComputeMeshArea(m) / par.initialApproximateFN);//m.bbox.Diag() * 0.025;//std::sqrt(vcg::tri::Stat<Mesh>::ComputeMeshArea(m) * 2 / par.initialApproximateFN);

    //        if (par.targetEdgeLen == 0)
    //            par.targetEdgeLen = edgeL;

    ////        ScalarType edgeHigh = edgeL * 2.;

    ////        para.SetTargetLen(par.targetEdgeLen);

    ////        vcg::tri::Append<Mesh, Mesh>::MeshCopy(*ret, m);

    ////        ret->UpdateDataStructures();

    ////        ret->InitSharpFeatures(par.creaseAngle);
    ////        ret->ErodeDilate(par.erodeDilate);


    //        vcg::tri::IsotropicRemeshing<Mesh>::Do(m, para);
    //        std::cout << "Iter: "<<  0 << " faces: " <<m.FN() << " quality: " <<  computeAR(m) << std::endl;


    //        //const ScalarType thr = 0.01;
    //        collapseSurvivingMicroEdges(m, 0.01);

    ////        ret->UpdateDataStructures();
    ////        ret->InitSharpFeatures(par.creaseAngle);
    ////        ret->ErodeDilate(par.erodeDilate);
    ////        MP.InitSharpFeatures(par.creaseAngle);
    ////        MP.ErodeDilate(par.erodeDilate);

    //        //para.SetTargetLen(par.targetEdgeLen * 0.85);
    //        para.adapt = true;
    //        para.smoothFlag   = true;
    //        para.maxSurfDist = m.bbox.Diag() / 2500.;

    //        vcg::tri::IsotropicRemeshing<Mesh>::Do(m, para);
    //        auto quality = computeAR(m);
    //        std::cout << "Iter: "<<  1 << " faces: " << m.FN() << " quality: " << quality << std::endl;

    //        std::cerr << "[REMESH] RemeshedFaces:" << m.FN() << std::endl;
    //        std::cerr << "[REMESH] RemeshedAspect:" << computeAR(m) << std::endl;

    //        const ScalarType thr = 0.01;
    //        collapseSurvivingMicroEdges(m, thr);

    //        vcg::tri::UpdateSelection<Mesh>::FaceClear(m);
    //        vcg::tri::UpdateTopology<Mesh>::FaceFace(m);

    //        vcg::tri::ForEachFace(m, [&] (FaceType & f) {
    //            if (!f.IsD() && vcg::QualityRadii(f.cP(0), f.cP(1), f.cP(2)) <= 0.01)
    //                f.SetS();
    //        });

    //        int zeroArea = 0;
    //        do
    //        {
    //            zeroArea = vcg::tri::Clean<Mesh>::RemoveZeroAreaFace(m);
    //            std::cout << "removed " << zeroArea << " zero area faces " << std::endl;
    //        }
    //        while(zeroArea != 0);

    //        vcg::tri::UpdateSelection<Mesh>::FaceDilate(m);
    //        //		vcg::tri::UpdateSelection<Mesh>::FaceDilate(*ret);
    //        //		vcg::tri::Smooth<Mesh>::VertexCoordLaplacian(*ret, 15, true);

    //        vcg::tri::Allocator<Mesh>::CompactEveryVector(m);
    //        std::cout << "[REMESH] remeshing ends.." << std::endl;
    //        //return ret;
    //    }
};

#endif // AUTOREMESHER_H
