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

#ifndef PATCH_MANAGER
#define PATCH_MANAGER

#include "mesh_type.h"
#include "patch_tracer.h"
#include <wrap/igl/arap_parametrization.h>
#include <wrap/igl/lscm_parametrization.h>
#include "candidate_path.h"
#include "vert_field_graph.h"
#include <vcg/complex/algorithms/parametrization/distortion.h>
#include <vcg/space/outline2_packer.h>
#include "catmull_clarkability.h"
#include <vcg/space/triangle3.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/attribute_seam.h>
//#include <vcg/complex/algorithms/crease_cut.h>
#include <vcg/space/intersection2.h>

//#define MIN_SAMPLES_HARD 4
//#define MAX_SAMPLES 1000
//#define MIN_SAMPLES 50
//#define MAX_NARROW_CONST 0.05
//#define NARROW_NEED 1
//#define MAX_TWIN_DIKSTRA 1
//#define MIN_ADMITTIBLE 3
//#define MAX_ADMITTIBLE 6
//#define MAX_BORDER_SAMPLE 8
//#define MIN_BORDER_SAMPLE 1

class PatchGeneralParameters
{

public:

    static size_t &MaxTwinDikstra()
    {
        static size_t NNeed=1;
        assert(NNeed>=0);
        return NNeed;
    }

    static size_t &NarrowNeed()
    {
        static size_t NNeed=1;
        assert((NNeed==1)||(NNeed==2));
        return NNeed;
    }

    static size_t &MinSampleHard()
    {
        static size_t MinSampleHard=4;
        return MinSampleHard;
    }

    static typename TraceMesh::ScalarType &MaxNarrowConst()
    {
        static typename TraceMesh::ScalarType MaxNarr=0.05;
        return MaxNarr;
    }

    static size_t &MaxSamples()
    {
        static size_t MaxSamples=1000;
        return MaxSamples;
    }

    static size_t &MinSamples()
    {
        static size_t MinSamples=50;
        return MinSamples;
    }

    static size_t &MinBorderSample()
    {
        static size_t MinSample=1;
        return MinSample;
    }

    static size_t &MaxBorderSample()
    {
        static size_t MaxSample=8;
        return MaxSample;
    }

    static size_t MinAdmittable()
    {
        return 3;
    }

    static size_t &MaxAdmittable()
    {
        static size_t MaxAdmit=6;
        return MaxAdmit;
    }
};

#define NARROW_NEED PatchGeneralParameters::NarrowNeed()
#define MIN_SAMPLES_HARD PatchGeneralParameters::MinSampleHard()
#define MAX_SAMPLES PatchGeneralParameters::MaxSamples()
#define MIN_SAMPLES PatchGeneralParameters::MinSamples()
#define MAX_ADMITTIBLE PatchGeneralParameters::MaxAdmittable()
#define MIN_ADMITTIBLE PatchGeneralParameters::MinAdmittable()
#define MAX_BORDER_SAMPLE PatchGeneralParameters::MaxBorderSample()
#define MIN_BORDER_SAMPLE PatchGeneralParameters::MinBorderSample()
#define MAX_NARROW_CONST PatchGeneralParameters::MaxNarrowConst()
#define MAX_TWIN_DIKSTRA PatchGeneralParameters::MaxTwinDikstra()

template <class MeshType>
typename MeshType::ScalarType MeshArea(MeshType &mesh)
{
    typedef typename MeshType::ScalarType ScalarType;

    ScalarType currA=0;
    for (size_t i=0;i<mesh.face.size();i++)
        currA+=(vcg::DoubleArea(mesh.face[i])/2.0);
    return currA;
}

//template <class MeshType>
//void SplitMeshAlongEdgeSel(MeshType &mesh)
//{
//   // typedef typename MeshType::ScalarType ScalarType;

//    MeshType dstMesh;
//	vcg::tri::AttributeSeam::SplitVertex(srcMesh, dstMesh, ExtractVertex, CompareVertex, CopyVertex);

////    ScalarType currA=0;
////    for (size_t i=0;i<mesh.face.size();i++)
////        currA+=(vcg::DoubleArea(mesh.face[i])/2.0);
////    return currA;
//}

enum ParamType{Arap,LSQMap};

template <class ScalarType>
struct PatchInfo
{
    int NumEmitters;
    int NumCorners;
    int Genus;
    std::vector<ScalarType> CurvedL;
    bool CClarkability;
    int ExpectedValence;
    int NumSing;
    ScalarType QualityFunctorValue;

    PatchInfo()
    {
        NumEmitters=0;
        NumCorners=0;
        Genus=0;
        CClarkability=false;//std::numeric_limits<ScalarType>::max();
        ExpectedValence=-1;
        NumSing=-1;
        QualityFunctorValue=0;
    }
};


template <class MeshType>
class PatchManager
{
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename vcg::face::Pos<FaceType> PosType;

public:

    static void SelectPatchFacesVert(MeshType &totMesh,const std::vector<size_t> &PatchFaces)
    {
        vcg::tri::UpdateSelection<MeshType>::Clear(totMesh);
        for (size_t i=0;i<PatchFaces.size();i++)
            totMesh.face[PatchFaces[i]].SetS();

        vcg::tri::UpdateSelection<MeshType>::VertexFromFaceLoose(totMesh);
    }

    static int ExpectedValence(MeshType &mesh,
                               const std::vector<size_t> &PatchFaces,
                               const std::vector<size_t> &PatchCorners,
                               bool &OnCorner)
    {
        int ExpVal=4;
        vcg::tri::UnMarkAll<MeshType>(mesh);
        //vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);
        for (size_t i=0;i<PatchFaces.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                VertexType *v=mesh.face[PatchFaces[i]].V(j);
                if (vcg::tri::IsMarked(mesh,v))continue;
                vcg::tri::Mark(mesh,v);
                //if (v->IsV())continue;
                //v->SetV();

                if (v->SingularityValence==4)continue;
                if (ExpVal!=4){OnCorner=false;return -1;}//multiple singularities

                //check if it is a corner or not
                size_t IndexV=vcg::tri::Index(mesh,v);
                std::vector<size_t>::const_iterator it = find (PatchCorners.begin(), PatchCorners.end(),IndexV);
                if (it != PatchCorners.end())
                    OnCorner=true;//it is a corner of the patch
                else
                    OnCorner=false;
                ExpVal=v->SingularityValence;
            }
        }
        //vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);
        return ExpVal;
    }

    static int NumSingularities(MeshType &mesh,
                                const std::vector<size_t> &PatchFaces,
                                const std::vector<size_t> &PatchCorners)
    {
        int NumSing=0;
        //vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);
        vcg::tri::UnMarkAll<MeshType>(mesh);
        for (size_t i=0;i<PatchFaces.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                VertexType *v=mesh.face[PatchFaces[i]].V(j);
                if (vcg::tri::IsMarked(mesh,v))continue;
                vcg::tri::Mark(mesh,v);
                //if (v->IsV())continue;
                //v->SetV();
                if (v->SingularityValence==4)continue;

                size_t IndexV=vcg::tri::Index(mesh,v);
                std::vector<size_t>::const_iterator it = find (PatchCorners.begin(), PatchCorners.end(),IndexV);
                if (it != PatchCorners.end())continue;
                NumSing++;
            }
        }
        //vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);
        return NumSing;
    }


    //    static int PatchGenus(MeshType &mesh,const std::vector<size_t> &PatchFaces)
    //    {

    //        //    SelectPatchFacesVert<MeshType>(mesh,PatchFaces);

    //        //    MeshType subM;
    //        //    vcg::tri::Append<MeshType,MeshType>::Mesh(subM,mesh,true);
    //        //    vcg::tri::Clean<MeshType>::RemoveDuplicateVertex(subM);
    //        //    vcg::tri::Allocator<MeshType>::CompactEveryVector(subM);

    //        //    vcg::tri::UnMarkAll(mesh);

    //        //std::unordered_set<std::pair<size_t,size_t> > EdgeSet;
    //        std::set<std::pair<CoordType,CoordType> > EdgeSet;
    //        std::set<CoordType> VertSet;
    //        size_t NumF=PatchFaces.size();
    //        size_t NumV=0;
    //        size_t NumE=0;
    //        for (size_t i=0;i<PatchFaces.size();i++)
    //        {
    //            FaceType *f=&mesh.face[PatchFaces[i]];
    //            for (size_t j=0;j<3;j++)
    //            {
    //                //count the vertex
    //                //if (!vcg::tri::IsMarked(mesh,f->V0(j)))
    //                if (!VertSet.count(f->P0(j)))
    //                {
    //                    //vcg::tri::Mark(mesh,f->V0(j));
    //                    VertSet.insert(f->P0(j));
    //                    NumV++;
    //                }
    //                //            size_t IndV0=vcg::tri::Index(mesh,f->V0(j));
    //                //            size_t IndV1=vcg::tri::Index(mesh,f->V1(j));
    //                //            EdgeSet.insert(std::pair<size_t,size_t>(std::min(IndV0,IndV1),std::max(IndV0,IndV1)));
    //                CoordType P0=f->P0(j);
    //                CoordType P1=f->P1(j);
    //                EdgeSet.insert(std::pair<CoordType,CoordType>(std::min(P0,P1),std::max(P0,P1)));
    //            }
    //        }
    //        //    for (size_t i=0;i<subM.vert.size();i++)
    //        //    {
    //        //        if (subM.vert[i].IsD())continue;
    //        //        if (subM.vert[i].IsS())NumV++;
    //        //    }

    //        NumE=EdgeSet.size();
    //        return ( NumV + NumF - NumE );
    //    }

    //    static int PatchGenus(MeshType &mesh,const std::vector<size_t> &PatchFaces)
    //    {

    //        vcg::tri::UnMarkAll(mesh);

    //        //std::unordered_set<std::pair<size_t,size_t> > EdgeSet;
    //        std::set<std::pair<CoordType,CoordType> > EdgeSet;
    //        //std::set<CoordType> VertSet;
    //        size_t NumF=PatchFaces.size();
    //        size_t NumV=0;
    //        size_t NumE=0;

    //        for (size_t i=0;i<PatchFaces.size();i++)
    //            vcg::tri::Mark(mesh,&mesh.face[PatchFaces[i]]);

    //        for (size_t i=0;i<PatchFaces.size();i++)
    //        {
    //            FaceType *f0=&mesh.face[PatchFaces[i]];
    //            for (size_t j=0;j<3;j++)
    //            {
    //                //count the vertex
    //                if (!vcg::tri::IsMarked(mesh,f0->V0(j)))
    //                {
    //                    NumV++;
    //                    vcg::tri::Mark(mesh,f0->V0(j));
    //                }
    //                bool AddEdge=vcg::face::IsBorder((*f0),j);
    //                FaceType *f1=f0->FFp(j);
    //                AddEdge|=(!vcg::tri::IsMarked(mesh,f1));
    //                AddEdge|=(f1<f0);
    //                if (AddEdge)
    //                    NumE++;
    //                //                if (f1==f)
    //                //                if (!VertSet.count(f->P0(j)))
    //                //                {
    //                //                    //vcg::tri::Mark(mesh,f->V0(j));
    //                //                    VertSet.insert(f->P0(j));
    //                //                    NumV++;
    //                //                }
    //                //            size_t IndV0=vcg::tri::Index(mesh,f->V0(j));
    //                //            size_t IndV1=vcg::tri::Index(mesh,f->V1(j));
    //                //            EdgeSet.insert(std::pair<size_t,size_t>(std::min(IndV0,IndV1),std::max(IndV0,IndV1)));
    //                //                CoordType P0=f->P0(j);
    //                //                CoordType P1=f->P1(j);
    //                //                EdgeSet.insert(std::pair<CoordType,CoordType>(std::min(P0,P1),std::max(P0,P1)));
    //            }
    //        }
    //        //    for (size_t i=0;i<subM.vert.size();i++)
    //        //    {
    //        //        if (subM.vert[i].IsD())continue;
    //        //        if (subM.vert[i].IsS())NumV++;
    //        //    }

    //        //NumE=EdgeSet.size();
    //        return ( NumV + NumF - NumE );
    //    }

    static int FastPatchGenus(MeshType &mesh,const std::vector<size_t> &PatchFaces)
    {

        //vcg::tri::UnMarkAll(mesh);

        //std::unordered_set<std::pair<size_t,size_t> > EdgeSet;
        std::set<std::pair<CoordType,CoordType> > EdgeSet;
        std::set<CoordType> VertSet;

        size_t NumF=PatchFaces.size();

        //        for (size_t i=0;i<PatchFaces.size();i++)
        //            vcg::tri::Mark(mesh,&mesh.face[PatchFaces[i]]);

        for (size_t i=0;i<PatchFaces.size();i++)
        {
            FaceType *f0=&mesh.face[PatchFaces[i]];
            for (size_t j=0;j<3;j++)
            {
                CoordType P0=f0->P0(j);
                CoordType P1=f0->P1(j);
                std::pair<CoordType,CoordType> EdgePos(std::min(P0,P1),std::max(P0,P1));
                EdgeSet.insert(EdgePos);
                VertSet.insert(P0);
                VertSet.insert(P1);
            }
        }
        size_t NumV=VertSet.size();
        size_t NumE=EdgeSet.size();

        return ( NumV + NumF - NumE );
    }

    static int PatchGenusCopyMesh(MeshType &mesh,const std::vector<size_t> &PatchFaces)
    {
        MeshType PatchMesh;
        GetMeshFromPatch(mesh,PatchFaces,PatchMesh,true);
        int GenusVal=vcg::tri::Clean<MeshType>::CountHoles(PatchMesh);
        return GenusVal;
    }

    static int PatchGenus(MeshType &mesh,
                          const std::vector<size_t> &PatchFaces,
                          bool allowDarts,
                          bool allowSelfGlued)
    {

        if ((allowDarts)||(allowSelfGlued))
            return (PatchGenusCopyMesh(mesh,PatchFaces));
        else
            return FastPatchGenus(mesh,PatchFaces);
    }


    static void ComputeUV(MeshType &mesh, ParamType &PType,bool fixSVert)
    {
        if (PType==Arap)
        {
            if (fixSVert)
                vcg::tri::OptimizeUV_ARAP(mesh,100,MeshType::VertexType::SELECTED,true);
            else
                vcg::tri::OptimizeUV_ARAP(mesh,100,0,true);
        }
        else
        {
            if (fixSVert)
                vcg::tri::OptimizeUV_LSCM(mesh,MeshType::VertexType::SELECTED);
            else
                vcg::tri::OptimizeUV_LSCM(mesh,0);
        }
    }


    static void GetCornersUV(size_t numV,
                             const std::vector<ScalarType> &EdgeL,
                             std::vector<vcg::Point2<ScalarType> > &CornerUV)
    {
        ScalarType MaxV=0.45;
        ScalarType totLen=0;
        for (size_t i=0;i<EdgeL.size();i++)
            totLen+=EdgeL[i];

        std::vector<ScalarType> RatioL(EdgeL.size(),0);
        for (size_t i=0;i<EdgeL.size();i++)
            RatioL[i]=EdgeL[i]/totLen;

        for (size_t i=0;i<EdgeL.size();i++)
            RatioL[i]=std::min(RatioL[i],MaxV);

        CornerUV.clear();
        assert(numV>=3);

        //derive regular borders
        //ScalarType angle_step=(2*M_PI)/(ScalarType)numV;
        //std::vector<ScalarType> angleVect;
        for (size_t i=0;i<RatioL.size();i++)
            RatioL[i]*=(2*M_PI);

        ScalarType curr_angle=0;
        for (size_t i=0;i<numV;i++)
        {
            ScalarType U=cos(curr_angle);
            ScalarType V=sin(curr_angle);
            CornerUV.push_back(vcg::Point2<ScalarType>(U,V));
            curr_angle+=RatioL[i];//angle_step;
        }
    }

    static void ParametrizeCorners(MeshType &mesh,bool scaleEdges,
                                   std::vector<std::vector<size_t> > BorderSeq,
                                   const std::vector<size_t> &CornersIDX)
    {
        typedef typename MeshType::ScalarType ScalarType;
        std::vector<vcg::Point2<ScalarType> > CornerUV;
        std::vector<ScalarType> BorderLen;

        assert(CornersIDX.size()==BorderSeq.size());
        if (scaleEdges)
            GetLenghts(mesh,BorderSeq,BorderLen);
        else
            BorderLen=std::vector<ScalarType>(CornersIDX.size(),1);

        GetCornersUV(CornersIDX.size(),BorderLen,CornerUV);
        for (size_t i=0;i<CornersIDX.size();i++)
        {
            size_t IndexV=CornersIDX[i];
            assert(IndexV>=0);
            assert(IndexV<mesh.vert.size());
            mesh.vert[IndexV].T().P()=CornerUV[i];
        }
    }


    static void SelectVertices(MeshType &mesh,const std::vector<size_t> &CornersIDX)
    {
        vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);
        for (size_t i=0;i<CornersIDX.size();i++)
        {
            size_t IndexV=CornersIDX[i];
            assert(IndexV>=0);
            assert(IndexV<mesh.vert.size());
            mesh.vert[IndexV].SetS();
        }
    }


    static void SelectVertices(MeshType &mesh,const std::vector<std::vector<size_t> > &CornersIDX)
    {
        vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);
        for (size_t i=0;i<CornersIDX.size();i++)
            for (size_t j=0;j<CornersIDX[i].size();j++)
            {
                size_t IndexV=CornersIDX[i][j];
                assert(IndexV>=0);
                assert(IndexV<mesh.vert.size());
                mesh.vert[IndexV].SetS();
            }
    }

    static void SelectFaces(MeshType &mesh,const std::vector<size_t> &FacesIDX)
    {
        vcg::tri::UpdateSelection<MeshType>::FaceClear(mesh);
        for (size_t i=0;i<FacesIDX.size();i++)
        {
            size_t IndexF=FacesIDX[i];
            assert(IndexF>=0);
            assert(IndexF<mesh.face.size());
            mesh.face[IndexF].SetS();
        }
    }

    static void SelectFaces(MeshType &mesh,const std::vector<std::vector<size_t> > &FacesIDX)
    {
        vcg::tri::UpdateSelection<MeshType>::FaceClear(mesh);
        for (size_t i=0;i<FacesIDX.size();i++)
            for (size_t j=0;j<FacesIDX[i].size();j++)
            {
                size_t IndexF=FacesIDX[i][j];
                assert(IndexF>=0);
                assert(IndexF<mesh.face.size());
                mesh.face[IndexF].SetS();
            }
    }

    static void GetIndexFromQ(MeshType &mesh,
                              const std::vector<size_t> &CornersIdxQ,
                              std::vector<size_t> &CornersIdx)
    {
        for (size_t i=0;i<CornersIdxQ.size();i++)
        {
            bool found=false;
            size_t currIdx=0;
            for (size_t j=0;j<mesh.vert.size();j++)
            {
                if (mesh.vert[j].Q()==CornersIdxQ[i])
                {
                    currIdx=j;
                    found=true;
                    break;
                }
            }
            assert(found);
            CornersIdx.push_back(currIdx);
        }
    }


    static void DeriveBorderSeq(MeshType &mesh,std::vector<std::vector<size_t> > &BorderSeq)
    {
        typedef typename vcg::face::Pos<typename MeshType::FaceType> PosType;
        BorderSeq.clear();
        //set the initial pos
        PosType InitPos;
        bool found=false;
        for (size_t i=0;i<mesh.face.size();i++)
        {
            if (found)break;
            for (size_t j=0;j<3;j++)
            {
                if (!vcg::face::IsBorder(mesh.face[i],j))continue;
                InitPos=PosType(&mesh.face[i],j);
                if (InitPos.VFlip()->IsS())
                {
                    found=true;
                    break;
                }
            }
        }
        //vcg::tri::io::ExporterPLY<MeshType>::Save(mesh,"testBorder.ply",vcg::tri::io::Mask::IOM_FLAGS);
        assert(found);
        assert(InitPos.IsBorder());
        //assert(InitPos.VFlip()->IsS());

        PosType CurrPos=InitPos;
        BorderSeq.resize(1);
        size_t IndexV=vcg::tri::Index(mesh,CurrPos.VFlip());
        BorderSeq.back().push_back(IndexV);
        //std::cout<<"A"<<std::endl;
        do{
            size_t IndexV=vcg::tri::Index(mesh,CurrPos.V());
            BorderSeq.back().push_back(IndexV);
            CurrPos.NextB();
            if ((CurrPos!=InitPos)&&(CurrPos.VFlip()->IsS()))
            {
                BorderSeq.resize(BorderSeq.size()+1);
                size_t IndexV=vcg::tri::Index(mesh,CurrPos.VFlip());
                BorderSeq.back().push_back(IndexV);
            }
        }while (CurrPos!=InitPos);
        //std::cout<<"B"<<std::endl;
    }

    static bool DeriveFirstSeqPos(MeshType &mesh,const std::vector<size_t> &Faces,
                                  PosType &InitPos,const VertexType *startingV=NULL)
    {
        for (size_t i=0;i<Faces.size();i++)
        {
            size_t IndexF=Faces[i];
            for (size_t j=0;j<3;j++)
            {
                //if (!vcg::face::IsBorder(mesh.face[IndexF],j))continue;
                if (mesh.face[IndexF].IsF(j))continue;
                InitPos=PosType(&mesh.face[IndexF],j);
                if ((InitPos.VFlip()->IsS())&&(startingV==NULL))
                {
                    //found=true;
                    return true;
                }
                if ((InitPos.VFlip()==startingV)&&(startingV!=NULL))
                {
                    return true;
                }
            }
        }
        return false;
    }

    static void DeriveFauxSeqPos(MeshType &mesh,const std::vector<size_t> &Faces,
                                 std::vector<std::vector<PosType> > &BorderSeqPos,
                                 const size_t &startingV)
    {
        BorderSeqPos.clear();

        assert(mesh.vert[startingV].IsS());
        //set the initial pos
        PosType InitPos;
        bool found=DeriveFirstSeqPos(mesh,Faces,InitPos,&mesh.vert[startingV]);

        assert(found);
        assert(!InitPos.IsFaux());

        PosType CurrPos=InitPos;
        BorderSeqPos.resize(1);
        do{
            BorderSeqPos.back().push_back(CurrPos);
            CurrPos.NextNotFaux();
            if (CurrPos==InitPos)return;
            if(CurrPos.VFlip()->IsS())
                BorderSeqPos.resize(BorderSeqPos.size()+1);
        }while (true);
    }

    //CHECK FOR CORRECTNESS
    static void DeriveFauxSeq(MeshType &mesh,
                              const std::vector<size_t> &Faces,
                              std::vector<std::vector<size_t> > &BorderSeq)
    {
        typedef typename vcg::face::Pos<typename MeshType::FaceType> PosType;
        BorderSeq.clear();
        //set the initial pos
        PosType InitPos;
        bool found=false;
        //for (size_t i=0;i<mesh.face.size();i++)
        for (size_t i=0;i<Faces.size();i++)
        {
            size_t IndexF=Faces[i];
            if (found)break;
            for (size_t j=0;j<3;j++)
            {
                //if (!vcg::face::IsBorder(mesh.face[IndexF],j))continue;
                if (mesh.face[IndexF].IsF(j))continue;
                InitPos=PosType(&mesh.face[IndexF],j);
                if (InitPos.VFlip()->IsS())
                {
                    found=true;
                    break;
                }
            }
        }
        //vcg::tri::io::ExporterPLY<MeshType>::Save(mesh,"testBorder.ply",vcg::tri::io::Mask::IOM_FLAGS);
        assert(found);
        assert(!InitPos.IsFaux());
        //assert(InitPos.VFlip()->IsS());

        PosType CurrPos=InitPos;
        BorderSeq.resize(1);

        size_t IndexV=vcg::tri::Index(mesh,CurrPos.VFlip());
        BorderSeq.back().push_back(IndexV);
        //std::cout<<"A"<<std::endl;
        do{
            size_t IndexV=vcg::tri::Index(mesh,CurrPos.V());
            BorderSeq.back().push_back(IndexV);

            CurrPos.NextNotFaux();
            if ((CurrPos!=InitPos)&&(CurrPos.VFlip()->IsS()))
            {
                BorderSeq.resize(BorderSeq.size()+1);

                size_t IndexV=vcg::tri::Index(mesh,CurrPos.VFlip());
                BorderSeq.back().push_back(IndexV);
            }
        }while (CurrPos!=InitPos);
        //std::cout<<"B"<<std::endl;
    }

    static void SetFacePartitionOnFaceQ(MeshType &mesh,const std::vector<std::vector<size_t> > &PatchFaces)
    {
        vcg::tri::UpdateQuality<MeshType>::FaceConstant(mesh,-1);
        for (size_t i=0;i<PatchFaces.size();i++)
            for (size_t j=0;j<PatchFaces[i].size();j++)
            {
                assert(PatchFaces[i][j]<mesh.face.size());
                mesh.face[PatchFaces[i][j]].Q()=i;
            }
    }


    static ScalarType FieldLenght(//const MeshType &mesh,
                                  const std::vector<size_t> &BorderSeq,
                                  //std::unordered_map<std::pair<size_t,size_t>,typename MeshType::ScalarType>  &EdgeMap)
                                  std::map<std::pair<size_t,size_t>,typename MeshType::ScalarType> &EdgeMap)
    {
        typedef typename MeshType::ScalarType ScalarType;
        ScalarType currL=0;
        for (size_t i=0;i<BorderSeq.size()-1;i++)
        {
            size_t IndexV0=BorderSeq[i];
            size_t IndexV1=BorderSeq[i+1];
            std::pair<size_t,size_t> Key(std::min(IndexV0,IndexV1),
                                         std::max(IndexV0,IndexV1));
            assert(EdgeMap.count(Key)>0);
            currL+=EdgeMap[Key];
        }
        return currL;
    }


    static ScalarType Lenght(const MeshType &mesh,
                             const std::vector<size_t> &BorderSeq)
    {
        typedef typename MeshType::ScalarType ScalarType;
        typedef typename MeshType::CoordType CoordType;
        ScalarType currL=0;
        for (size_t i=0;i<BorderSeq.size()-1;i++)
        {
            CoordType P0=mesh.vert[BorderSeq[i]].P();
            CoordType P1=mesh.vert[BorderSeq[i+1]].P();
            currL+=(P1-P0).Norm();
        }
        return currL;
    }


    static void GetLenghts(const MeshType &mesh,
                           const std::vector<std::vector<size_t> > &BorderSeq,
                           std::vector<ScalarType> &BorderLen)
    {
        for (size_t i=0;i<BorderSeq.size();i++)
            BorderLen.push_back(Lenght(mesh,BorderSeq[i]));
    }


    static void ParametrizeSeq(MeshType &mesh,const std::vector<size_t> &BorderSeq)
    {

        vcg::Point2<ScalarType> T0=mesh.vert[BorderSeq[0]].T().P();
        vcg::Point2<ScalarType> T1=mesh.vert[BorderSeq.back()].T().P();
        ScalarType sideL=Lenght(mesh,BorderSeq);
        //ScalarType ScaleFact=(T0-T1).Norm()/SideL;
        ScalarType currL=0;
        for (size_t i=0;i<BorderSeq.size()-2;i++)
        {
            size_t Idx0=BorderSeq[i];
            size_t Idx1=BorderSeq[i+1];
            CoordType P0=mesh.vert[Idx0].P();
            CoordType P1=mesh.vert[Idx1].P();
            currL+=(P1-P0).Norm();
            ScalarType Interp=currL/sideL;
            assert(Interp<1);
            mesh.vert[Idx1].T().P()=T1*Interp+T0*(1-Interp);
        }
    }


    static void ParametrizeSeq(MeshType &mesh,const std::vector<std::vector<size_t> > &BorderSeq)
    {
        for (size_t i=0;i<BorderSeq.size();i++)
            ParametrizeSeq(mesh,BorderSeq[i]);
    }


    static void ComputeUV(MeshType &mesh, ParamType PType,
                          bool fixBorders,bool scaleEdges,
                          const std::vector<size_t> &CornersIdx,
                          std::vector<size_t> &SortedCorner)
    {
        if (fixBorders)
        {
            SortedCorner.clear();
            SelectVertices(mesh,CornersIdx);
            std::vector<std::vector<size_t> > BorderSeq;
            DeriveBorderSeq(mesh,BorderSeq);
            for (size_t i=0;i<BorderSeq.size();i++)
            {
                assert(BorderSeq[i].size()>=2);
                SortedCorner.push_back(BorderSeq[i][0]);
            }

            ParametrizeCorners(mesh,scaleEdges,BorderSeq,SortedCorner);

            ParametrizeSeq(mesh,BorderSeq);
            vcg::tri::UpdateSelection<MeshType>::VertexFromBorderFlag(mesh);
        }

        ComputeUV(mesh,PType,fixBorders);
    }



    static void ComputeParametrizedSubMesh(MeshType &mesh,
                                           MeshType &subMesh,
                                           const std::vector<size_t> &PatchFaces,
                                           const std::vector<size_t> &PatchCorners,
                                           std::vector<size_t> &SortedCorners,
                                           ParamType PType,
                                           bool FixBorders,
                                           bool ScaleEdges,
                                           bool InternalCuts)
    {
        for (size_t i=0;i<mesh.vert.size();i++)
            mesh.vert[i].Q()=i;
        for (size_t i=0;i<mesh.face.size();i++)
            mesh.face[i].Q()=i;

        GetMeshFromPatch(mesh,PatchFaces,subMesh,InternalCuts);
        std::vector<size_t> CornersIdx;
        GetIndexFromQ(subMesh,PatchCorners,CornersIdx);
        ComputeUV(subMesh,PType,FixBorders,ScaleEdges,CornersIdx,SortedCorners);
    }


    static void SetUVtoPos(MeshType &mesh)
    {
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            mesh.vert[i].P().X()=mesh.vert[i].T().P().X();
            mesh.vert[i].P().Y()=mesh.vert[i].T().P().Y();
            mesh.vert[i].P().Z()=0;
        }
    }

    static void GetMeshFromPatch(MeshType &mesh,
                                 const std::vector<size_t> &Partition,
                                 MeshType &PatchMesh,
                                 bool InternalCuts)
    {
        //size_t t0=clock();

        //std::cout<<"Test C there are :"<<mesh.vert[0].FramePos.size()<<" frames"<<std::endl;


        PatchMesh.Clear();
        PatchMesh.face.reserve(Partition.size());
        PatchMesh.vert.reserve(Partition.size());

        //vcg::tri::UnMarkAll(mesh);
        std::map<size_t,size_t> GlobalToLocal;
        std::vector<size_t> GlobalV;

        std::vector<std::vector<int> > LocalFaceV;
        LocalFaceV.resize(Partition.size(),std::vector<int>(3,-1));
        for (size_t i=0;i<Partition.size();i++)
        {
            size_t IndexF=Partition[i];
            //std::vector<int> FaceV(3,-1);
            for (int j=0;j<mesh.face[IndexF].VN();j++)
            {
                VertexType *v=mesh.face[IndexF].V(j);
                size_t IndexV=vcg::tri::Index(mesh,v);
                //allocate in case it is not already there
                if (GlobalToLocal.count(IndexV)==0)
                {
                    GlobalV.push_back(IndexV);
                    GlobalToLocal[IndexV]=GlobalV.size()-1;
                    LocalFaceV[i][j]=(GlobalV.size()-1);
                }else
                {
                    LocalFaceV[i][j]=GlobalToLocal[IndexV];
                }
            }
            assert(LocalFaceV[i][0]>=0);
            assert(LocalFaceV[i][0]<(int)GlobalV.size());
            assert(LocalFaceV[i][1]>=0);
            assert(LocalFaceV[i][1]<(int)GlobalV.size());
            assert(LocalFaceV[i][2]>=0);
            assert(LocalFaceV[i][2]<(int)GlobalV.size());
        }

        //allocate vertices
        vcg::tri::Allocator<MeshType>::AddVertices(PatchMesh,GlobalV.size());
        assert(PatchMesh.vert.size()==GlobalV.size());

        //copy values from vertices
        for (size_t i=0;i<PatchMesh.vert.size();i++)
        {
            size_t IndexV=GlobalV[i];
            VertexType *OrigV=&mesh.vert[IndexV];

            PatchMesh.vert[i].ImportData(*OrigV);

        }

        //then add faces
        vcg::tri::Allocator<MeshType>::AddFaces(PatchMesh,LocalFaceV.size());
        assert(PatchMesh.face.size()==LocalFaceV.size());
        assert(PatchMesh.face.size()==Partition.size());

        //add vertices and copy values from original ace
        for (size_t i=0;i<PatchMesh.face.size();i++)
        {
            VertexType *v0=&PatchMesh.vert[LocalFaceV[i][0]];
            VertexType *v1=&PatchMesh.vert[LocalFaceV[i][1]];
            VertexType *v2=&PatchMesh.vert[LocalFaceV[i][2]];

            PatchMesh.face[i].V(0)=v0;
            PatchMesh.face[i].V(1)=v1;
            PatchMesh.face[i].V(2)=v2;

            FaceType *f=&mesh.face[Partition[i]];
            PatchMesh.face[i].ImportData(*f);
        }
        PatchMesh.UpdateAttributes();


        //std::cout<<"Test F there are :"<<PatchMesh.vert[0].FramePos.size()<<" frames"<<std::endl;

        if (InternalCuts)
        {
            //bool HasToSplit=false;
            //copy the edge sel

            for (size_t i=0;i<Partition.size();i++)
            {
                size_t IndexF=Partition[i];
                for (int j=0;j<mesh.face[IndexF].VN();j++)
                {
                    if (!mesh.face[IndexF].IsFaceEdgeS(j))continue;
                    PatchMesh.face[i].SetFaceEdgeS(j);
                }
            }
            //std::cout<<"Splitting Internal Cuts"<<std::endl;
            std::map<CoordType,ScalarType> VertQ;
#ifdef MULTI_FRAME
            std::map<CoordType,std::vector<CoordType> > PosFrameMap;
#endif
            for (size_t i=0;i<PatchMesh.vert.size();i++)
            {
                VertQ[PatchMesh.vert[i].P()]=PatchMesh.vert[i].Q();

#ifdef MULTI_FRAME
                PosFrameMap[PatchMesh.vert[i].P()]=PatchMesh.vert[i].FramePos;
#endif
            }
            VertSplitter<MeshType>::SplitAlongEdgeSel(PatchMesh);
            for (size_t i=0;i<PatchMesh.vert.size();i++)
            {
                assert(VertQ.count(PatchMesh.vert[i].P())>0);
                PatchMesh.vert[i].Q()=VertQ[PatchMesh.vert[i].P()];
#ifdef MULTI_FRAME
                PatchMesh.vert[i].FramePos=PosFrameMap[PatchMesh.vert[i].P()];
#endif
            }
            vcg::tri::Allocator<MeshType>::CompactEveryVector(PatchMesh);

            for (size_t i=0;i<PatchMesh.face.size();i++)
            {
                FaceType *f=&mesh.face[Partition[i]];
                PatchMesh.face[i].ImportData(*f);
            }

            PatchMesh.UpdateAttributes();
        }
        //std::cout<<"Test G there are :"<<PatchMesh.vert[0].FramePos.size()<<" frames"<<std::endl;

        //        size_t t1=clock();
        //        Time_InitSubPatches0_0+=t1-t0;
    }



    //    static void MeshDifference(MeshType &mesh0,
    //                               MeshType &mesh1)
    //    {
    //        assert(mesh0.vert.size()==mesh1.vert.size());
    //        assert(mesh0.face.size()==mesh1.face.size());
    //        std::vector<CoordType> Pos0;
    //        std::vector<CoordType> Pos1;
    //        for (size_t i=0;i<mesh0.vert.size();i++)
    //            Pos0.push_back(mesh0.vert[i].P());
    //        for (size_t i=0;i<mesh1.vert.size();i++)
    //            Pos1.push_back(mesh1.vert[i].P());

    ////        std::sort(Pos0.begin(),Pos0.end());
    ////        std::sort(Pos1.begin(),Pos1.end());
    ////        assert(Pos0==Pos1);
    //    }

    //    static void SortByQ(MeshType &mesh)
    //    {
    //        std::vector<std::pair<int,int> > QVertx;
    //        for (size_t i=0;i<mesh.vert.size();i++)
    //            QVertx.push_back(std::pair<int,int>(mesh.vert[i].Q(),i));

    //        std::vector<std::pair<int,int> > QFaces;
    //        for (size_t i=0;i<mesh.face.size();i++)
    //            QFaces.push_back(std::pair<int,int>(mesh.face[i].Q(),i));

    //        std::sort(QVertx.begin(),QVertx.end());
    //        std::sort(QFaces.begin(),QFaces.end());

    //        MeshType Mnew;
    //        for (size_t i=0;i<QVertx.size();i++)
    //        {

    //        }
    ////        assert(mesh0.vert.size()==mesh1.vert.size());
    ////        assert(mesh0.face.size()==mesh1.face.size());
    ////        std::vector<CoordType> Pos0;
    ////        std::vector<CoordType> Pos1;
    ////        for (size_t i=0;i<mesh0.vert.size();i++)
    ////            Pos0.push_back(mesh0.vert[i].P());
    ////        for (size_t i=0;i<mesh1.vert.size();i++)
    ////            Pos1.push_back(mesh1.vert[i].P());

    ////        std::sort(Pos0.begin(),Pos0.end());
    ////        std::sort(Pos1.begin(),Pos1.end());
    ////        assert(Pos0==Pos1);
    //    }

    //    static void GetMeshFromPatch(MeshType &mesh,
    //                                 const std::vector<size_t> &Partition,
    //                                 MeshType &PatchMesh)
    //    {
    ////        size_t t0=clock();
    ////        //        for (size_t i=0;i<mesh.vert.size();i++)
    ////        //            mesh.vert[i].Q()=i;
    ////        //        for (size_t i=0;i<mesh.face.size();i++)
    ////        //            mesh.face[i].Q()=i;

    ////        vcg::tri::UpdateSelection<MeshType>::Clear(mesh);

    ////        size_t ExpNumF=0;
    ////        for (size_t i=0;i<Partition.size();i++)
    ////        {
    ////            mesh.face[Partition[i]].SetS();
    ////            ExpNumF++;
    ////        }

    ////        //        size_t ExpNumV=vcg::tri::UpdateSelection<MeshType>::VertexFromFaceLoose(mesh);

    ////        size_t t1=clock();

    ////        PatchMesh.Clear();
    ////        //        PatchMesh.face.reserve(ExpNumF);
    ////        //        PatchMesh.vert.reserve(ExpNumV);
    ////        vcg::tri::Append<MeshType,MeshType>::Mesh(PatchMesh,mesh,true);
    ////        PatchMesh.UpdateAttributes();

    ////        size_t t2=clock();

    ////        Time_InitSubPatches0_0+=t1-t0;
    ////        Time_InitSubPatches0_1+=t2-t1;

    //////        vcg::tri::UpdateSelection<MeshType>::Clear(PatchMesh);
    //////        vcg::tri::UpdateSelection<MeshType>::Clear(mesh);

    ////        MeshType PatchMesh2;
    ////        GetMeshFromPatch2(mesh,Partition,PatchMesh2);
    ////        MeshDifference(PatchMesh,PatchMesh2);

    //     }

    static void GetMeshFromPatch(MeshType &mesh,
                                 const size_t &IndexPatch,
                                 const std::vector<std::vector<size_t> > &Partitions,
                                 MeshType &PatchMesh,
                                 bool InternalCuts)
    {
        GetMeshFromPatch(mesh,Partitions[IndexPatch],PatchMesh,InternalCuts);
    }


    static void GetUVOutline(MeshType &testM,
                             std::vector<vcg::Point2d > &Poly2D,
                             //vcg::Box2d &box2D,
                             typename MeshType::ScalarType borders)
    {
        vcg::tri::UpdateTopology<MeshType>::FaceFace(testM);
        vcg::tri::UpdateFlags<MeshType>::VertexBorderFromFaceAdj(testM);
        //box2D.SetNull();
        vcg::Box2d PolyBox;
        for (size_t i=0;i<testM.vert.size();i++)
        {
            if (!testM.vert[i].IsB())continue;
            Poly2D.push_back(testM.vert[i].T().P());
            PolyBox.Add(testM.vert[i].T().P());
            //box2D.Add(testM.vert[i].T().P());
        }
        if (borders==0)return;

        vcg::Box2d PolyBoxTarget=PolyBox;
        PolyBoxTarget.Offset(borders);

        ScalarType ScaleVal=PolyBoxTarget.Diag()/PolyBox.Diag();
        for (size_t i=0;i<Poly2D.size();i++)
        {
            //translate to the center
            Poly2D[i]-=PolyBox.Center();
            //then scale to match borders
            Poly2D[i]*=ScaleVal;
            //translate back
            Poly2D[i]+=PolyBox.Center();
        }

        //box2D.Offset(borders);
    }


    static void ArrangeUVPatches(std::vector<MeshType*> &ParamPatches,
                                 typename MeshType::ScalarType borders=0,
                                 bool allowRot=false)
    {
        typedef typename MeshType::ScalarType ScalarType;
        std::vector<std::vector<vcg::Point2<ScalarType> > > ParaOutlines;
        ParaOutlines.resize(ParamPatches.size());

        //vcg::Box2d box2D,totalBox;
        ScalarType AreaUV0=0;
        for (size_t i=0;i<ParamPatches.size();i++)
        {
            //GetUVOutline(*ParamPatches[i],ParaOutlines[i],box2D,borders);
            GetUVOutline(*ParamPatches[i],ParaOutlines[i],borders);
            AreaUV0+=vcg::tri::UV_Utils<MeshType>::PerVertUVArea(*ParamPatches[i]);//(box2D.DimX()*box2D.DimY());
            //totalBox.Add(box2D);
        }

        int EdgeSize=floor(1000);//sqrt(AreaUV)+0.5);
        EdgeSize=std::max(EdgeSize,1);
        vcg::Point2i siz(EdgeSize*2,EdgeSize*2);
        vcg::Point2<ScalarType> coveredA;
        std::vector<vcg::Similarity2<ScalarType> > trVec;
        if (allowRot)
            vcg::PolyPacker<ScalarType>::PackAsObjectOrientedRect(ParaOutlines,siz,trVec,coveredA,borders);
        else
            vcg::PolyPacker<ScalarType>::PackAsAxisAlignedRect(ParaOutlines,siz,trVec,coveredA);//,borders);

        //find the final position
        ScalarType AreaUV1=0;
        for (size_t i=0;i<ParamPatches.size();i++)
        {
            for (size_t j=0;j<(*ParamPatches[i]).vert.size();j++)
            {
                vcg::Point2<ScalarType> UVVert=(*ParamPatches[i]).vert[j].T().P();
                (*ParamPatches[i]).vert[j].T().P()=trVec[i]*UVVert;
            }
            AreaUV1+=vcg::tri::UV_Utils<MeshType>::PerVertUVArea(*ParamPatches[i]);
        }

        ScalarType ratio=sqrt(AreaUV0/AreaUV1);

        //rescale it back
        ScalarType AreaUV2=0;
        for (size_t i=0;i<ParamPatches.size();i++)
        {
            for (size_t j=0;j<(*ParamPatches[i]).vert.size();j++)
                (*ParamPatches[i]).vert[j].T().P()*=ratio;

            AreaUV2+=vcg::tri::UV_Utils<MeshType>::PerVertUVArea(*ParamPatches[i]);
        }

        //        std::cout<<"Area UV0:"<<AreaUV0<<std::endl;
        //        std::cout<<"Area UV2:"<<AreaUV2<<std::endl;

    }



    static void SelectMeshPatchBorders(MeshType &mesh,
                                       const std::vector<int>  &FacePatches,
                                       bool SetF=false)
    {
        //std::set<std::pair<size_t,size_t> > BorderPatches;
        if (!SetF)
            vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(mesh);
        else
            vcg::tri::UpdateFlags<MeshType>::FaceSetF(mesh);

        assert(FacePatches.size()==mesh.face.size());
        //first add borders
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (mesh.face[i].IsB(j))
                {
                    if (!SetF)
                        mesh.face[i].SetFaceEdgeS(j);
                    else
                        mesh.face[i].ClearF(j);
                }
                else
                {
                    size_t FOpp=vcg::tri::Index(mesh,mesh.face[i].FFp(j));
                    assert(FOpp!=i);
                    if (FacePatches[i]<0)assert(FacePatches[i]==-1);
                    //assert(FacePatches[i]>=0);

                    if (FacePatches[i]!=FacePatches[FOpp])
                    {
                        if (!SetF)
                            mesh.face[i].SetFaceEdgeS(j);
                        else
                            mesh.face[i].ClearF(j);
                    }
                }
            }
    }


    static void SelectMeshPatchBorders(MeshType &mesh,
                                       const std::vector<std::vector<size_t> >  &PatchFaces,
                                       bool SetF=false)
    {
        std::vector<int > FacePartitions;
        DerivePerFacePartition(mesh,PatchFaces,FacePartitions);
        SelectMeshPatchBorders(mesh,FacePartitions,SetF);
    }


    static void SelectVertOnMeshPatchBorders(MeshType &mesh,const std::vector<int>  &FacePatches)
    {
        //std::set<std::pair<size_t,size_t> > BorderPatches;
        vcg::tri::UpdateFlags<MeshType>::VertexClearS(mesh);
        assert(FacePatches.size()==mesh.face.size());
        //first add borders
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (mesh.face[i].IsB(j))
                {
                    mesh.face[i].V0(j)->SetS();
                    mesh.face[i].V1(j)->SetS();
                }
                else
                {
                    size_t FOpp=vcg::tri::Index(mesh,mesh.face[i].FFp(j));
                    assert(FOpp!=i);
                    //assert(FacePatches[i]>=0);

                    if (FacePatches[i]!=FacePatches[FOpp])
                    {
                        mesh.face[i].V0(j)->SetS();
                        mesh.face[i].V1(j)->SetS();
                    }
                }
            }
    }


    static void DerivePerFacePartition(const MeshType &mesh,
                                       const std::vector<std::vector<size_t> > &Partitions,
                                       std::vector<int > &FacePartitions)
    {
        FacePartitions.clear();
        FacePartitions.resize(mesh.face.size(),-1);

        for (size_t i=0;i<Partitions.size();i++)
            for (size_t j=0;j<Partitions[i].size();j++)
            {
                size_t IndexF=Partitions[i][j];
                FacePartitions[IndexF]=i;
            }
    }

    static void DerivePerPartitionFaces(const std::vector<int> &FacePartitions,
                                        std::vector<std::vector<size_t> > &PartitionFaces)
    {
        PartitionFaces.clear();
        //get the maximum partition index
        size_t maxP=(*std::max_element(FacePartitions.begin(),FacePartitions.end()));
        PartitionFaces.resize(maxP+1);

        for (size_t i=0;i<FacePartitions.size();i++)
        {
            size_t IndexP=FacePartitions[i];
            PartitionFaces[IndexP].push_back(i);
        }
    }

    static void SelectVertOnMeshPatchBorders(MeshType &mesh,const std::vector<std::vector<int> >  &PatchFaces)
    {
        std::vector<int> FacePartitions;
        DerivePerFacePartition(mesh,PatchFaces,FacePartitions);
        SelectVertOnMeshPatchBorders(mesh,FacePartitions);
    }


    static void SelectTJunctionVertices(MeshType &mesh,const std::vector<int>  &FacePatches)
    {
        std::vector<std::set<int> > VertPatch(mesh.vert.size());
        for (size_t i=0;i<FacePatches.size();i++)
            for (size_t j=0;j<3;j++)
            {
                size_t IndexV=vcg::tri::Index(mesh,mesh.face[i].V(j));
                VertPatch[IndexV].insert(FacePatches[i]);
            }

        for (size_t i=0;i<mesh.vert.size();i++)
            if (VertPatch[i].size()==3)
                mesh.vert[i].SetS();
    }



    static int SelectFolds(MeshType &mesh,
                           int dilateStep=0,
                           typename MeshType::ScalarType MinQ=0.2)
    {
        //typedef typename MeshType::FaceType FaceType;
        typedef typename MeshType::CoordType CoordType;
        typedef typename MeshType::ScalarType ScalarType;

        size_t numF=0;
        vcg::tri::UpdateSelection<MeshType>::Clear(mesh);
        //first select faces with bad aspect ratio
        for (size_t i=0;i<mesh.face.size();i++)
        {
            CoordType P0=mesh.face[i].P(0);
            CoordType P1=mesh.face[i].P(1);
            CoordType P2=mesh.face[i].P(2);
            ScalarType QFace=vcg::QualityRadii(P0,P1,P2);
            if (QFace<MinQ){mesh.face[i].SetS();numF++;}
        }

        //save smooth normals
        std::vector<CoordType> SmoothNorms;
        std::vector<CoordType> SmoothPos;
        vcg::tri::UpdateNormal<MeshType>::PerFaceNormalized(mesh);
        for (size_t i=0;i<mesh.vert.size();i++)
            SmoothPos.push_back(mesh.vert[i].P());
        for (size_t i=0;i<mesh.face.size();i++)
            SmoothNorms.push_back(mesh.face[i].N());

        //compute rest normals
        for (size_t i=0;i<mesh.vert.size();i++)
            mesh.vert[i].P()=mesh.vert[i].RPos;
        std::vector<CoordType> RestNorms;
        vcg::tri::UpdateNormal<MeshType>::PerFaceNormalized(mesh);
        for (size_t i=0;i<mesh.vert.size();i++)
            RestNorms.push_back(mesh.face[i].N());

        for (size_t i=0;i<RestNorms.size();i++)
            if (RestNorms[i]*SmoothNorms[i]<0){mesh.face[i].SetS();numF++;}

        vcg::tri::UpdateSelection<MeshType>::VertexFromFaceLoose(mesh);

        //restore the position
        for (size_t i=0;i<mesh.vert.size();i++)
            mesh.vert[i].P()=SmoothPos[i];
        vcg::tri::UpdateNormal<MeshType>::PerFaceNormalized(mesh);

        for (int i=0;i<dilateStep;i++)
        {
            vcg::tri::UpdateSelection<MeshType>::FaceFromVertexLoose(mesh);
            vcg::tri::UpdateSelection<MeshType>::VertexFromFaceLoose(mesh);
        }
        return numF;
    }

    static void SolveFolds(MeshType &mesh,
                           typename MeshType::ScalarType MinQ)
    {
        size_t dilateS=0;
        size_t numF0=SelectFolds(mesh);
        size_t numF1=numF0;
        if (numF0==0)return;
        do{
            std::cout<<"There are "<<numF0<<" folds"<<std::endl;

            //unselect the boundaries
            for (size_t i=0;i<mesh.vert.size();i++)
                if (mesh.vert[i].IsB())mesh.vert[i].ClearS();

            vcg::tri::Smooth<MeshType>::VertexCoordLaplacian(mesh,10,true);
            dilateS++;
            numF1=SelectFolds(mesh,dilateS,MinQ);
            if (numF1==0)return;
        }while ((numF0>=numF1)&&(dilateS>3));
    }


    static void LaplacianPos(MeshType &mesh,
                             std::vector<typename MeshType::CoordType> &AvPos,
                             bool OnlySContribute)
    {

        AvPos=std::vector<CoordType>(mesh.vert.size(),CoordType(0,0,0));
        std::vector<size_t> NumDiv(mesh.vert.size(),0);

        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
                CoordType Pos0=mesh.vert[IndexV0].P();
                CoordType Pos1=mesh.vert[IndexV1].P();
                bool IsSel0=mesh.vert[IndexV0].IsS();
                bool IsSel1=mesh.vert[IndexV1].IsS();
                bool addContr0=(!OnlySContribute)||((OnlySContribute)&&(IsSel0));
                bool addContr1=(!OnlySContribute)||((OnlySContribute)&&(IsSel1));
                if (addContr0)
                {
                    AvPos[IndexV1]+=Pos0;
                    NumDiv[IndexV1]++;
                }
                if (addContr1)
                {
                    AvPos[IndexV0]+=Pos1;
                    NumDiv[IndexV0]++;
                }
            }

        for (size_t i=0;i<mesh.vert.size();i++)
        {
            //no contributes
            if (NumDiv[i]==0)continue;
            AvPos[i]/=(ScalarType)NumDiv[i];
        }
    }


//    static void LaplacianInternalStep(MeshType &mesh,//const std::vector<int>  &FacePatches,
//                                      typename MeshType::ScalarType Damp)//,bool FixV)
//    {
//        //SelectVertOnMeshPatchBorders(mesh,FacePatches);
//        std::vector<typename MeshType::CoordType> AvPos;
//        LaplacianPos(mesh,AvPos,false);
//        for (size_t i=0;i<mesh.vert.size();i++)
//        {
//            if (mesh.vert[i].IsS())continue;
//            //if (mesh.vert[i].IsV() && FixV)continue;
//            mesh.vert[i].P()=mesh.vert[i].P()*Damp+AvPos[i]*(1-Damp);
//        }
//    }
//    static void SetVCornersFromEdgeSel()
//    {

//        for (size_t i=0;i<mesh.face.size();i++)
//            for (size_t j=0;j<3;j++)
//            {
//                if (!mesh.face[i].IsFaceEdgeS(j))continue;
//                size_t VInd0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
//                size_t VInd1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
//                mesh.vert[VInd0].SetV();
//                mesh.vert[VInd1].SetV();
//            }
//    }

    static void LaplacianNonSelectedEdgeStep(MeshType &mesh,//const std::vector<int>  &FacePatches,
                                         typename MeshType::ScalarType Damp)//,bool FixV)
    {
        //SelectVertOnMeshPatchBorders(mesh,FacePatches);
        std::vector<typename MeshType::CoordType> AvPos;
        LaplacianPos(mesh,AvPos,false);
        //mark the one not selected
        vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!mesh.face[i].IsFaceEdgeS(j))continue;
                size_t VInd0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                size_t VInd1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
                mesh.vert[VInd0].SetV();
                mesh.vert[VInd1].SetV();
            }

        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (mesh.vert[i].IsV())continue;
            //if (mesh.vert[i].IsV() && FixV)continue;
            mesh.vert[i].P()=mesh.vert[i].P()*Damp+AvPos[i]*(1-Damp);
        }
    }

//    static void LaplacianSelectedEdgeStep(MeshType &mesh,
//                                    //const std::vector<int>  &FacePatches,
//                                      typename MeshType::ScalarType Damp)//,
//    //bool FixV)
//    {

//        //save previous selection

//        typedef typename MeshType::CoordType CoordType;
//        typedef typename MeshType::ScalarType ScalarType;
//        //SelectMeshPatchBorders(mesh,FacePatches);

//        std::vector<size_t> NumDiv(mesh.vert.size(),0);
//        std::vector<typename MeshType::CoordType> AvPos(mesh.vert.size(),
//                                                        CoordType(0,0,0));

//        //smooth path
//        for (size_t i=0;i<mesh.face.size();i++)
//            for (size_t j=0;j<3;j++)
//            {
//                if (!mesh.face[i].IsFaceEdgeS(j))continue;
//                size_t VInd0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
//                size_t VInd1=vcg::tri::Index(mesh,mesh.face[i].V1(j));

//                CoordType Pos0=mesh.vert[VInd0].P();
//                CoordType Pos1=mesh.vert[VInd1].P();

//                AvPos[VInd0]+=Pos1;
//                NumDiv[VInd0]++;
//                AvPos[VInd1]+=Pos0;
//                NumDiv[VInd1]++;
//            }

//        for (size_t i=0;i<mesh.vert.size();i++)
//        {
//            //no contributes
//            if (NumDiv[i]==0)continue;
//            CoordType TargetPos=AvPos[i]/(ScalarType)NumDiv[i];
//            if (mesh.vert[i].IsB())continue;
//            //if (mesh.vert[i].IsV()&&FixV)continue;
//            mesh.vert[i].P()=mesh.vert[i].P()*Damp+TargetPos*(1-Damp);
//        }
//        //vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(mesh);
//    }


    static void LaplacianSelectedEdgeStep(MeshType &mesh,
                                    //const std::vector<int>  &FacePatches,
                                      typename MeshType::ScalarType Damp)//,
    //bool FixV)
    {

        //save previous selection

        typedef typename MeshType::CoordType CoordType;
        typedef typename MeshType::ScalarType ScalarType;
        //SelectMeshPatchBorders(mesh,FacePatches);

        std::vector<size_t> NumDiv(mesh.vert.size(),0);
        std::vector<typename MeshType::CoordType> AvPos(mesh.vert.size(),
                                                        CoordType(0,0,0));

        //smooth path
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!mesh.face[i].IsFaceEdgeS(j))continue;
                size_t VInd0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                size_t VInd1=vcg::tri::Index(mesh,mesh.face[i].V1(j));

                bool AddV=false;
                AddV|=(vcg::face::IsBorder(mesh.face[i],j));
                AddV|=(VInd0<VInd1);

                if (!AddV)continue;
                CoordType Pos0=mesh.vert[VInd0].P();
                CoordType Pos1=mesh.vert[VInd1].P();

                AvPos[VInd0]+=Pos1;
                NumDiv[VInd0]++;
                AvPos[VInd1]+=Pos0;
                NumDiv[VInd1]++;
            }

        for (size_t i=0;i<mesh.vert.size();i++)
        {
            //no contributes
            //if (NumDiv[i]==0)continue;
            if (NumDiv[i]!=2)continue;
            CoordType TargetPos=AvPos[i]/(ScalarType)NumDiv[i];
            if (mesh.vert[i].IsB())continue;
            //if (mesh.vert[i].IsV()&&FixV)continue;
            mesh.vert[i].P()=mesh.vert[i].P()*Damp+TargetPos*(1-Damp);
        }
        //vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(mesh);
    }

    static void ReprojectOn(MeshType &mesh,MeshType &target,
                            vcg::GridStaticPtr<typename MeshType::FaceType,typename MeshType::ScalarType> &Gr)
    {

        //reproject everything
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            //fixed
            CoordType Pos;
            ScalarType MinD;
            vcg::tri::GetClosestFaceBase(target,Gr,
                                         mesh.vert[i].P(),
                                         mesh.bbox.Diag(),MinD,Pos);
            mesh.vert[i].P()=Pos;
        }
    }

    static void SaveEdgeSel(MeshType &mesh,std::vector<std::vector<bool> > &EdgeSel)
    {
        assert(EdgeSel.size()==mesh.face.size());
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!mesh.face[i].IsFaceEdgeS(j))
                    EdgeSel[i][j]=false;
                else
                    EdgeSel[i][j]=true;
            }
    }



    static void RestoreEdgeSel(MeshType &mesh,
                               const std::vector<std::vector<bool> > &EdgeSel)
    {
        vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(mesh);
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!EdgeSel[i][j])continue;
                mesh.face[i].SetFaceEdgeS(j);
            }
    }

//    static void SmoothMeshPatches(MeshType &mesh,
//                                  const std::vector<int>  &FacePatches,
//                                  bool UsePredefinedSelE,
//                                  size_t Steps,
//                                  typename MeshType::ScalarType Damp,
//                                  typename MeshType::ScalarType MinQ)
//    {

//        std::vector<std::vector<bool> > EdgeSel;

//        //select borders
//        if (!UsePredefinedSelE)
//        {
//            EdgeSel.resize(mesh.face.size(),std::vector<bool>(3,false));
//            SaveEdgeSel(mesh,EdgeSel);
//            SelectMeshPatchBorders(mesh,FacePatches);
//        }

//        SelectVertOnMeshPatchBorders(mesh,FacePatches);

//        //init reprojection grid
//        MeshType TargetMesh;
//        vcg::tri::Append<MeshType,MeshType>::Mesh(TargetMesh,mesh);
//        TargetMesh.UpdateAttributes();
//        vcg::GridStaticPtr<FaceType,ScalarType> Gr;
//        Gr.Set(TargetMesh.face.begin(),TargetMesh.face.end());

//        //then for each smooth step
//        for (size_t s=0;s<Steps;s++)
//        {
//            //save old position in case quality check is needed
//            std::vector<CoordType> OldPos;
//            if (MinQ>0)
//            {
//                for (size_t i=0;i<mesh.vert.size();i++)
//                    OldPos.push_back(mesh.vert[i].P());
//            }

//            //save old normals
//            vcg::tri::UpdateNormal<MeshType>::PerFaceNormalized(mesh);
//            std::vector<CoordType> OldNorm;
//            if (MinQ>0)
//            {
//                for (size_t i=0;i<mesh.face.size();i++)
//                    OldNorm.push_back(mesh.face[i].N());
//            }

//            //PERFORM SMOOTHING
//            LaplacianBorderStep(mesh,Damp);//,true);
//            ReprojectOn(mesh,TargetMesh,Gr);
//            LaplacianInternalStep(mesh,Damp);//,true);
//            ReprojectOn(mesh,TargetMesh,Gr);

//            //no quality check we are fine!
//            if (MinQ<=0)continue;

//            //check each face
//            for (size_t i=0;i<mesh.face.size();i++)
//            {
//                //find the one that has 2 boder edges
//                int indexE=-1;
//                for (size_t j=0;j<3;j++)
//                {
//                    if (mesh.face[i].IsFaceEdgeS(j) &&
//                            mesh.face[i].IsFaceEdgeS((j+1)%3))
//                        indexE=j;
//                }
//                if (indexE==-1)continue;
//                size_t IndexV=vcg::tri::Index(mesh,mesh.face[i].V(indexE));

//                //if border not change
//                if (mesh.vert[IndexV].IsB())continue;

//                //check quality of the face
//                bool IsOk=true;
//                CoordType P0=mesh.face[i].P(0);
//                CoordType P1=mesh.face[i].P(1);
//                CoordType P2=mesh.face[i].P(2);
//                CoordType NewNormF=vcg::Normal(P0,P1,P2);
//                NewNormF.Normalize();
//                CoordType OldNormF=OldNorm[i];
//                if ((NewNormF*OldNormF)<0)
//                    IsOk=false;
//                ScalarType QFace=vcg::QualityRadii(P0,P1,P2);
//                if (QFace<MinQ)
//                    IsOk=false;

//                //restore if not ok
//                if (!IsOk)
//                {
//                    mesh.vert[IndexV].P()=OldPos[IndexV];
//                }
//            }
//        }
//        vcg::tri::UpdateSelection<MeshType>::Clear(mesh);
//        if (MinQ>0)
//            SolveFolds(mesh,MinQ);

//        if (!UsePredefinedSelE)
//            RestoreEdgeSel(mesh,EdgeSel);
//    }

    static void SmoothMeshPatchesFromEdgeSel(MeshType &mesh,
                                  size_t Steps,
                                  typename MeshType::ScalarType Damp,
                                  typename MeshType::ScalarType MinQ)
    {

        //init reprojection grid
        MeshType TargetMesh;
        vcg::tri::Append<MeshType,MeshType>::Mesh(TargetMesh,mesh);
        TargetMesh.UpdateAttributes();
        vcg::GridStaticPtr<FaceType,ScalarType> Gr;
        Gr.Set(TargetMesh.face.begin(),TargetMesh.face.end());

        //then for each smooth step
        for (size_t s=0;s<Steps;s++)
        {
            //save old position in case quality check is needed
            std::vector<CoordType> OldPos;
            if (MinQ>0)
            {
                for (size_t i=0;i<mesh.vert.size();i++)
                    OldPos.push_back(mesh.vert[i].P());
            }

            //save old normals
            vcg::tri::UpdateNormal<MeshType>::PerFaceNormalized(mesh);
            std::vector<CoordType> OldNorm;
            if (MinQ>0)
            {
                for (size_t i=0;i<mesh.face.size();i++)
                    OldNorm.push_back(mesh.face[i].N());
            }

            //PERFORM SMOOTHING
            LaplacianSelectedEdgeStep(mesh,Damp);//,true);
            ReprojectOn(mesh,TargetMesh,Gr);
            LaplacianNonSelectedEdgeStep(mesh,Damp);//,true);
            ReprojectOn(mesh,TargetMesh,Gr);

            //no quality check we are fine!
            if (MinQ<=0)continue;

            //check each face
            for (size_t i=0;i<mesh.face.size();i++)
            {
                //find the one that has 2 boder edges
                int indexE=-1;
                for (size_t j=0;j<3;j++)
                {
                    if (mesh.face[i].IsFaceEdgeS(j) &&
                            mesh.face[i].IsFaceEdgeS((j+1)%3))
                        indexE=j;
                }
                if (indexE==-1)continue;
                size_t IndexV=vcg::tri::Index(mesh,mesh.face[i].V(indexE));

                //if border not change
                if (mesh.vert[IndexV].IsB())continue;

                //check quality of the face
                bool IsOk=true;
                CoordType P0=mesh.face[i].P(0);
                CoordType P1=mesh.face[i].P(1);
                CoordType P2=mesh.face[i].P(2);
                CoordType NewNormF=vcg::Normal(P0,P1,P2);
                NewNormF.Normalize();
                CoordType OldNormF=OldNorm[i];
                if ((NewNormF*OldNormF)<0)
                    IsOk=false;
                ScalarType QFace=vcg::QualityRadii(P0,P1,P2);
                if (QFace<MinQ)
                    IsOk=false;

                //restore if not ok
                if (!IsOk)
                {
                    mesh.vert[IndexV].P()=OldPos[IndexV];
                }
            }
        }
        vcg::tri::UpdateSelection<MeshType>::Clear(mesh);
        if (MinQ>0)
            SolveFolds(mesh,MinQ);
    }

    static void PatchSideLenght(MeshType &mesh,
                                const std::vector<size_t>  &PatchFaces,
                                const std::vector<size_t>  &PatchCorners,
                                std::vector<typename MeshType::ScalarType>  &CurvedL,
                                std::vector<typename MeshType::ScalarType>  &EuclL,
                                //std::unordered_map<std::pair<size_t,size_t>,typename MeshType::ScalarType>  &EdgeMap)
                                std::map<std::pair<size_t,size_t>,typename MeshType::ScalarType> &EdgeMap)
    {

        //    int t0=clock();
        //then get the border sequence
        std::vector<std::vector<size_t> > BorderSeq;
        BorderSeq.resize(PatchFaces.size());

        if (PatchCorners.size()==0)return;
        SelectVertices(mesh,PatchCorners);
        //DerivePatchBorderSeq(mesh,PatchFaces,PatchCorners,BorderSeq);
        DeriveFauxSeq(mesh,PatchFaces,BorderSeq);//(mesh,PatchFaces,BorderSeq);
        //    std::cout<<"Size 0 "<<PatchCorners.size()<<std::endl;
        //    std::cout<<"Size 1 "<<BorderSeq.size()<<std::endl;
        //    std::cout<<"Size 2 "<<BorderSeq[0].size()<<std::endl;
        //    assert(BorderSeq.size()==PatchCorners.size());

        //    int t1=clock();
        for (size_t j=0;j<BorderSeq.size();j++)
        {
            //CurvedL.push_back(FieldLenght(mesh,BorderSeq[j],EdgeMap));
            CurvedL.push_back(FieldLenght(BorderSeq[j],EdgeMap));
            CoordType Pos0=mesh.vert[BorderSeq[j][0]].P();
            CoordType Pos1=mesh.vert[BorderSeq[j].back()].P();
            EuclL.push_back((Pos1-Pos0).Norm());
        }
        //    int t2=clock();
        //        std::cout<<"** Timing Patch Length **"<<std::endl;
        //        std::cout<<"Time Derive Border Seq "<<t1-t0<<std::endl;
        //        std::cout<<"Time Derive Border Len "<<t2-t1<<std::endl;
    }


    static void PatchesSideLenght(MeshType &mesh,
                                  const std::vector<std::vector<size_t> > &PatchFaces,
                                  const std::vector<std::vector<size_t> > &PatchCorners,
                                  std::vector<std::vector<typename MeshType::ScalarType> > &CurvedL,
                                  std::vector<std::vector<typename MeshType::ScalarType> > &EuclL,
                                  //std::unordered_map<std::pair<size_t,size_t>,typename MeshType::ScalarType>  &EdgeMap)
                                  std::map<std::pair<size_t,size_t>,typename MeshType::ScalarType> &EdgeMap)
    {

        CurvedL.clear();
        EuclL.clear();
        CurvedL.resize(PatchFaces.size());
        EuclL.resize(PatchFaces.size());

        //    SetFacePartitionOnFaceQ(mesh,PatchFaces);
        //    SelectCorners(mesh,PatchCorners);
        //    SelectMeshPatchBorders(mesh,PatchFaces);

        //    std::vector<std::vector<bool> > EdgeSel;
        //    SaveEdgeSel(mesh,EdgeSel);
        //    SelectMeshPatchBorders(mesh,PatchFaces);
        //    RestoreEdgeSel(mesh,EdgeSel);

        SelectMeshPatchBorders(mesh,PatchFaces,true);
        SelectVertices(mesh,PatchCorners);

        for (size_t i=0;i<PatchFaces.size();i++)
        {
            PatchSideLenght(mesh,PatchFaces[i],PatchCorners[i],CurvedL[i],EuclL[i],EdgeMap);
        }

    }


    static void PatchesLenghtRatios(MeshType &mesh,
                                    const std::vector<std::vector<size_t> > &PatchFaces,
                                    const std::vector<std::vector<size_t> > &PatchCorners,
                                    std::vector<typename MeshType::ScalarType> &Variance,
                                    std::vector<typename MeshType::ScalarType> &LenghtDist,
                                    //std::unordered_map<std::pair<size_t,size_t>,typename MeshType::ScalarType>  &EdgeMap)
                                    std::map<std::pair<size_t,size_t>,typename MeshType::ScalarType> &EdgeMap)
    {
        std::vector<std::vector<ScalarType> > SideL,EuclL;
        PatchesSideLenght(mesh,PatchFaces,PatchCorners,SideL,EuclL,EdgeMap);
        for (size_t i=0;i<SideL.size();i++)
        {
            ScalarType MinL=(*std::min_element(SideL[i].begin(), SideL[i].end()));
            ScalarType MaxL=(*std::max_element(EuclL[i].begin(), EuclL[i].end()));
            Variance.push_back(MaxL/MinL);

            ScalarType MaxR=0;
            for (size_t j=0;j<SideL[i].size();j++)
                MaxR=std::max(MaxR,SideL[i][j]/EuclL[i][j]);

            LenghtDist.push_back(MaxR);
        }
    }


    static void ColorByVarianceLenght(MeshType &mesh,
                                      const std::vector<std::vector<size_t> > &PatchFaces,
                                      const std::vector<std::vector<size_t> > &PatchCorners,
                                      //std::unordered_map<std::pair<size_t,size_t>,typename MeshType::ScalarType>  &EdgeMap)
                                      std::map<std::pair<size_t,size_t>,typename MeshType::ScalarType> &EdgeMap)
    {
        std::vector<ScalarType> VarianceL,DistL;
        PatchesLenghtRatios(mesh,PatchFaces,PatchCorners,VarianceL,DistL,EdgeMap);
        for (size_t i=0;i<PatchFaces.size();i++)
        {
            for (size_t j=0;j<PatchFaces[i].size();j++)
            {
                size_t IndexF=PatchFaces[i][j];
                mesh.face[IndexF].Q()=VarianceL[i];
            }
        }
        vcg::tri::UpdateColor<MeshType>::PerFaceQualityRamp(mesh);
    }


    static void ColorByDistortionLenght(MeshType &mesh,
                                        const std::vector<std::vector<size_t> > &PatchFaces,
                                        const std::vector<std::vector<size_t> > &PatchCorners,
                                        //std::unordered_map<std::pair<size_t,size_t>,typename MeshType::ScalarType>  &EdgeMap)
                                        std::map<std::pair<size_t,size_t>,typename MeshType::ScalarType> &EdgeMap)
    {
        std::vector<ScalarType> VarianceL,DistL;

        PatchesLenghtRatios(mesh,PatchFaces,PatchCorners,VarianceL,DistL,EdgeMap);
        for (size_t i=0;i<PatchFaces.size();i++)
        {
            for (size_t j=0;j<PatchFaces[i].size();j++)
            {
                size_t IndexF=PatchFaces[i][j];
                mesh.face[IndexF].Q()=DistL[i];
            }
        }
        std::pair<ScalarType, ScalarType> minmax = vcg::tri::Stat<MeshType>::ComputePerFaceQualityMinMax(mesh);
        std::cout<<"Min Dist: "<<minmax.first<<std::endl;
        std::cout<<"Max Dist: "<<minmax.second<<std::endl;
        vcg::tri::UpdateColor<MeshType>::PerFaceQualityRamp(mesh,minmax.second,minmax.first);
    }


    static void ColorByCatmullClarkability(MeshType &mesh,
                                           const std::vector<std::vector<size_t> > &PatchFaces,
                                           const std::vector<std::vector<size_t> > &PatchCorners,
                                           //std::unordered_map<std::pair<size_t,size_t>,typename MeshType::ScalarType>  &EdgeMap,
                                           std::map<std::pair<size_t,size_t>,typename MeshType::ScalarType>  &EdgeMap,
                                           const typename MeshType::ScalarType CClarkability,
                                           const typename MeshType::ScalarType avEdge,
                                           bool SkipValence4)
    {
        std::vector<std::vector<ScalarType> > SideL,EuclL;
        PatchesSideLenght(mesh,PatchFaces,PatchCorners,SideL,EuclL,EdgeMap);
        for (size_t i=0;i<SideL.size();i++)
        {
            //bool CC=IsCatmullClarkable(PatchFaces[i].size(),SideL[i],Thr,SkipValence4);
            //            size_t addedS=AddedSingularities(PatchFaces[i].size(),SideL[i],
            //                                             CClarkability*avEdge,SkipValence4);

            size_t addedS=AddedSingularities(SideL[i],CClarkability*avEdge,SkipValence4);

            for (size_t j=0;j<PatchFaces[i].size();j++)
            {
                size_t IndexF=PatchFaces[i][j];
                if (addedS==0)
                    mesh.face[IndexF].C()=vcg::Color4b::Green;
                else
                    mesh.face[IndexF].C()=vcg::Color4b::Red;
            }
        }
    }


    static void ParametrizePatches(MeshType &mesh,
                                   MeshType &splittedUV,
                                   std::vector<std::vector<size_t> > &PatchFaces,
                                   std::vector<std::vector<size_t> > &PatchCorners,
                                   ParamType PType,
                                   bool FixBorders,
                                   bool ScaleEdges,
                                   bool Arrange=true,
                                   bool InternalCuts=false,
                                   typename MeshType::ScalarType borderArrange=0)
    {
        std::vector<MeshType*> ParamPatches;
        for (size_t i=0;i<PatchFaces.size();i++)
        {
            ParamPatches.push_back(new MeshType());
            std::vector<size_t> SortedCorners;
            ComputeParametrizedSubMesh(mesh,*ParamPatches.back(),
                                       PatchFaces[i],PatchCorners[i],
                                       SortedCorners,PType,
                                       FixBorders,ScaleEdges,
                                       InternalCuts);
        }
        if (Arrange)
            ArrangeUVPatches(ParamPatches,borderArrange);

        splittedUV.Clear();
        for (size_t i=0;i<PatchFaces.size();i++)
        {
            vcg::tri::Append<MeshType,MeshType>::Mesh(splittedUV,*ParamPatches[i]);
            delete(ParamPatches[i]);
        }
        //            if (Save)
        //            {
        SetUVtoPos(splittedUV);
        vcg::tri::io::ExporterPLY<MeshType>::Save(splittedUV,"parametrize.ply",vcg::tri::io::Mask::IOM_FACECOLOR);
        //            }
    }



    static void ColorByUVDistortionFaces(MeshType &mesh,
                                         std::vector<std::vector<size_t> > &PatchFaces,
                                         std::vector<std::vector<size_t> > &PatchCorners,
                                         ParamType PType,
                                         bool FixCorners,
                                         bool FixBorders)
    {
        MeshType splittedUV;
        ParametrizePatches(mesh,splittedUV,PatchFaces,PatchCorners,PType,FixCorners,FixBorders,false);
        std::vector<size_t> OrigIndex;
        for (size_t i=0;i<splittedUV.face.size();i++)
            OrigIndex.push_back(splittedUV.face[i].Q());

        vcg::tri::Distortion<MeshType,false>::SetQasDistorsion(splittedUV,vcg::tri::Distortion<MeshType,false>::ARAPDist);
        for (size_t i=0;i<splittedUV.face.size();i++)
        {
            size_t IndexF=OrigIndex[i];
            mesh.face[IndexF].Q()=splittedUV.face[i].Q();
        }
        vcg::tri::UpdateColor<MeshType>::PerFaceQualityRamp(mesh);
    }


    static size_t RemainingEmitters(MeshType &mesh,
                                    std::vector<size_t> &PatchFaces,
                                    std::vector<size_t> &VerticesNeeds)
    {
        size_t ret=0;
        //vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);
        vcg::tri::UnMarkAll(mesh);
        for (size_t i=0;i<PatchFaces.size();i++)
            for (size_t j=0;j<3;j++)
            {
                size_t IndexV=vcg::tri::Index(mesh,mesh.face[PatchFaces[i]].V(j));
                //if (mesh.vert[IndexV].IsV())continue;
                if (vcg::tri::IsMarked(mesh,&mesh.vert[IndexV]))continue;
                vcg::tri::Mark(mesh,&mesh.vert[IndexV]);
                //mesh.vert[IndexV].SetV();
                ret+=VerticesNeeds[IndexV];
            }
        return ret;
    }


    static void GetBorderSequenceFrom(MeshType &mesh,
                                      const vcg::face::Pos<typename MeshType::FaceType> &StartPos,
                                      std::vector<vcg::face::Pos<typename MeshType::FaceType> > &SeqPos,
                                      std::vector<size_t> &BorderSequences)
    {
        typedef typename MeshType::FaceType FaceType;
        typedef typename vcg::face::Pos<FaceType> PosType;

        BorderSequences.clear();
        PosType CurrPos=StartPos;
        SeqPos.clear();
        //CurrPos.F()->SetS();
        //search if there is a selected vertex
        //std::cout<<"a"<<std::endl;
        do
        {
            if (CurrPos.VFlip()->IsS())
                CurrPos.FlipV();

            if (!CurrPos.V()->IsS())
                CurrPos.NextB();

        }while((!CurrPos.V()->IsS())&&(CurrPos!=StartPos));
        //std::cout<<"b"<<std::endl;
        //CurrPos.F()->SetS();

        //then start
        size_t IndexV=vcg::tri::Index(mesh,CurrPos.V());
        BorderSequences.push_back(IndexV);

        //go on the other side
        CurrPos.FlipV();
        assert(CurrPos.IsBorder());
        PosType Pos0=CurrPos;
        SeqPos.push_back(CurrPos);
        do
        {
            //CurrPos.F()->SetS();
            size_t IndexV=vcg::tri::Index(mesh,CurrPos.V());
            BorderSequences.push_back(IndexV);
            if (!CurrPos.V()->IsS())
            {
                CurrPos.NextB();
                SeqPos.push_back(CurrPos);
            }
        }while((!CurrPos.V()->IsS())&&(CurrPos!=Pos0));
        //    std::cout<<"c"<<std::endl;
        if (CurrPos.V()->IsS())
        {
            size_t IndexV=vcg::tri::Index(mesh,CurrPos.V());
            BorderSequences.push_back(IndexV);
        }
    }

    static void GetSequencesLenghtAngle(MeshType &mesh,
                                        std::vector<size_t> &BorderSequence,
                                        typename MeshType::ScalarType &angle,
                                        typename MeshType::ScalarType &lenght)
    {
        typedef typename MeshType::CoordType CoordType;
        typedef typename MeshType::ScalarType ScalarType;

        angle=0;lenght=0;
        assert(BorderSequence.size()>1);
        CoordType Dir0=mesh.vert[BorderSequence[1]].P()-mesh.vert[BorderSequence[0]].P();
        Dir0.Normalize();
        for (size_t i=0;i<BorderSequence.size()-1;i++)
        {
            int IndxV0=BorderSequence[i];
            int IndxV1=BorderSequence[(i+1)];
            CoordType Dir1=(mesh.vert[IndxV1].P()-mesh.vert[IndxV0].P());
            lenght+=Dir0.Norm();
            Dir0.Normalize();
            Dir1.Normalize();
            ScalarType angleInt=fabs(vcg::Angle(Dir0,Dir1));
            angle+=angleInt;
            Dir0=Dir1;
        }
    }

    //get individual sequences of sharp features

    static void GetBorderSequences(MeshType &mesh,const std::vector<size_t> &Corners,
                                   std::vector<std::vector<size_t> > &BorderSequences,
                                   bool DebugMsg=false)
    {

        if (DebugMsg)
            std::cout<<"Deriving Sequences"<<std::endl;
        vcg::tri::UpdateSelection<MeshType>::Clear(mesh);
        SelectVertices(mesh,Corners);

        vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(mesh);
        bool found_start=false;

        for (size_t i=0;i<mesh.face.size();i++)
        {
            //find first border edge which has not been already processed
            for (size_t j=0;j<3;j++)
            {
                if (!vcg::face::IsBorder(mesh.face[i],j))continue;
                if (mesh.face[i].IsFaceEdgeS(j))continue;
                found_start=true;
                PosType StartPos=PosType(&mesh.face[i],j);

                std::vector<size_t> CurrBorderSequence;
                //std::cout<<"a"<<std::endl;
                assert(StartPos.IsBorder());
                std::vector<PosType> PosBorderSeq;
                GetBorderSequenceFrom(mesh,StartPos,PosBorderSeq,CurrBorderSequence);

                for (size_t j=0;j<PosBorderSeq.size();j++)
                    PosBorderSeq[j].F()->SetFaceEdgeS(PosBorderSeq[j].E());

                BorderSequences.push_back(CurrBorderSequence);

            }
        }

        if (DebugMsg)
            std::cout<<"Making Globally Consistent "<<std::endl;
        //then sort them globally
        for (size_t i=0;i<BorderSequences.size();i++)
        {
            bool IsLoop=(BorderSequences[i][0]==BorderSequences[i].back());
            if (IsLoop)
            {
                //erase last repeated element
                BorderSequences[i].pop_back();
                //get the smallest
                size_t smallestI=0;
                CoordType smallestPos=mesh.vert[BorderSequences[i][0]].P();
                for (size_t j=1;j<BorderSequences[i].size();j++)
                {
                    CoordType currPos=mesh.vert[BorderSequences[i][j]].P();
                    if (currPos>smallestPos)continue;
                    smallestPos=currPos;
                    smallestI=j;
                }

                //then re-order
                std::vector<size_t> NewBorderSeq;
                size_t sizeSeq=BorderSequences[i].size();
                for (size_t j=0;j<sizeSeq;j++)
                {
                    size_t IndexV=(j+smallestI)%sizeSeq;
                    NewBorderSeq.push_back(BorderSequences[i][IndexV]);
                }
                BorderSequences[i]=NewBorderSeq;
                //then reverse if needed
                CoordType Pos0=mesh.vert[BorderSequences[i][1]].P();
                CoordType Pos1=mesh.vert[BorderSequences[i].back()].P();
                if (Pos1<Pos0)
                    std::reverse(BorderSequences[i].begin()+1,BorderSequences[i].end());
            }else
            {
                //then reverse if needed
                CoordType Pos0=mesh.vert[BorderSequences[i][0]].P();
                CoordType Pos1=mesh.vert[BorderSequences[i].back()].P();
                if (Pos1<Pos0)
                    std::reverse(BorderSequences[i].begin(),BorderSequences[i].end());
            }
        }
        if (DebugMsg)
            std::cout<<"Done with Initializing Border Sequences"<<std::endl;
    }



    static void GetPatchInfo(MeshType &mesh,std::vector<std::vector<size_t> > &PatchFaces,
                             std::vector<std::vector<size_t> > &PatchCorners,std::vector<size_t> &VerticesNeeds,
                             //std::unordered_map<std::pair<size_t,size_t>,typename MeshType::ScalarType>  &EdgeMap,
                             std::map<std::pair<size_t,size_t>,typename MeshType::ScalarType>  &EdgeMap,
                             std::vector<PatchInfo<typename MeshType::ScalarType> > &PInfo,
                             const typename MeshType::ScalarType Thr,bool SkipValence4,
                             bool AllowDarts,bool AllowSelfGluedPatches,
                             bool CheckQuadrangulationLimits)
    {
        //typedef typename MeshType::ScalarType ScalarType;
        PInfo.clear();
        PInfo.resize(PatchFaces.size());

        SelectMeshPatchBorders(mesh,PatchFaces,true);
        size_t int0=0;
        size_t int1=0;
        size_t int2=0;
        size_t int3=0;
        size_t int4=0;
        for (size_t i=0;i<PatchFaces.size();i++)
        {
            size_t t0=clock();
            PInfo[i].NumCorners=PatchCorners[i].size();
            PInfo[i].NumEmitters=RemainingEmitters(mesh,PatchFaces[i],VerticesNeeds);
            size_t t1=clock();
            int0+=t1-t0;

            //if ((!AllowDarts)&&(!AllowSelfGluedPatches))
            //initialize genus by default
            PInfo[i].Genus=PatchGenus(mesh,PatchFaces[i],AllowDarts,AllowSelfGluedPatches);

//            //check if there is internal path not connected
//            if ((AllowDarts)&&(!AllowSelfGluedPatches))
//            {
//                //in such case check if there is no internal holes
//                //(but only if the mesh is not glued to itself)
//                if (PInfo[i].Genus==1)
//                {
//                    bool SingleHoles=SingleConnectedBorderFromEdgeS(mesh,PatchFaces[i]);
//                    if (SingleHoles)
//                        PInfo[i].Genus=1;
//                    else
//                        PInfo[i].Genus=2;
//                }
//            }
//            if ((!AllowDarts)&&(AllowSelfGluedPatches))
//            {
//                //only if self glue is allowed
//                bool SingleHoles=SingleConnectedBorderFromEdgeS(mesh,PatchFaces[i]);
//                if (SingleHoles)
//                    PInfo[i].Genus=1;
//                else
//                    PInfo[i].Genus=2;
//                //PInfo[i].Genus=PatchGenusCopyMesh(mesh,PatchFaces[i]);
//            }

//            if ((AllowDarts)&&(AllowSelfGluedPatches))
//            {
//                //only if self glue is allowed
//                bool SingleHoles=SingleConnectedBorderFromEdgeS(mesh,PatchFaces[i]);
//                if (SingleHoles)
//                    PInfo[i].Genus=1;
//                else
//                    PInfo[i].Genus=2;
//                //PInfo[i].Genus=PatchGenusCopyMesh(mesh,PatchFaces[i]);
//            }


            size_t t2=clock();
            int1=t2-t1;

            bool SingOnCorner;
            int ExpVal=ExpectedValence(mesh,PatchFaces[i],PatchCorners[i],SingOnCorner);

            //tolerant wrt corner singularities
            if (!SingOnCorner)
                PInfo[i].ExpectedValence=ExpVal;
            else
                PInfo[i].ExpectedValence=PInfo[i].NumCorners;

            size_t t3=clock();

            int2+=t3-t2;

            PInfo[i].NumSing=NumSingularities(mesh,PatchFaces[i],PatchCorners[i]);
            size_t t4=clock();
            int3+=t4-t3;


            //            if ((PInfo[i].NumCorners<(int)MIN_ADMITTIBLE)||
            //                    (PInfo[i].NumCorners>(int)MAX_ADMITTIBLE)||
            //                    (PInfo[i].Genus!=1)||
            //                    (PInfo[i].NumEmitters>0)||
            //                    (Thr<=0))
            if ((CheckQuadrangulationLimits && (PInfo[i].NumCorners<(int)MIN_ADMITTIBLE))||
                    (CheckQuadrangulationLimits && (PInfo[i].NumCorners>(int)MAX_ADMITTIBLE))||
                    (PInfo[i].Genus!=1)||
                    (PInfo[i].NumEmitters>0)||
                    (Thr<=0))
            {
                PInfo[i].CClarkability=false;
            }
            else
            {
                std::vector<typename MeshType::ScalarType > CurvedL,EuclL;
                PatchSideLenght(mesh,PatchFaces[i],PatchCorners[i],CurvedL,EuclL,EdgeMap);
                PInfo[i].CurvedL=CurvedL;
                //PInfo[i].CornerL=EuclL;
                //PInfo[i].CClarkability=IsCatmullClarkable(PatchFaces[i].size(),CurvedL,Thr,SkipValence4);
                PInfo[i].CClarkability=IsCatmullClarkable(CurvedL,Thr,SkipValence4);
                //            PInfo[i].CurvedL=CurvedL[i];
                //            PInfo[i].CornerL=EuclL[i];
                //            PInfo[i].CClarkability=IsCatmullClarkable(PatchFaces[i].size(),PInfo[i].CurvedL,Thr);
            }


            size_t t5=clock();
            int4+=t5-t4;
        }

        //        std::cout<<"** Timing Patch Type Update **"<<std::endl;
        //        std::cout<<"Time Step Emitters "<<int0<<std::endl;
        //        std::cout<<"Time Step Genus "<<int1<<std::endl;
        //        std::cout<<"Time Step Exp Val "<<int2<<std::endl;
        //        std::cout<<"Time Step Num Sing "<<int3<<std::endl;
        //        std::cout<<"Time Step CCability "<<int4<<std::endl<<std::endl;
    }



    static ScalarType PatchArea(MeshType &mesh,const std::vector<size_t> &PatchFaces)
    {
        MeshType patch_mesh;
        GetMeshFromPatch(mesh,PatchFaces,patch_mesh);
        return(MeshArea(patch_mesh));
    }


    typename MeshType::ScalarType NonOkArea(MeshType &mesh,
                                            const std::vector<std::vector<size_t> > &PatchFaces,
                                            const std::vector<PatchInfo<typename MeshType::ScalarType> > &PInf,
                                            size_t MinSides,
                                            size_t MaxSides,
                                            std::vector<typename MeshType::ScalarType> &QThresold)
    {
        typename MeshType::ScalarType CurrA=0;
        assert(PInf.size()==PatchFaces.size());
        for (size_t i=0;i<PatchFaces.size();i++)
        {
            if (IsOKPatch(PInf[i],MinSides,MaxSides,QThresold))continue;
            CurrA+=PatchArea(mesh,PatchFaces[i]);
        }
        return CurrA;
    }

    //    static bool BetterConfiguration(//MeshType &mesh,
    //                                    //                                    const std::vector<std::vector<size_t> > &PatchFaces0,
    //                                    //                                    const std::vector<std::vector<size_t> > &PatchFaces1,
    //                                    const std::vector<PatchInfo<typename MeshType::ScalarType> > &PInf0,
    //                                    const std::vector<PatchInfo<typename MeshType::ScalarType> > &PInf1,
    //                                    size_t MinSides,size_t MaxSides,
    //                                    const typename MeshType::ScalarType CClarkability,
    //                                    const typename MeshType::ScalarType avgEdge,
    //                                    bool match_sing_valence,
    //                                    bool print_debug=false)
    static bool BetterConfiguration(const std::vector<PatchInfo<typename MeshType::ScalarType> > &PInf0,
                                    const std::vector<PatchInfo<typename MeshType::ScalarType> > &PInf1,
                                    size_t MinSides,size_t MaxSides,
                                    const typename MeshType::ScalarType CClarkability,
                                    const typename MeshType::ScalarType avgEdge,
                                    bool match_sing_valence,
                                    bool CountEmitters,
                                    bool print_debug)
    {
        //        std::cout<<"CClarkability"<<CClarkability<<std::endl;
        //        std::cout<<"AvgEdge"<<avgEdge<<std::endl;
        //        std::cout<<"Match Sing"<<match_sing_valence<<std::endl;
        size_t NonOKGenus0=0;
        size_t NonOKGenus1=0;
        size_t NonOKEmitters0=0;
        size_t NonOKEmitters1=0;
        size_t NonOKSize0=0;
        size_t NonOKSize1=0;
        size_t Sing0=0;
        size_t Sing1=0;
        ScalarType NonQArea0=0;
        ScalarType NonQArea1=0;

        //at least one sing inside
        int MatchSing0=0;
        int MatchSing1=0;
        int MaxInternalSing0=1;
        int MaxInternalSing1=1;

        if (print_debug)
            std::cout<<"*** REMOVAL STATS ***"<<std::endl;

        if (print_debug)
            std::cout<<"Num Patches 0:"<<PInf0.size()<<std::endl;
        for (size_t i=0;i<PInf0.size();i++)
        {
            if (print_debug)
                std::cout<<"Num Sides 0:"<<PInf0[i].NumCorners<<std::endl;

            if (PInf0[i].Genus!=1)NonOKGenus0++;
            if ((CountEmitters)&&(PInf0[i].NumEmitters>0))NonOKEmitters0++;
            //(!allow_darts)&&
            if (PInf0[i].NumCorners<(int)MinSides)NonOKSize0++;
            if (PInf0[i].NumCorners>(int)MaxSides)NonOKSize0++;
            NonQArea0+=PInf0[i].QualityFunctorValue;
            if (match_sing_valence)
            {
                MaxInternalSing0=std::max(MaxInternalSing0,PInf0[i].NumSing);

                if ((PInf0[i].ExpectedValence!=4)&&
                        (PInf0[i].NumCorners==(int)PInf0[i].ExpectedValence))
                    MatchSing0++;
            }
            //if match singularity than no need to check CCability for valence 4
            if (CClarkability>0)
            {
                //Sing0+=AddedSingularities(PatchFaces0[i].size(),PInf0[i].CurvedL,
                //                          CClarkability*avgEdge,match_sing_valence,print_debug);
                Sing0+=AddedSingularities(PInf0[i].CurvedL,CClarkability*avgEdge,match_sing_valence,print_debug);
            }
        }

        if (print_debug)
            std::cout<<"Num Patches 1:"<<PInf1.size()<<std::endl;
        for (size_t i=0;i<PInf1.size();i++)
        {
            if (print_debug)
                std::cout<<"Num Sides 1:"<<PInf1[i].NumCorners<<std::endl;
            if (PInf1[i].Genus!=1)NonOKGenus1++;
            if ((CountEmitters)&&(PInf1[i].NumEmitters>0))NonOKEmitters1++;
            //(!allow_darts)&&
            if (PInf1[i].NumCorners<(int)MinSides)NonOKSize1++;
            if (PInf1[i].NumCorners>(int)MaxSides)NonOKSize1++;
            NonQArea1+=PInf1[i].QualityFunctorValue;
            if (match_sing_valence)
            {
                MaxInternalSing1=std::max(MaxInternalSing1,PInf1[i].NumSing);

                if ((PInf1[i].ExpectedValence!=4)&&
                        (PInf1[i].NumCorners==(int)PInf1[i].ExpectedValence))
                    MatchSing1++;
            }
            //if match singularity than no need to check CCability for valence 4
            if (CClarkability>0)
            {
                //                Sing1+=AddedSingularities(PatchFaces1[i].size(),PInf1[i].CurvedL,
                //                                          CClarkability*avgEdge,match_sing_valence,
                //                                          print_debug);
                Sing1+=AddedSingularities(PInf1[i].CurvedL,CClarkability*avgEdge,
                                          match_sing_valence,print_debug);
            }

        }

        if (print_debug)
        {
            std::cout<<"Non Ok Genus 0:"<<NonOKGenus0<<std::endl;
            std::cout<<"Non Ok Genus 1:"<<NonOKGenus1<<std::endl;
            std::cout<<"Non Ok Emitters 0:"<<NonOKEmitters0<<std::endl;
            std::cout<<"Non Ok Emitters 1:"<<NonOKEmitters1<<std::endl;
            std::cout<<"Non Ok NonOKSize 0:"<<NonOKSize0<<std::endl;
            std::cout<<"Non Ok NonOKSize 1:"<<NonOKSize1<<std::endl;
            std::cout<<"MaxInternalSing 0:"<<MaxInternalSing0<<std::endl;
            std::cout<<"MaxInternalSing 1:"<<MaxInternalSing1<<std::endl;
            std::cout<<"MatchSing 0:"<<MatchSing0<<std::endl;
            std::cout<<"MatchSing 1:"<<MatchSing1<<std::endl;
            std::cout<<"Sing 0:"<<Sing0<<std::endl;
            std::cout<<"Sing 1:"<<Sing1<<std::endl;
        }
        //    std::cout<<"CClarkability:"<<CClarkability<<std::endl;
        //    std::cout<<"avgEdge:"<<avgEdge<<std::endl;
        if (print_debug)
            std::cout<<"*** END REMOVAL STATS ***"<<std::endl;

        if (NonOKGenus1!=NonOKGenus0)
            return (NonOKGenus1<NonOKGenus0);

//        if (NonQArea0!=NonQArea1)
//            return (NonQArea1<NonQArea0);

        if (NonOKEmitters1!=NonOKEmitters0)
            return (NonOKEmitters1<NonOKEmitters0);

        //    if (NonOKSize0>0)
        //        return false;

        //    if (NonOKSize1>0)
        //        return false;
        if (NonOKSize1!=NonOKSize0)
            return (NonOKSize1<NonOKSize0);

        if (MaxInternalSing1!=MaxInternalSing0)
            return (MaxInternalSing1<MaxInternalSing0);

        if (MatchSing1!=MatchSing0)
            return(MatchSing1>MatchSing0);

        //if ((Sing0==0)&&(Sing1==0))return true;

        if (Sing1!=Sing0)
            return (Sing1<Sing0);


//        if (NonQArea0!=NonQArea1)
//            return (NonQArea1<NonQArea0);
        if (fabs(NonQArea0-NonQArea1)>(NonQArea0*0.000001))
            return (NonQArea1<NonQArea0);
        return true;
        //return (Sing1<Sing0);
    }


    static void MeshTrace(const VertexFieldGraph<MeshType> &VFGraph,
                          const CandidateTrace &CurrTrace,
                          typename MeshType::ScalarType Width,
                          MeshType &outMesh)
    {
        typedef typename MeshType::CoordType CoordType;
        //typedef typename MeshType::ScalarType ScalarType;
        outMesh.Clear();
        size_t Limit=CurrTrace.PathNodes.size();
        assert(Limit>=2);
        //std::cout<<"Limit "<<Limit<<std::endl;

        size_t Size=CurrTrace.PathNodes.size();
        if (CurrTrace.IsLoop)Limit++;
        for (size_t i=0;i<Limit-1;i++)
        {
            size_t N0=CurrTrace.PathNodes[i];
            size_t N1=CurrTrace.PathNodes[(i+1)%Size];
            //       std::cout<<"N0 "<<N0<<std::endl;
            //       std::cout<<"N1 "<<N1<<std::endl;
            CoordType P0=VFGraph.NodePos(N0);
            CoordType P1=VFGraph.NodePos(N1);
            //      std::cout<<"got position"<<std::endl;
            MeshType TempMesh;
            vcg::tri::OrientedCylinder(TempMesh,P0,P1,Width,true,8,1);
            vcg::tri::Append<MeshType,MeshType>::Mesh(outMesh,TempMesh);
        }
    }

    static void MeshTraces(const VertexFieldGraph<MeshType> &VFGraph,
                           const std::vector<CandidateTrace> &TraceSet,
                           const std::vector<bool> &Selected,
                           MeshType &outMesh)
    {
        assert(Selected.size()==TraceSet.size());
        typedef typename MeshType::ScalarType ScalarType;
        ScalarType Width=VFGraph.Mesh().bbox.Diag()*0.002;
        outMesh.Clear();
        for (size_t i=0;i<TraceSet.size();i++)
        {
            MeshType traceMesh;
            MeshTrace(VFGraph,TraceSet[i],Width,traceMesh);
            vcg::Color4b currCol=vcg::Color4b::Scatter(TraceSet.size(),i);
            if (Selected[i])
                currCol=vcg::Color4b::Red;

            vcg::tri::UpdateColor<MeshType>::PerFaceConstant(traceMesh,currCol);
            vcg::tri::Append<MeshType,MeshType>::Mesh(outMesh,traceMesh);
        }
    }

    //NOTE if Dest==0 then the function check if is a single connected component
    static bool IsConnectedTo(size_t SourceN, std::vector<size_t> &Dest,
                              std::map<size_t,std::vector<size_t> > &Adj)
    {
        //if (Dest.size()==0)return false;
        //if ((Dest.size()==1)&&(Dest[0]==SourceN))return false;
        if (Dest.size()==1)
        {
            if (Dest[0]==SourceN)return false;
        }

        std::vector<size_t> to_explore;
        std::set<size_t> has_explored;
        std::set<size_t> DestSet(Dest.begin(),Dest.end());

        to_explore.push_back(SourceN);
        has_explored.insert(SourceN);

        do{
            assert(to_explore.size()>0);
            size_t curr_node=to_explore.back();
            assert(has_explored.count(curr_node)>0);
            to_explore.pop_back();
            for (size_t i=0;i<Adj[curr_node].size();i++)
            {
                size_t adj_node=Adj[curr_node][i];
                if (has_explored.count(adj_node)>0)continue;

                to_explore.push_back(adj_node);
                has_explored.insert(adj_node);

                if (DestSet.count(adj_node)>0)return true;
            }
        }while (to_explore.size()>0);

        if (Dest.size()==0)
            return (Adj.size()==has_explored.size());

        return false;
    }

    static void MarkFaces(MeshType &mesh,std::vector<size_t> &PatchFaces)
    {
        //mark the faces
        vcg::tri::UnMarkAll(mesh);
        for (size_t i=0;i<PatchFaces.size();i++)
        {
            assert(PatchFaces[i]>=0);
            assert(PatchFaces[i]<mesh.face.size());
            vcg::tri::Mark(mesh,&mesh.face[PatchFaces[i]]);
        }
    }

    static bool HasInternalFeaturesFromEdgeSel(MeshType &mesh,
                                               std::vector<size_t> &PatchFaces)
    {
        MarkFaces(mesh,PatchFaces);

        for (size_t i=0;i<PatchFaces.size();i++)
        {
            size_t IndexF=PatchFaces[i];
            assert(IndexF>=0);
            assert(IndexF<mesh.face.size());
            assert(vcg::tri::IsMarked(mesh,&mesh.face[IndexF]));

            for (size_t j=0;j<3;j++)
            {
                if (!mesh.face[IndexF].IsFaceEdgeS(j))continue;

                bool IsOnBorder=vcg::face::IsBorder(mesh.face[IndexF],j);
                FaceType *fOpp=mesh.face[IndexF].FFp(j);
                assert(fOpp!=NULL);
                IsOnBorder|=(!vcg::tri::IsMarked(mesh,fOpp));

                if (!IsOnBorder)
                    return true;
            }
        }
        return false;
    }

//    static void GetPatchBorderV(MeshType &mesh,
//                                std::vector<size_t> &PatchFaces,
//                                std::vector<size_t> &BorderV)
//    {
//        MarkFaces(mesh,PatchFaces);
//        BorderV.clear();
//        for (size_t i=0;i<PatchFaces.size();i++)
//        {
//            size_t IndexF=PatchFaces[i];
//            assert(IndexF<mesh.face.size());

//            assert(vcg::tri::IsMarked(mesh,&mesh.face[IndexF]));

//            for (size_t j=0;j<3;j++)
//            {
//                bool IsOnBorder=vcg::face::IsBorder(mesh.face[IndexF],j);
//                FaceType *fOpp=mesh.face[IndexF].FFp(j);
//                IsOnBorder|=(!vcg::tri::IsMarked(mesh,fOpp));

//                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[IndexF].V0(j));
//                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[IndexF].V1(j));

//                if (IsOnBorder)
//                {
//                    BorderV.push_back(IndexV0);
//                    BorderV.push_back(IndexV1);
//                }
//            }
//        }
//        std::sort(BorderV.begin(),BorderV.end());
//        std::vector<size_t>::iterator it;
//        it = std::unique (BorderV.begin(),BorderV.end());
//        BorderV.resize( std::distance(BorderV.begin(),it) );
//    }

    static void GetPatchBorderV(MeshType &mesh,
                                std::vector<size_t> &PatchFaces,
                                std::vector<size_t> &BorderV)
    {
        BorderV.clear();
        std::map<std::pair<CoordType,CoordType>,size_t> EdgeNum;

        for (size_t i=0;i<PatchFaces.size();i++)
        {
            size_t IndexF=PatchFaces[i];
            assert(IndexF<mesh.face.size());
            for (size_t j=0;j<3;j++)
            {
                CoordType Pe0=mesh.face[i].P0(j);
                CoordType Pe1=mesh.face[i].P1(j);
                std::pair<CoordType,CoordType> E(std::min(Pe0,Pe1),std::max(Pe0,Pe1));
                if (EdgeNum.count(E)==0)
                    EdgeNum[E]=1;
                else
                    EdgeNum[E]++;
            }
        }
        for (size_t i=0;i<PatchFaces.size();i++)
        {
            size_t IndexF=PatchFaces[i];
            assert(IndexF<mesh.face.size());
            for (size_t j=0;j<3;j++)
            {
                CoordType Pe0=mesh.face[i].P0(j);
                CoordType Pe1=mesh.face[i].P1(j);
                std::pair<CoordType,CoordType> E(std::min(Pe0,Pe1),std::max(Pe0,Pe1));
                int NumE=EdgeNum[E];
                assert((NumE==1)||(NumE==2));

                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[IndexF].V0(j));
                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[IndexF].V1(j));

                bool IsOnBorder=(NumE==1);
                if (IsOnBorder)
                {
                    BorderV.push_back(IndexV0);
                    BorderV.push_back(IndexV1);
                }
            }
        }
        std::sort(BorderV.begin(),BorderV.end());
        std::vector<size_t>::iterator it;
        it = std::unique (BorderV.begin(),BorderV.end());
        BorderV.resize( std::distance(BorderV.begin(),it) );
    }

    static void GetSelectedEdgesGraph(MeshType &mesh,
                                      std::vector<size_t> &PatchFaces,
                                      std::map<size_t,std::vector<size_t> > &Adj,
                                      bool OnlyInternal)
    {
        MarkFaces(mesh,PatchFaces);
        Adj.clear();

        for (size_t i=0;i<PatchFaces.size();i++)
        {
            size_t IndexF=PatchFaces[i];

            assert(IndexF<mesh.face.size());
            assert(vcg::tri::IsMarked(mesh,&mesh.face[IndexF]));

            for (size_t j=0;j<3;j++)
            {
                bool IsOnBorder=vcg::face::IsBorder(mesh.face[IndexF],j);
                FaceType *fOpp=mesh.face[IndexF].FFp(j);
                assert(fOpp!=NULL);

                IsOnBorder|=(!vcg::tri::IsMarked(mesh,fOpp));

                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[IndexF].V0(j));
                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[IndexF].V1(j));

                bool AddEdge=false;
                bool IsEdgeS=mesh.face[IndexF].IsFaceEdgeS(j);

                if (IsOnBorder && (!OnlyInternal))
                {
                    //assert(IsEdgeS);
                    AddEdge=true;
                }

                if ((!IsOnBorder) && (IsEdgeS))
                    AddEdge=true;

                if (AddEdge)
                {
                    Adj[IndexV0].push_back(IndexV1);
                    Adj[IndexV1].push_back(IndexV0);
                }
            }
        }
    }


    static void GetSelectedEdgesGraph(MeshType &mesh,
                                      std::map<size_t,std::vector<size_t> > &Adj,
                                      bool OnlyInternal)
    {
        Adj.clear();
        std::set<std::pair<size_t,size_t> > VisitedE;
        for (size_t i=0;i<mesh.face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                bool IsOnBorder=vcg::face::IsBorder(mesh.face[i],j);

                if (OnlyInternal && IsOnBorder)continue;

                bool IsEdgeS=mesh.face[i].IsFaceEdgeS(j);

                if (!IsEdgeS)continue;

                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));

                std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),
                                             std::max(IndexV0,IndexV1));

                if (VisitedE.count(key)>0)continue;

                    Adj[IndexV0].push_back(IndexV1);
                    Adj[IndexV1].push_back(IndexV0);
            }
        }
    }

//    static void GetCornerValuesFromInternalFeatures(MeshType &mesh,
//                                                    std::vector<size_t> &PatchFaces,
//                                                    std::vector<size_t> &Corners,
//                                                    std::vector<size_t> &CornerValence)
//    {
//        //mark the faces
//        vcg::tri::UnMarkAll(mesh);
//        for (size_t i=0;i<PatchFaces.size();i++)
//            vcg::tri::Mark(mesh,&mesh.face[PatchFaces[i]]);

//        //create the adjacency graph
//        std::map<size_t,std::vector<size_t> > Adj;
//        std::vector<size_t> BorderV;
//        for (size_t i=0;i<PatchFaces.size();i++)
//        {
//            size_t IndexF=PatchFaces[i];
//            assert(vcg::tri::IsMarked(mesh,&mesh.face[IndexF]));

//            for (size_t j=0;j<3;j++)
//            {
//                bool IsOnBorder=vcg::face::IsBorder(mesh.face[IndexF],j);
//                FaceType *fOpp=mesh.face[IndexF].FFp(j);
//                IsOnBorder|=(!vcg::tri::IsMarked(mesh,fOpp));

//                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[IndexF].V0(j));
//                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[IndexF].V1(j));

//                if (IsOnBorder)
//                {
//                    BorderV.push_back(IndexV0);
//                    BorderV.push_back(IndexV1);
//                }
//                else
//                {
//                    if (mesh.face[IndexF].IsFaceEdgeS(j))
//                    {
//                        Adj[IndexV0].push_back(IndexV1);
//                        Adj[IndexV1].push_back(IndexV0);
//                    }
//                }
//            }
//        }

//        std::sort(BorderV.begin(),BorderV.end());
//        std::vector<size_t>::iterator it;
//        it = std::unique (BorderV.begin(),BorderV.end());
//        BorderV.resize( std::distance(BorderV.begin(),it) );

//        //collect sources
//        std::vector<size_t> Sources;
//        for (size_t i=0;i<BorderV.size();i++)
//        {
//            size_t IndexV=BorderV[i];
//            if (Adj.count(IndexV)>0)
//                Sources.push_back(IndexV);
//        }
//        //std::cout<<"There are:"<<Sources.size()<<" sources"<<std::endl;
//        CornerValence=std::vector<size_t>(Corners.size(),1);

//        //then check for each source if connected to another source
//        for (size_t i=0;i<Sources.size();i++)
//        {
//            size_t IndexV=Sources[i];
//            bool HasConn=IsConnectedTo(IndexV,Sources,Adj);

//            std::vector<size_t>::iterator IteC;
//            IteC=std::find(Corners.begin(),Corners.end(),IndexV);
//            assert(IteC!=Corners.end());
//            size_t Index=distance(Corners.begin(),IteC);

//            if (HasConn)
//                CornerValence[Index]=2;
//            else
//                CornerValence[Index]=0;
//        }
//    }

//    static void GetCornerValuesFromInternalFeatures(MeshType &mesh,
//                                                    std::vector<size_t> &PatchFaces,
//                                                    std::vector<size_t> &Corners,
//                                                    std::vector<size_t> &CornerValence)
//    {
//        CornerValence=std::vector<size_t>(Corners.size(),1);
//        if (!HasInternalFeaturesFromEdgeSel(mesh,PatchFaces))
//            return;

//        std::vector<size_t> BorderV;
//        GetPatchBorderV(mesh,PatchFaces,BorderV);


//        std::map<size_t,std::vector<size_t> > Adj;
//        GetSelectedEdgesGraph(mesh,PatchFaces,Adj,true);

//        //collect sources that are also on internal features
//        std::vector<size_t> Sources;
//        for (size_t i=0;i<BorderV.size();i++)
//        {
//            size_t IndexV=BorderV[i];
//            if (Adj.count(IndexV)>0)
//                Sources.push_back(IndexV);
//        }


//        //then check for each source if connected to another source only with internal path
//        // in such case the vertex is duplicated as is a patch glued to itself
//        for (size_t i=0;i<Sources.size();i++)
//        {
//            size_t IndexV=Sources[i];
//            bool HasConn=IsConnectedTo(IndexV,Sources,Adj);

//            std::vector<size_t>::iterator IteC;
//            IteC=std::find(Corners.begin(),Corners.end(),IndexV);
//            if(IteC==Corners.end())continue;//it might be a concave and so not marked as a corner from that side
////            if (IteC==Corners.end())
////            {
////                std::cout<<"There are:"<<Sources.size()<<" sources"<<std::endl;
////                std::cout<<"There are:"<<Corners.size()<<" corners"<<std::endl;

////                std::cout<<"Sources"<<std::endl;
////                for (size_t i=0;i<Sources.size();i++)
////                    std::cout<<Sources[i]<<std::endl;

////                std::cout<<"Corners"<<std::endl;
////                std::set<CoordType> CornersSet;
////                for (size_t i=0;i<Corners.size();i++)
////                {
////                    CornersSet.insert(mesh.vert[Corners[i]].P());
////                    //std::cout<<Corners[i]<<std::endl;
////                }

////                MeshType curr_patch;
////                GetMeshFromPatch(mesh,PatchFaces,curr_patch,true);

////                for (size_t i=0;i<curr_patch.vert.size();i++)
////                {
////                    if (CornersSet.count(curr_patch.vert[i].P())==0)continue;
////                    curr_patch.vert[i].SetS();
////                }
////                vcg::tri::io::ExporterPLY<MeshType>::Save(curr_patch,"test_patch.ply",
////                                                          vcg::tri::io::Mask::IOM_VERTFLAGS);
////                //assert(0);
////            }
//            size_t Index=distance(Corners.begin(),IteC);
//            assert(Index>=0);
//            assert(Index<CornerValence.size());
//            if (HasConn)
//                CornerValence[Index]=2;
//            else
//                CornerValence[Index]=0;
//        }
//    }

//    static void GetCornerValuesFromInternalFeatures(MeshType &mesh,
//                                                    std::vector<size_t> &PatchFaces,
//                                                    std::vector<size_t> &Corners,
//                                                    std::vector<size_t> &CornerValence)
//    {
//        CornerValence=std::vector<size_t>(Corners.size(),1);
//        if (!HasInternalFeaturesFromEdgeSel(mesh,PatchFaces))
//            return;

//        //GetMeshFromPatch(mesh,PatchFaces,PatchM,false);
//        std::vector<size_t> BorderV;
//        GetPatchBorderV(mesh,PatchFaces,BorderV);


//        std::map<size_t,std::vector<size_t> > Adj;
//        GetSelectedEdgesGraph(mesh,PatchFaces,Adj,true);

//        //collect sources that are also on internal features
//        std::vector<size_t> Sources;
//        for (size_t i=0;i<BorderV.size();i++)
//        {
//            size_t IndexV=BorderV[i];
//            if (Adj.count(IndexV)>0)
//                Sources.push_back(IndexV);
//        }
//        std::cout<<"Sources:"<<Sources.size()<<std::endl;

//        //then check for each source if connected to another source only with internal path
//        // in such case the vertex is duplicated as is a patch glued to itself
//        for (size_t i=0;i<Sources.size();i++)
//        {
//            size_t IndexV=Sources[i];
//            bool HasConn=IsConnectedTo(IndexV,Sources,Adj);

//            std::vector<size_t>::iterator IteC;
//            IteC=std::find(Corners.begin(),Corners.end(),IndexV);
//            if(IteC==Corners.end())continue;//it might be a concave and so not marked as a corner from that side
//            size_t Index=distance(Corners.begin(),IteC);
//            assert(Index>=0);
//            assert(Index<CornerValence.size());
//            if (HasConn)
//                CornerValence[Index]=2;
//            else
//                CornerValence[Index]=0;
//        }
//    }

    static void GetCornerValuesFromInternalFeatures(MeshType &mesh,
                                                    std::vector<size_t> &PatchFaces,
                                                    std::vector<size_t> &Corners,
                                                    std::vector<size_t> &CornerValence)
    {
        CornerValence=std::vector<size_t>(Corners.size(),1);
        if (!HasInternalFeaturesFromEdgeSel(mesh,PatchFaces))
            return;

        //get the mesh without internal cuts
        MeshType PatchM;
        GetMeshFromPatch(mesh,PatchFaces,PatchM,false);

        //this in case of features, so internal cuts
        vcg::tri::Clean<MeshType>::RemoveDuplicateVertex(PatchM);
        vcg::tri::Clean<MeshType>::RemoveUnreferencedVertex(PatchM);
        vcg::tri::Allocator<MeshType>::CompactEveryVector(PatchM);
        PatchM.UpdateAttributes();

//        std::vector<size_t> BorderV;
//        for (size_t i=0;i<PatchM.vert.s)
        //GetPatchBorderV(mesh,PatchFaces,BorderV);


        std::map<size_t,std::vector<size_t> > Adj;
        GetSelectedEdgesGraph(PatchM,Adj,true);

        //collect sources that are also on internal features
        std::vector<size_t> Sources;
        for (size_t i=0;i<PatchM.vert.size();i++)
        {
           if (!PatchM.vert[i].IsB())continue;
            if (Adj.count(i)>0)
                Sources.push_back(i);
        }
        //std::cout<<"Sources:"<<Sources.size()<<std::endl;

        std::map<CoordType,size_t> CornerIdx;
        for (size_t i=0;i<Corners.size();i++)
        {
            CoordType CurrP=mesh.vert[Corners[i]].P();
            CornerIdx[CurrP]=i;
        }

        //then check for each source if connected to another source only with internal path
        // in such case the vertex is duplicated as is a patch glued to itself
        for (size_t i=0;i<Sources.size();i++)
        {
            size_t IndexV=Sources[i];
            bool HasConn=IsConnectedTo(IndexV,Sources,Adj);
            CoordType CurrP=PatchM.vert[IndexV].P();

            if (CornerIdx.count(CurrP)==0)continue;

            size_t IndexC=CornerIdx[CurrP];

            assert(IndexC>=0);
            assert(IndexC<CornerValence.size());
            if (HasConn)
                CornerValence[IndexC]=2;
            else
                CornerValence[IndexC]=0;
//            std::vector<size_t>::iterator IteC;
//            IteC=std::find(Corners.begin(),Corners.end(),IndexV);
//            if(IteC==Corners.end())continue;//it might be a concave and so not marked as a corner from that side
//            size_t Index=distance(Corners.begin(),IteC);
//            assert(Index>=0);
//            assert(Index<CornerValence.size());
//            if (HasConn)
//                CornerValence[Index]=2;
//            else
//                CornerValence[Index]=0;
        }
    }

    static bool SingleConnectedBorderFromEdgeS(MeshType &mesh,std::vector<size_t> &PatchFaces)
    {

        std::vector<size_t> BorderV;
        GetPatchBorderV(mesh,PatchFaces,BorderV);

        //there is no border
        if (BorderV.size()==0)return false;

        //get a graph of all borders
        std::map<size_t,std::vector<size_t> > Adj;
        GetSelectedEdgesGraph(mesh,PatchFaces,Adj,false);

        //then check if single connected omponents
        std::vector<size_t> Dest;
        bool SingleConn=IsConnectedTo(BorderV[0],Dest,Adj);

        return SingleConn;
    }


   static void GetFacesSurroundingUsingEdgeSel(MeshType &mesh,
                                               std::vector<vcg::face::Pos<FaceType> > &StartPos,
                                               std::vector<size_t> &IndedF)
   {
       vcg::tri::UpdateFlags<MeshType>::FaceClearV(mesh);

       IndedF.clear();
       std::vector<size_t> StackF;
       for (size_t i=0;i<StartPos.size();i++)
       {
           FaceType *f0=StartPos[i].F();
           FaceType *f1=StartPos[i].FFlip();
           int IndexF0=vcg::tri::Index(mesh,f0);
           int IndexF1=vcg::tri::Index(mesh,f1);
           StackF.push_back(IndexF0);
           StackF.push_back(IndexF1);
           f0->SetV();
           f1->SetV();
       }

       IndedF=StackF;
       do{
           size_t currF=StackF.back();
           StackF.pop_back();
           FaceType *f=&mesh.face[currF];
           for (size_t i=0;i<3;i++)
           {
               if (vcg::face::IsBorder(*f,i))continue;
               if (f->IsFaceEdgeS(i))continue;//in this case cannot go over

               FaceType *fopp=f->FFp(i);
               if (fopp->IsV())continue;

               int IndexFopp=vcg::tri::Index(mesh,fopp);
               IndedF.push_back(IndexFopp);
               StackF.push_back(IndexFopp);
               fopp->SetV();
           }
       }while (!StackF.empty());

       std::sort(IndedF.begin(),IndedF.end());
       std::vector<size_t>::iterator it;
       it = std::unique (IndedF.begin(),IndedF.end());
       IndedF.resize( std::distance(IndedF.begin(),it) );
   }

   static void GetSubMeshSurroundingUsingEdgeSel(MeshType &mesh,
                                                 std::vector<vcg::face::Pos<FaceType> > &StartPos,
                                                 MeshType &subMesh)
   {
       subMesh.Clear();

       std::vector<size_t> IndedF;
       GetFacesSurroundingUsingEdgeSel(mesh,StartPos,IndedF);
       vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);
       vcg::tri::UpdateSelection<MeshType>::FaceClear(mesh);

       for (size_t i=0;i<IndedF.size();i++)
           mesh.face[IndedF[i]].SetS();

       vcg::tri::UpdateSelection<MeshType>::VertexFromFaceLoose(mesh);

       //then copy the submesh
       vcg::tri::Append<MeshType,MeshType>::Mesh(subMesh,mesh,true);

       subMesh.UpdateAttributes();

       //open if needed
       //vcg::tri::CutMeshAlongSelectedFaceEdges(subMesh);
       VertSplitter<MeshType>::SplitAlongEdgeSel(subMesh);

       vcg::tri::Clean<MeshType>::RemoveUnreferencedVertex(subMesh);
       vcg::tri::Allocator<MeshType>::CompactEveryVector(subMesh);
       subMesh.UpdateAttributes();
   }

   static bool IntersectUVBorderE(FaceType &f0,FaceType &f1)
   {
       for (size_t i=0;i<3;i++)
       {
           if (!vcg::face::IsBorder(f0,i))continue;

           vcg::Point2<ScalarType> P0UV=f0.V0(i)->T().P();
           vcg::Point2<ScalarType> P1UV=f0.V1(i)->T().P();
           vcg::Segment2<ScalarType> S0(P0UV,P1UV);

           for (size_t j=0;j<3;j++)
           {
               if (!vcg::face::IsBorder(f1,j))continue;

               vcg::Point2<ScalarType> P0UV=f1.V0(j)->T().P();
               vcg::Point2<ScalarType> P1UV=f1.V1(j)->T().P();
               vcg::Segment2<ScalarType> S1(P0UV,P1UV);

               if (S0.P0()==S1.P0())continue;
               if (S0.P0()==S1.P1())continue;
               if (S0.P1()==S1.P0())continue;
               if (S0.P1()==S1.P1())continue;

               vcg::Point2<ScalarType> UVInt;
               if (vcg::SegmentSegmentIntersection(S0,S1,UVInt))
                   return true;
           }
       }
       return false;
   }

   static bool SelfOverlapUV(MeshType &mesh)
   {
       for (size_t i=0;i<mesh.face.size()-1;i++)
       {
           FaceType *f0=&mesh.face[i];
           for (size_t j=(i+1);j<mesh.face.size();j++)
           {
              FaceType *f1=&mesh.face[j];
              if (IntersectUVBorderE(*f0,*f1))
                  return true;
           }
       }
       return false;
   }
};

#endif
