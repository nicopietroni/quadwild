#ifndef PATCH_MANAGER
#define PATCH_MANAGER

#include <wrap/igl/arap_parametrization.h>
#include <wrap/igl/lscm_parametrization.h>
#include <tracing/candidate_path.h>
#include "vert_field_graph.h"
#include <vcg/complex/algorithms/parametrization/distortion.h>
#include <vcg/space/outline2_packer.h>
#include "catmull_clarkability.h"

template <class MeshType>
typename MeshType::ScalarType MeshArea(MeshType &mesh)
{
    typedef typename MeshType::ScalarType ScalarType;

    ScalarType currA=0;
    for (size_t i=0;i<mesh.face.size();i++)
        currA+=(vcg::DoubleArea(mesh.face[i])/2.0);
    return currA;
}

template <class ScalarType>
struct PatchInfo
{
    int NumEmitters;
    int NumCorners;
    int Genus;
    std::vector<ScalarType> CornerL;
    std::vector<ScalarType> CurvedL;
    ScalarType CClarkability;
    int PossibleSing;
    std::vector<ScalarType> Q;

    PatchInfo()
    {
        NumEmitters=0;
        NumCorners=0;
        Genus=0;
        CClarkability=std::numeric_limits<ScalarType>::max();
        //PossibleSing=-1;
    }
};

template <class MeshType>
void SelectPatchFacesVert(MeshType &totMesh,const std::vector<size_t> &PatchFaces)
{
    vcg::tri::UpdateSelection<MeshType>::Clear(totMesh);
    for (size_t i=0;i<PatchFaces.size();i++)
        totMesh.face[PatchFaces[i]].SetS();

    vcg::tri::UpdateSelection<MeshType>::VertexFromFaceLoose(totMesh);
}

template <class MeshType>
int PatchGenus(MeshType &testMesh,const std::vector<size_t> &PatchFaces)
{
    SelectPatchFacesVert<MeshType>(testMesh,PatchFaces);

    std::set<std::pair<size_t,size_t> > EdgeSet;
    size_t NumF=0;
    size_t NumV=0;
    size_t NumE=0;
    for (size_t i=0;i<testMesh.face.size();i++)
    {
        if (!testMesh.face[i].IsS())continue;
        NumF++;
        for (size_t j=0;j<3;j++)
        {
            size_t IndV0=vcg::tri::Index(testMesh,testMesh.face[i].V0(j));
            size_t IndV1=vcg::tri::Index(testMesh,testMesh.face[i].V1(j));
            EdgeSet.insert(std::pair<size_t,size_t>(std::min(IndV0,IndV1),std::max(IndV0,IndV1)));
        }
    }
    for (size_t i=0;i<testMesh.vert.size();i++)
    {
        if (testMesh.vert[i].IsD())continue;
        if (testMesh.vert[i].IsS())NumV++;
    }

    NumE=EdgeSet.size();
    return ( NumV + NumF - NumE );
}

enum ParamType{Arap,LSQMap};

template <class MeshType>
void ComputeUV(MeshType &mesh, ParamType &PType,bool fixSVert)
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

template <class ScalarType>
void GetCornersUV(size_t numV,
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

template <class MeshType>
void ParametrizeCorners(MeshType &mesh,bool scaleEdges,
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

template <class MeshType>
void SelectCorners(MeshType &mesh,const std::vector<size_t> &CornersIDX)
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


template <class MeshType>
void GetIndexFromQ(MeshType &mesh,
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

template <class MeshType>
void DeriveBorderSeq(MeshType &mesh,std::vector<std::vector<size_t> > &BorderSeq)
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
    //assert(found);
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

template <class MeshType>
void DerivePatchBorderSeq(MeshType &mesh,
                          const std::vector<std::vector<size_t> > &PatchFaces,
                          const std::vector<std::vector<size_t> > &PatchCorners,
                          std::vector<std::vector<std::vector<size_t> > > &BorderSeq)
{
    BorderSeq.clear();
    for (size_t i=0;i<PatchFaces.size();i++)
    {
        BorderSeq.resize(BorderSeq.size()+1);
        MeshType subMesh;
        GetMeshFromPatch(mesh,PatchFaces[i],subMesh);
        std::vector<size_t> CornersIdx;
        GetIndexFromQ(subMesh,PatchCorners[i],CornersIdx);


        SelectCorners(subMesh,CornersIdx);
        DeriveBorderSeq(subMesh,BorderSeq[i]);

        //then set to original Index
        for (size_t j=0;j<BorderSeq[i].size();j++)
            for (size_t k=0;k<BorderSeq[i][j].size();k++)
                BorderSeq[i][j][k]=subMesh.vert[BorderSeq[i][j][k]].Q();
    }
}

template <class MeshType>
typename MeshType::ScalarType FieldLenght(const MeshType &mesh,
                                          const std::vector<size_t> &BorderSeq,
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

template <class MeshType>
typename MeshType::ScalarType Lenght(const MeshType &mesh,
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

template <class MeshType>
void GetLenghts(const MeshType &mesh,
                const std::vector<std::vector<size_t> > &BorderSeq,
                std::vector<typename MeshType::ScalarType> &BorderLen)
{
    for (size_t i=0;i<BorderSeq.size();i++)
        BorderLen.push_back(Lenght(mesh,BorderSeq[i]));
}

template <class MeshType>
void ParametrizeSeq(MeshType &mesh,const std::vector<size_t> &BorderSeq)
{
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::CoordType CoordType;

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

template <class MeshType>
void ParametrizeSeq(MeshType &mesh,const std::vector<std::vector<size_t> > &BorderSeq)
{
    for (size_t i=0;i<BorderSeq.size();i++)
        ParametrizeSeq(mesh,BorderSeq[i]);
}

template <class MeshType>
void ComputeUV(MeshType &mesh, ParamType PType,
               bool fixBorders,bool scaleEdges,
               const std::vector<size_t> &CornersIdx,
               std::vector<size_t> &SortedCorner)
{
    if (fixBorders)
    {
        SortedCorner.clear();
        SelectCorners(mesh,CornersIdx);
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


template <class MeshType>
void ComputeParametrizedSubMesh(MeshType &mesh,
                                MeshType &subMesh,
                                const std::vector<size_t> &PatchFaces,
                                const std::vector<size_t> &PatchCorners,
                                std::vector<size_t> &SortedCorners,
                                ParamType PType,
                                bool FixBorders,
                                bool ScaleEdges)
{
    GetMeshFromPatch(mesh,PatchFaces,subMesh);
    std::vector<size_t> CornersIdx;
    GetIndexFromQ(subMesh,PatchCorners,CornersIdx);
    ComputeUV(subMesh,PType,FixBorders,ScaleEdges,CornersIdx,SortedCorners);
}

template <class MeshType>
void SetUVtoPos(MeshType &mesh)
{
    for (size_t i=0;i<mesh.vert.size();i++)
    {
        mesh.vert[i].P().X()=mesh.vert[i].T().P().X();
        mesh.vert[i].P().Y()=mesh.vert[i].T().P().Y();
        mesh.vert[i].P().Z()=0;
    }
}

template <class MeshType>
void GetMeshFromPatch(MeshType &mesh,
                      const std::vector<size_t> &Partition,
                      MeshType &PatchMesh)
{
    for (size_t i=0;i<mesh.vert.size();i++)
        mesh.vert[i].Q()=i;
    for (size_t i=0;i<mesh.face.size();i++)
        mesh.face[i].Q()=i;

    vcg::tri::UpdateSelection<MeshType>::Clear(mesh);

    for (size_t i=0;i<Partition.size();i++)
        mesh.face[Partition[i]].SetS();

    vcg::tri::UpdateSelection<MeshType>::VertexFromFaceLoose(mesh);

    PatchMesh.Clear();
    vcg::tri::Append<MeshType,MeshType>::Mesh(PatchMesh,mesh,true);
    PatchMesh.UpdateAttributes();
}

template <class MeshType>
void GetMeshFromPatch(MeshType &mesh,
                      const size_t &IndexPatch,
                      const std::vector<std::vector<size_t> > &Partitions,
                      MeshType &PatchMesh)
{
    GetMeshFromPatch(mesh,Partitions[IndexPatch],PatchMesh);
}

template <class MeshType>
void GetUVOutline(MeshType &testM,
                  std::vector<vcg::Point2d > &Poly2D,
                  vcg::Box2d &box2D)
{
    vcg::tri::UpdateTopology<MeshType>::FaceFace(testM);
    vcg::tri::UpdateFlags<MeshType>::VertexBorderFromFaceAdj(testM);
    box2D.SetNull();
    for (size_t i=0;i<testM.vert.size();i++)
    {
        if (!testM.vert[i].IsB())continue;
        Poly2D.push_back(testM.vert[i].T().P());
        box2D.Add(testM.vert[i].T().P());
    }
}

template <class MeshType>
void ArrangeUVPatches(std::vector<MeshType*> &ParamPatches)
{
    typedef typename MeshType::ScalarType ScalarType;
    std::vector<std::vector<vcg::Point2<ScalarType> > > ParaOutlines;
    ParaOutlines.resize(ParamPatches.size());

    vcg::Box2d box2D;
    ScalarType AreaUV=0;
    for (size_t i=0;i<ParamPatches.size();i++)
    {
        GetUVOutline<MeshType>(*ParamPatches[i],ParaOutlines[i],box2D);
        AreaUV+=(box2D.DimX()*box2D.DimY());
    }
    int EdgeSize=floor(sqrt(AreaUV)+0.5);
    EdgeSize=std::max(EdgeSize,1);
    vcg::Point2i siz(EdgeSize*2,EdgeSize*2);
    vcg::Point2<ScalarType> coveredA;
    std::vector<vcg::Similarity2<ScalarType> > trVec;
    vcg::PolyPacker<ScalarType>::PackAsObjectOrientedRect(ParaOutlines,siz,trVec,coveredA);

    for (size_t i=0;i<ParamPatches.size();i++)
    {
        for (size_t j=0;j<(*ParamPatches[i]).vert.size();j++)
        {
            vcg::Point2<ScalarType> UVVert=(*ParamPatches[i]).vert[j].T().P();
            (*ParamPatches[i]).vert[j].T().P()=trVec[i]*UVVert;
        }
    }
}


template <class MeshType>
void SelectMeshPatchBorders(MeshType &mesh,const std::vector<int>  &FacePatches)
{
    //std::set<std::pair<size_t,size_t> > BorderPatches;
    vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(mesh);
    assert(FacePatches.size()==mesh.face.size());
    //first add borders
    for (size_t i=0;i<mesh.face.size();i++)
        for (size_t j=0;j<3;j++)
        {
            if (mesh.face[i].IsB(j))
                mesh.face[i].SetFaceEdgeS(j);
            else
            {
                size_t FOpp=vcg::tri::Index(mesh,mesh.face[i].FFp(j));
                assert(FOpp!=i);
                assert(FacePatches[i]>=0);

                if (FacePatches[i]!=FacePatches[FOpp])
                    mesh.face[i].SetFaceEdgeS(j);
            }
        }
}

template <class MeshType>
void SelectTJunctionVertices(MeshType &mesh,const std::vector<int>  &FacePatches)
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

template <class MeshType>
void SmoothMeshPatches(MeshType &mesh,
                       const std::vector<int>  &FacePatches,
                       size_t Steps=3,
                       typename MeshType::ScalarType Damp=0.5)
{
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;

    MeshType TargetMesh;
    vcg::tri::Append<MeshType,MeshType>::Mesh(TargetMesh,mesh);
    TargetMesh.UpdateAttributes();
    vcg::GridStaticPtr<FaceType,ScalarType> Gr;
    Gr.Set(TargetMesh.face.begin(),TargetMesh.face.end());

    for (size_t i=0;i<mesh.vert.size();i++)
        if (mesh.vert[i].IsB())mesh.vert[i].SetS();

    SelectTJunctionVertices(mesh,FacePatches);

    SelectMeshPatchBorders(mesh,FacePatches);

    for (size_t s=0;s<Steps;s++)
    {
        std::vector<CoordType> AvPos(mesh.vert.size(),CoordType(0,0,0));
        std::vector<size_t> NumDiv(mesh.vert.size(),0);

        //smooth path
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!mesh.face[i].IsFaceEdgeS(j))continue;
                size_t VInd0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                size_t VInd1=vcg::tri::Index(mesh,mesh.face[i].V1(j));

                CoordType Pos0=mesh.vert[VInd0].P();
                CoordType Pos1=mesh.vert[VInd1].P();

                AvPos[VInd0]+=Pos1;
                NumDiv[VInd0]++;
                AvPos[VInd1]+=Pos0;
                NumDiv[VInd1]++;
            }

        vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            //fixed
            if (mesh.vert[i].IsS())continue;
            //no contributes
            if (NumDiv[i]==0)continue;

            CoordType TargetPos=(AvPos[i]/(ScalarType)NumDiv[i]);
            mesh.vert[i].P()=(mesh.vert[i].P()*Damp) +
                    ((TargetPos)*(1-Damp));

            mesh.vert[i].SetV();
        }

        //smooth the rest
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
                CoordType Pos0=mesh.vert[IndexV0].P();
                CoordType Pos1=mesh.vert[IndexV1].P();
                if (!mesh.vert[IndexV0].IsV())
                {
                    AvPos[IndexV0]+=Pos1;
                    NumDiv[IndexV0]++;
                }
                if (!mesh.vert[IndexV1].IsV())
                {
                    AvPos[IndexV1]+=Pos0;
                    NumDiv[IndexV1]++;
                }
            }
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            //fixed
            if (mesh.vert[i].IsS())continue;
            if (mesh.vert[i].IsV())continue;
            //no contributes
            if (NumDiv[i]==0)continue;
            mesh.vert[i].P()=(mesh.vert[i].P()*Damp) +
                    ((AvPos[i]/NumDiv[i])*(1-Damp));
        }

        //reproject everything
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            //fixed
            CoordType Pos;
            ScalarType MinD;
            vcg::tri::GetClosestFaceBase(TargetMesh,Gr,
                                         mesh.vert[i].P(),
                                         mesh.bbox.Diag(),MinD,Pos);
            if (mesh.vert[i].IsS())continue;
            mesh.vert[i].P()=Pos;
        }
    }
    vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);
    mesh.UpdateAttributes();
}

template <class MeshType>
void SmoothMeshPatches(MeshType &mesh,
                       const std::vector<std::vector<size_t> >  &PatchFaces,
                       size_t Steps=3,
                       typename MeshType::ScalarType Damp=0.5)
{
    std::vector<int>  FacePatches(mesh.face.size(),-1);
    for (size_t i=0;i<PatchFaces.size();i++)
        for (size_t j=0;j<PatchFaces[i].size();j++)
        {
            size_t CurrF=PatchFaces[i][j];
            FacePatches[CurrF]=i;
        }
    SmoothMeshPatches(mesh,FacePatches,Steps,Damp);
}

//template <class MeshType>
//void PatchesSideLenght(MeshType &mesh,
//                       const std::vector<std::vector<size_t> > &PatchFaces,
//                       const std::vector<std::vector<size_t> > &PatchCorners,
//                       std::vector<std::vector<typename MeshType::ScalarType> > &CurvedL,
//                       std::vector<std::vector<typename MeshType::ScalarType> > &EuclL,
//                       int SmoothNum=10)
//{
//    typedef typename MeshType::ScalarType ScalarType;
//    typedef typename MeshType::CoordType CoordType;
//    std::vector<CoordType> OldPos;

//    if (SmoothNum>0)
//    {
//        for (size_t i=0;i<mesh.vert.size();i++)
//            OldPos.push_back(mesh.vert[i].P());

//        SmoothMeshPatches(mesh,PatchFaces,SmoothNum);
//    }

//    //then get the border sequence
//    std::vector<std::vector<std::vector<size_t> > > BorderSeq;
//    DerivePatchBorderSeq(mesh,PatchFaces,PatchCorners,BorderSeq);

//    CurvedL.resize(BorderSeq.size());
//    EuclL.resize(BorderSeq.size());

//    for (size_t i=0;i<BorderSeq.size();i++)
//        for (size_t j=0;j<BorderSeq[i].size();j++)
//        {
//            CurvedL[i].push_back(Lenght(mesh,BorderSeq[i][j]));
//            CoordType Pos0=mesh.vert[BorderSeq[i][j][0]].P();
//            CoordType Pos1=mesh.vert[BorderSeq[i][j].back()].P();
//            EuclL[i].push_back((Pos1-Pos0).Norm());
//        }

//    if (SmoothNum>0)
//    {
//        for (size_t i=0;i<mesh.vert.size();i++)
//            mesh.vert[i].P()=OldPos[i];
//    }
//}

template <class MeshType>
void PatchesSideLenght(MeshType &mesh,
                       const std::vector<std::vector<size_t> > &PatchFaces,
                       const std::vector<std::vector<size_t> > &PatchCorners,
                       std::vector<std::vector<typename MeshType::ScalarType> > &CurvedL,
                       std::vector<std::vector<typename MeshType::ScalarType> > &EuclL,
                       std::map<std::pair<size_t,size_t>,typename MeshType::ScalarType> &EdgeMap)
{
    typedef typename MeshType::CoordType CoordType;


    //then get the border sequence
    std::vector<std::vector<std::vector<size_t> > > BorderSeq;
    DerivePatchBorderSeq(mesh,PatchFaces,PatchCorners,BorderSeq);

    CurvedL.resize(BorderSeq.size());
    EuclL.resize(BorderSeq.size());

    for (size_t i=0;i<BorderSeq.size();i++)
        for (size_t j=0;j<BorderSeq[i].size();j++)
        {
            CurvedL[i].push_back(FieldLenght(mesh,BorderSeq[i][j],EdgeMap));
            CoordType Pos0=mesh.vert[BorderSeq[i][j][0]].P();
            CoordType Pos1=mesh.vert[BorderSeq[i][j].back()].P();
            EuclL[i].push_back((Pos1-Pos0).Norm());
        }

}

template <class MeshType>
void PatchesLenghtRatios(MeshType &mesh,
                         const std::vector<std::vector<size_t> > &PatchFaces,
                         const std::vector<std::vector<size_t> > &PatchCorners,
                         std::vector<typename MeshType::ScalarType> &Variance,
                         std::vector<typename MeshType::ScalarType> &LenghtDist,
                         std::map<std::pair<size_t,size_t>,typename MeshType::ScalarType> &EdgeMap)
{
    typedef typename MeshType::ScalarType ScalarType;
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

template <class MeshType>
void ColorByVarianceLenght(MeshType &mesh,
                           const std::vector<std::vector<size_t> > &PatchFaces,
                           const std::vector<std::vector<size_t> > &PatchCorners,
                           std::map<std::pair<size_t,size_t>,typename MeshType::ScalarType> &EdgeMap)
{
    typedef typename MeshType::ScalarType ScalarType;
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

template <class MeshType>
void ColorByDistortionLenght(MeshType &mesh,
                             const std::vector<std::vector<size_t> > &PatchFaces,
                             const std::vector<std::vector<size_t> > &PatchCorners,
                             std::map<std::pair<size_t,size_t>,typename MeshType::ScalarType> &EdgeMap)
{
    typedef typename MeshType::ScalarType ScalarType;
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
    std::pair<ScalarType, ScalarType> minmax = Stat<MeshType>::ComputePerFaceQualityMinMax(mesh);
    std::cout<<"Min Dist: "<<minmax.first<<std::endl;
    std::cout<<"Max Dist: "<<minmax.second<<std::endl;
    vcg::tri::UpdateColor<MeshType>::PerFaceQualityRamp(mesh,minmax.second,minmax.first);
}

template <class MeshType>
void ColorByCatmullClarkability(MeshType &mesh,
                                const std::vector<std::vector<size_t> > &PatchFaces,
                                const std::vector<std::vector<size_t> > &PatchCorners,
                                std::map<std::pair<size_t,size_t>,typename MeshType::ScalarType> &EdgeMap,
                                typename MeshType::ScalarType ScaleVal)
{
    typedef typename MeshType::ScalarType ScalarType;
    std::vector<std::vector<ScalarType> > SideL,EuclL;
    PatchesSideLenght(mesh,PatchFaces,PatchCorners,SideL,EuclL,EdgeMap);

    for (size_t i=0;i<SideL.size();i++)
    {
        ScalarType CC=CatmullClarkability(SideL[i]);

        for (size_t j=0;j<PatchFaces[i].size();j++)
        {
            size_t IndexF=PatchFaces[i][j];
            mesh.face[IndexF].Q()=CC;
        }
    }
    //    std::pair<ScalarType, ScalarType> minmax = Stat<MeshType>::ComputePerFaceQualityMinMax(mesh);
    //    std::cout<<"Min CC: "<<minmax.first<<std::endl;
    //    std::cout<<"Max CC: "<<minmax.second<<std::endl;

    std::cout<<"TEST CC"<<std::endl;
    ScalarType CC=CatmullClarkability(SideL[11]);
    for (size_t i=0;i<SideL[11].size();i++)
        std::cout<<"Value L: "<<SideL[11][i]<<std::endl;

    std::cout<<"Value CC: "<<CC<<std::endl;

    std::pair<ScalarType, ScalarType> minmax(0,ScaleVal);
    if (ScaleVal<0)ScaleVal=3;
    vcg::tri::UpdateColor<MeshType>::PerFaceQualityRamp(mesh,minmax.second,minmax.first);

//    for (size_t i=11;i<12;i++)
//    {
//    size_t i=11;
//    for (size_t j=0;j<PatchFaces[i].size();j++)
//    {
//        size_t IndexF=PatchFaces[i][j];
//        mesh.face[IndexF].C()=vcg::Color4b::Gray;
//    }
 //   }
}

template <class MeshType>
void ParametrizePatches(MeshType &mesh,
                        MeshType &splittedUV,
                        std::vector<std::vector<size_t> > &PatchFaces,
                        std::vector<std::vector<size_t> > &PatchCorners,
                        ParamType PType,
                        bool FixBorders,
                        bool ScaleEdges,
                        bool Arrange=true,
                        bool Save=false)
{
    std::vector<MeshType*> ParamPatches;
    for (size_t i=0;i<PatchFaces.size();i++)
    {
        ParamPatches.push_back(new MeshType());
        std::vector<size_t> SortedCorners;
        ComputeParametrizedSubMesh(mesh,*ParamPatches.back(),
                                   PatchFaces[i],PatchCorners[i],
                                   SortedCorners,PType,
                                   FixBorders,ScaleEdges);
    }
    if (Arrange)
        ArrangeUVPatches(ParamPatches);

    splittedUV.Clear();
    for (size_t i=0;i<PatchFaces.size();i++)
    {
        vcg::tri::Append<MeshType,MeshType>::Mesh(splittedUV,*ParamPatches[i]);
        delete(ParamPatches[i]);
    }
    if (Save)
    {
        SetUVtoPos(splittedUV);
        vcg::tri::io::ExporterPLY<MeshType>::Save(splittedUV,"parametrize.ply",vcg::tri::io::Mask::IOM_FACECOLOR);
    }
}


template <class MeshType>
void ColorByUVDistortionFaces(MeshType &mesh,
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

template <class MeshType>
size_t NumEmitters(MeshType &mesh,
                   std::vector<size_t> &PatchFaces,
                   std::vector<size_t> &VerticesNeeds)
{
    size_t ret=0;
    vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);
    for (size_t i=0;i<PatchFaces.size();i++)
        for (size_t j=0;j<3;j++)
        {
            size_t IndexV=vcg::tri::Index(mesh,mesh.face[PatchFaces[i]].V(j));
            if (mesh.vert[IndexV].IsV())continue;
            mesh.vert[IndexV].SetV();
            ret+=VerticesNeeds[IndexV];
        }
    return ret;
}

//bool UpdatePatchQuality(MeshType &mesh,
//                  std::vector<std::vector<size_t> > &PatchFaces,
//                  std::vector<std::vector<size_t> > &PatchCorners)
//{

//}

template <class MeshType>
void GetPatchInfo(MeshType &mesh,
                  std::vector<std::vector<size_t> > &PatchFaces,
                  std::vector<std::vector<size_t> > &PatchCorners,
                  std::vector<size_t> &VerticesNeeds,
                  std::map<std::pair<size_t,size_t>,typename MeshType::ScalarType> &EdgeMap,
                  std::vector<PatchInfo<typename MeshType::ScalarType> > &PInfo)
{
    typedef typename MeshType::ScalarType ScalarType;
    PInfo.clear();
    PInfo.resize(PatchFaces.size());
    std::vector<std::vector<typename MeshType::ScalarType> > CurvedL,EuclL;
    PatchesSideLenght(mesh,PatchFaces,PatchCorners,CurvedL,EuclL,EdgeMap);
    for (size_t i=0;i<PatchFaces.size();i++)
    {
        PInfo[i].NumCorners=PatchCorners[i].size();
        PInfo[i].NumEmitters=NumEmitters(mesh,PatchFaces[i],VerticesNeeds);
        PInfo[i].Genus=PatchGenus(mesh,PatchFaces[i]);
        PInfo[i].CurvedL=CurvedL[i];
        PInfo[i].CornerL=EuclL[i];

        if ((PInfo[i].NumCorners<3)||(PInfo[i].NumCorners>6))
        {
            PInfo[i].CClarkability=std::numeric_limits<ScalarType>::max();

        }
        else
        {
            PInfo[i].CClarkability=CatmullClarkability(CurvedL[i]);
        }
    }
}


template <class MeshType>
typename MeshType::ScalarType PatchArea(MeshType &mesh,const std::vector<size_t> &PatchFaces)
{
    MeshType patch_mesh;
    GetMeshFromPatch(mesh,PatchFaces,patch_mesh);
    return(MeshArea(patch_mesh));
}

//template <class ScalarType>
//bool IsOKPatch(const PatchInfo<ScalarType> &PInfo,
//               size_t MinSides,
//               size_t MaxSides,
//               ScalarType CClarkability,
//               std::vector<ScalarType> &QThresold)
//{
//    if (PInfo.Genus!=1)return false;
//    if (PInfo.NumEmitters>0)return false;
//    if (PInfo.NumCorners<MinSides)return false;
//    if (PInfo.NumCorners>MaxSides)return false;
//    if ((CClarkability>0)&&(PInfo.CClarkability>CClarkability))return false;
//    for (size_t i=0;i<QThresold.size();i++)
//        if (PInfo.Q[i]<QThresold[i])return false;
//    return true;
//}

template <class MeshType>
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

//template <class MeshType>
//size_t NonOkPartitions(MeshType &mesh,
//                       const std::vector<PatchInfo<typename MeshType::ScalarType> > &PInf,
//                       size_t MinSides,
//                       size_t MaxSides,
//                       typename MeshType::ScalarType CClarkability,
//                       std::vector<typename MeshType::ScalarType> &QThresold)
//{
//    size_t CurrN=0;
//    for (size_t i=0;i<PInf.size();i++)
//    {
//        if (IsOKPatch(PInf[i],MinSides,MaxSides,CClarkability,QThresold))continue;
//        CurrN++;
//    }
//    return CurrN;
//}

template <class MeshType>
bool BetterConfiguaration(MeshType &mesh,
                          const std::vector<std::vector<size_t> > &PatchFaces0,
                          const std::vector<std::vector<size_t> > &PatchFaces1,
                          const std::vector<PatchInfo<typename MeshType::ScalarType> > &PInf0,
                          const std::vector<PatchInfo<typename MeshType::ScalarType> > &PInf1,
                          size_t MinSides,size_t MaxSides,
                          typename MeshType::ScalarType CClarkability,
                          std::vector<typename MeshType::ScalarType> &QThresold)
{
    size_t NonOKGenus0=0;
    size_t NonOKGenus1=0;
    size_t NonOKEmitters0=0;
    size_t NonOKEmitters1=0;
    size_t NonOKSize0=0;
    size_t NonOKSize1=0;
//    size_t NonOKCC0=0;
//    size_t NonOKCC1=0;
    size_t Sing0=0;
    size_t Sing1=0;
    for (size_t i=0;i<PInf0.size();i++)
    {
        if (PInf0[i].Genus!=1)NonOKGenus0++;
        if (PInf0[i].NumEmitters>0)NonOKEmitters0++;
        if (PInf0[i].NumCorners<(int)MinSides)NonOKSize0++;
        if (PInf0[i].NumCorners>(int)MaxSides)NonOKSize0++;
        if (CClarkability>1)
            Sing0=AddedSingularities(PInf0[i].CurvedL,CClarkability);
        //if ((CClarkability>0)&&(PInf0[i].CClarkability>CClarkability))NonOKCC0++;
    }

    for (size_t i=0;i<PInf1.size();i++)
    {
        if (PInf1[i].Genus!=1)NonOKGenus1++;
        if (PInf1[i].NumEmitters>0)NonOKEmitters1++;
        if (PInf1[i].NumCorners<(int)MinSides)NonOKSize1++;
        if (PInf1[i].NumCorners>(int)MaxSides)NonOKSize1++;
        if (CClarkability>1)
            Sing1=AddedSingularities(PInf1[i].CurvedL,CClarkability);
        //if ((CClarkability>0)&&(PInf1[i].CClarkability>CClarkability))NonOKCC1++;
    }
    if (NonOKGenus1!=NonOKGenus0)return (NonOKGenus1<NonOKGenus0);
    if (NonOKEmitters1!=NonOKEmitters0)return (NonOKEmitters1<NonOKEmitters0);
    if (NonOKSize1!=NonOKSize0)return (NonOKSize1<NonOKSize0);
    //if (NonOKCC1!=NonOKCC0)return (NonOKCC1<NonOKCC0);
    if (Sing0!=Sing1)return (Sing1<Sing0);
    return true;
    //    size_t Num0=NonOkPartitions(mesh,PInf0,MinSides,MaxSides,CatmullClarkability,QThresold);
    //    size_t Num1=NonOkPartitions(mesh,PInf1,MinSides,MaxSides,CatmullClarkability,QThresold);
    //    return (Num1<=Num0);
    //    typename MeshType::ScalarType Area0=NonOkArea(mesh,PatchFaces0,PInf0,MinSides,MaxSides,QThresold);
    //    typename MeshType::ScalarType Area1=NonOkArea(mesh,PatchFaces1,PInf1,MinSides,MaxSides,QThresold);
    //    std::cout<<"Area 0:"<<Area0;
    //    std::cout<<"Area 1:"<<Area1;
    //    return (Area1<=Area0);
}


#endif
