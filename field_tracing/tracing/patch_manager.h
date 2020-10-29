#ifndef PATCH_MANAGER
#define PATCH_MANAGER

#include <wrap/igl/arap_parametrization.h>
#include <wrap/igl/lscm_parametrization.h>
#include <tracing/candidate_path.h>
#include "vert_field_graph.h"
#include <vcg/complex/algorithms/parametrization/distortion.h>
#include <vcg/space/outline2_packer.h>
#include "catmull_clarkability.h"
#include <vcg/space/triangle3.h>
#include <vcg/complex/algorithms/create/platonic.h>

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
    bool CClarkability;
    int PossibleSing;
    std::vector<ScalarType> Q;
    int ExpectedValence;
    int NumSing;
    PatchInfo()
    {
        NumEmitters=0;
        NumCorners=0;
        Genus=0;
        CClarkability=false;//std::numeric_limits<ScalarType>::max();
        ExpectedValence=-1;
        NumSing=-1;
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
int ExpectedValence(MeshType &mesh,const std::vector<size_t> &PatchFaces)
{
    typedef typename MeshType::VertexType VertexType;
    int ExpVal=4;
    vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);
    for (size_t i=0;i<PatchFaces.size();i++)
    {
        for (size_t j=0;j<3;j++)
        {
            VertexType *v=mesh.face[PatchFaces[i]].V(j);
            if (v->IsV())continue;
            v->SetV();
            if (v->SingularityValence==4)continue;
            if (ExpVal!=4)return -1;//multiple singularities
            ExpVal=v->SingularityValence;
        }
    }
    vcg:tri::UpdateFlags<MeshType>::VertexClearV(mesh);
    return ExpVal;
}

template <class MeshType>
int NumSingularities(MeshType &mesh,const std::vector<size_t> &PatchFaces)
{
    typedef typename MeshType::VertexType VertexType;
    int NumSing=0;
    vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);
    for (size_t i=0;i<PatchFaces.size();i++)
    {
        for (size_t j=0;j<3;j++)
        {
            VertexType *v=mesh.face[PatchFaces[i]].V(j);
            if (v->IsV())continue;
            v->SetV();
            if (v->SingularityValence==4)continue;
            NumSing++;
        }
    }
    vcg:tri::UpdateFlags<MeshType>::VertexClearV(mesh);
    return NumSing;
}

template <class MeshType>
int PatchGenus(MeshType &mesh,const std::vector<size_t> &PatchFaces)
{
    SelectPatchFacesVert<MeshType>(mesh,PatchFaces);

    MeshType subM;
    vcg::tri::Append<MeshType,MeshType>::Mesh(subM,mesh,true);
    vcg::tri::Clean<MeshType>::RemoveDuplicateVertex(subM);
    vcg::tri::Allocator<MeshType>::CompactEveryVector(subM);

    std::set<std::pair<size_t,size_t> > EdgeSet;
    size_t NumF=0;
    size_t NumV=0;
    size_t NumE=0;
    for (size_t i=0;i<subM.face.size();i++)
    {
        if (!subM.face[i].IsS())continue;
        NumF++;
        for (size_t j=0;j<3;j++)
        {
            size_t IndV0=vcg::tri::Index(subM,subM.face[i].V0(j));
            size_t IndV1=vcg::tri::Index(subM,subM.face[i].V1(j));
            EdgeSet.insert(std::pair<size_t,size_t>(std::min(IndV0,IndV1),std::max(IndV0,IndV1)));
        }
    }
    for (size_t i=0;i<subM.vert.size();i++)
    {
        if (subM.vert[i].IsD())continue;
        if (subM.vert[i].IsS())NumV++;
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
    vcg::tri::io::ExporterPLY<MeshType>::Save(mesh,"testBorder.ply",vcg::tri::io::Mask::IOM_FLAGS);
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

template <class MeshType>
void DerivePatchBorderSeq(MeshType &mesh,
                          const std::vector<size_t>  &PatchFaces,
                          const std::vector<size_t> &PatchCorners,
                          std::vector<std::vector<size_t> > &BorderSeq)
{
    MeshType subMesh;
    GetMeshFromPatch(mesh,PatchFaces,subMesh);
    std::vector<size_t> CornersIdx;
    GetIndexFromQ(subMesh,PatchCorners,CornersIdx);


    SelectCorners(subMesh,CornersIdx);
    if (CornersIdx.size()==0)return;
    DeriveBorderSeq(subMesh,BorderSeq);

    //then set to original Index
    for (size_t j=0;j<BorderSeq.size();j++)
        for (size_t k=0;k<BorderSeq[j].size();k++)
            BorderSeq[j][k]=subMesh.vert[BorderSeq[j][k]].Q();
}

template <class MeshType>
void DerivePatchBorderSeq(MeshType &mesh,
                          const std::vector<std::vector<size_t> > &PatchFaces,
                          const std::vector<std::vector<size_t> > &PatchCorners,
                          std::vector<std::vector<std::vector<size_t> > > &BorderSeq)
{
    BorderSeq.clear();
    //std::cout<<"Patches Num "<<PatchFaces.size()<<std::endl;
    BorderSeq.clear();
    BorderSeq.resize(PatchFaces.size());
    for (size_t i=0;i<PatchFaces.size();i++)
    {
        DerivePatchBorderSeq(mesh,PatchFaces[i],PatchCorners[i],BorderSeq[i]);
        //        BorderSeq.resize(BorderSeq.size()+1);
        //        MeshType subMesh;
        //        GetMeshFromPatch(mesh,PatchFaces[i],subMesh);
        //        std::vector<size_t> CornersIdx;
        //        GetIndexFromQ(subMesh,PatchCorners[i],CornersIdx);


        //        SelectCorners(subMesh,CornersIdx);
        //        if (CornersIdx.size()==0)continue;
        //        DeriveBorderSeq(subMesh,BorderSeq[i]);

        //        //then set to original Index
        //        for (size_t j=0;j<BorderSeq[i].size();j++)
        //            for (size_t k=0;k<BorderSeq[i][j].size();k++)
        //                BorderSeq[i][j][k]=subMesh.vert[BorderSeq[i][j][k]].Q();
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
                if (FacePatches[i]<0)assert(FacePatches[i]==-1);
                //assert(FacePatches[i]>=0);

                if (FacePatches[i]!=FacePatches[FOpp])
                    mesh.face[i].SetFaceEdgeS(j);
            }
        }
}

template <class MeshType>
void SelectVertOnMeshPatchBorders(MeshType &mesh,const std::vector<int>  &FacePatches)
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

template <class MeshType>
void DerivePerFacePartition(const MeshType &mesh,
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

template <class MeshType>
void SelectVertOnMeshPatchBorders(MeshType &mesh,const std::vector<std::vector<int> >  &PatchFaces)
{
    std::vector<int> FacePartitions;
    DerivePerFacePartition(mesh,PatchFaces,FacePartitions);
    SelectVertOnMeshPatchBorders(mesh,FacePartitions);
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
int SelectFolds(MeshType &mesh,
                int dilateStep=0,
                typename MeshType::ScalarType MinQ=0.2)
{
    typedef typename MeshType::FaceType FaceType;
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

    for (size_t i=0;i<dilateStep;i++)
    {
        vcg::tri::UpdateSelection<MeshType>::FaceFromVertexLoose(mesh);
        vcg::tri::UpdateSelection<MeshType>::VertexFromFaceLoose(mesh);
    }
    return numF;
}

//template <class MeshType>
//void Reproject(MeshType &mesh,bool fixS=true)
//{
//    typedef typename MeshType::FaceType FaceType;
//    typedef typename MeshType::CoordType CoordType;
//    typedef typename MeshType::ScalarType ScalarType;

//    MeshType TargetMesh;
//    vcg::tri::Append<MeshType,MeshType>::Mesh(TargetMesh,mesh);
//    for (size_t i=0;i<mesh.vert.size();i++)
//        TargetMesh.vert[i].P()=TargetMesh.vert[i].RPos;
//    TargetMesh.UpdateAttributes();
//    vcg::GridStaticPtr<FaceType,ScalarType> Gr;
//    Gr.Set(TargetMesh.face.begin(),TargetMesh.face.end());

//    //reproject everything
//    for (size_t i=0;i<mesh.vert.size();i++)
//    {
//        //fixed
//        CoordType Pos;
//        ScalarType MinD;
//        vcg::tri::GetClosestFaceBase(TargetMesh,Gr,
//                                     mesh.vert[i].P(),
//                                     mesh.bbox.Diag(),MinD,Pos);
//        if ((fixS)&&(mesh.vert[i].IsS()))continue;
//        mesh.vert[i].P()=Pos;
//    }
//}

template <class MeshType>
void SolveFolds(MeshType &mesh,
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

template <class MeshType>
void LaplacianPos(MeshType &mesh,
                  std::vector<typename MeshType::CoordType> &AvPos,
                  bool OnlySContribute)
{
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
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

template <class MeshType>
void LaplacianInternalStep(MeshType &mesh,const std::vector<int>  &FacePatches,
                           typename MeshType::ScalarType Damp=0.5)
{
    SelectVertOnMeshPatchBorders(mesh,FacePatches);
    std::vector<typename MeshType::CoordType> AvPos;
    LaplacianPos(mesh,AvPos,false);
    for (size_t i=0;i<mesh.vert.size();i++)
    {
        if (mesh.vert[i].IsS())continue;
        mesh.vert[i].P()=mesh.vert[i].P()*Damp+AvPos[i]*(1-Damp);
    }
}

template <class MeshType>
void LaplacianBorderStep(MeshType &mesh,
                         const std::vector<int>  &FacePatches,
                         typename MeshType::ScalarType Damp=0.5)
{
    //    SelectVertOnMeshPatchBorders(mesh,FacePatches);
    //    std::vector<typename MeshType::CoordType> AvPos;
    //    LaplacianPos(mesh,AvPos,true);
    //    for (size_t i=0;i<mesh.vert.size();i++)
    //    {
    //        if (!mesh.vert[i].IsS())continue;
    //        mesh.vert[i].P()=mesh.vert[i].P()*Damp+AvPos[i]*(1-Damp);
    //    }
//    for (size_t i=0;i<mesh.face.size();i++)
//        for (size_t j=0;j<3;j++)
//        {
//            if (!mesh.face[i].IsFaceEdgeS(j))continue;
//            std::cout<<i<<","<<j<<std::endl;
//        }
//    vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(mesh);
//    return;

    //save previous selection

    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    SelectMeshPatchBorders(mesh,FacePatches);

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
        if (NumDiv[i]==0)continue;
        CoordType TargetPos=AvPos[i]/(ScalarType)NumDiv[i];
        mesh.vert[i].P()=mesh.vert[i].P()*Damp+TargetPos*(1-Damp);
    }
    vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(mesh);
}

template <class MeshType>
void ReprojectOn(MeshType &mesh,MeshType &target,
                 vcg::GridStaticPtr<typename MeshType::FaceType,typename MeshType::ScalarType> &Gr)
{
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
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

//template <class MeshType>
//void SmoothMeshPatches(MeshType &mesh,
//                       const std::vector<int>  &FacePatches,
//                       size_t Steps=3,
//                       typename MeshType::ScalarType Damp=0.5)
//{
//    typedef typename MeshType::FaceType FaceType;
//    typedef typename MeshType::CoordType CoordType;
//    typedef typename MeshType::ScalarType ScalarType;

//    MeshType TargetMesh;
//    vcg::tri::Append<MeshType,MeshType>::Mesh(TargetMesh,mesh);
//    TargetMesh.UpdateAttributes();
//    vcg::GridStaticPtr<FaceType,ScalarType> Gr;
//    Gr.Set(TargetMesh.face.begin(),TargetMesh.face.end());

//    for (size_t i=0;i<mesh.vert.size();i++)
//        if (mesh.vert[i].IsB())mesh.vert[i].SetS();

//    SelectTJunctionVertices(mesh,FacePatches);

//    SelectMeshPatchBorders(mesh,FacePatches);

//    for (size_t s=0;s<Steps;s++)
//    {
//        std::vector<CoordType> AvPos(mesh.vert.size(),CoordType(0,0,0));
//        std::vector<size_t> NumDiv(mesh.vert.size(),0);

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

//        vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);
//        for (size_t i=0;i<mesh.vert.size();i++)
//        {
//            //fixed
//            if (mesh.vert[i].IsS())continue;
//            //no contributes
//            if (NumDiv[i]==0)continue;

//            CoordType TargetPos=(AvPos[i]/(ScalarType)NumDiv[i]);
//            mesh.vert[i].P()=(mesh.vert[i].P()*Damp) +
//                    ((TargetPos)*(1-Damp));

//            mesh.vert[i].SetV();
//        }

//        //smooth the rest
//        for (size_t i=0;i<mesh.face.size();i++)
//            for (size_t j=0;j<3;j++)
//            {
//                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
//                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
//                CoordType Pos0=mesh.vert[IndexV0].P();
//                CoordType Pos1=mesh.vert[IndexV1].P();
//                if (!mesh.vert[IndexV0].IsV())
//                {
//                    AvPos[IndexV0]+=Pos1;
//                    NumDiv[IndexV0]++;
//                }
//                if (!mesh.vert[IndexV1].IsV())
//                {
//                    AvPos[IndexV1]+=Pos0;
//                    NumDiv[IndexV1]++;
//                }
//            }
//        for (size_t i=0;i<mesh.vert.size();i++)
//        {
//            //fixed
//            if (mesh.vert[i].IsS())continue;
//            if (mesh.vert[i].IsV())continue;
//            //no contributes
//            if (NumDiv[i]==0)continue;
//            mesh.vert[i].P()=(mesh.vert[i].P()*Damp) +
//                    ((AvPos[i]/NumDiv[i])*(1-Damp));
//        }

//        //reproject everything
//        for (size_t i=0;i<mesh.vert.size();i++)
//        {
//            //fixed
//            CoordType Pos;
//            ScalarType MinD;
//            vcg::tri::GetClosestFaceBase(TargetMesh,Gr,
//                                         mesh.vert[i].P(),
//                                         mesh.bbox.Diag(),MinD,Pos);
//            if (mesh.vert[i].IsS())continue;
//            mesh.vert[i].P()=Pos;
//        }
//        //Reproject(mesh,true);
//    }
//    vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);
//    mesh.UpdateAttributes();

//    SolveFolds(mesh);
//}

template <class MeshType>
void SaveEdgeSel(MeshType &mesh,std::vector<std::vector<bool> > &EdgeSel)
{
   EdgeSel.resize(mesh.face.size(),std::vector<bool>(3,false));
   for (size_t i=0;i<mesh.face.size();i++)
       for (size_t j=0;j<3;j++)
       {
            if (!mesh.face[i].IsFaceEdgeS(j))continue;
            EdgeSel[i][j]=true;
       }
}


template <class MeshType>
void RestoreEdgeSel(MeshType &mesh,
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

template <class MeshType>
void SmoothMeshPatches(MeshType &mesh,
                       const std::vector<int>  &FacePatches,
                       size_t Steps=3,
                       typename MeshType::ScalarType Damp=0.5,
                       typename MeshType::ScalarType MinQ=0.2)
{
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::ScalarType ScalarType;

    std::vector<std::vector<bool> > EdgeSel;
    SaveEdgeSel(mesh,EdgeSel);

    MeshType TargetMesh;
    vcg::tri::Append<MeshType,MeshType>::Mesh(TargetMesh,mesh);
    TargetMesh.UpdateAttributes();
    vcg::GridStaticPtr<FaceType,ScalarType> Gr;
    Gr.Set(TargetMesh.face.begin(),TargetMesh.face.end());

    for (size_t s=0;s<Steps;s++)
    {
        //save old position
        std::vector<CoordType> OldPos;
        if (MinQ>0)
        for (size_t i=0;i<mesh.vert.size();i++)
            OldPos.push_back(mesh.vert[i].P());

        //save old normals
        vcg::tri::UpdateNormal<MeshType>::PerFaceNormalized(mesh);
        std::vector<CoordType> OldNorm;
        if (MinQ>0)
        for (size_t i=0;i<mesh.face.size();i++)
            OldNorm.push_back(mesh.face[i].N());

        //PERFORM SMOOTHING
        LaplacianBorderStep(mesh,FacePatches,Damp);
        ReprojectOn(mesh,TargetMesh,Gr);
        LaplacianInternalStep(mesh,FacePatches,Damp);
        ReprojectOn(mesh,TargetMesh,Gr);
        if (MinQ<=0)continue;

        //save new position
        std::vector<CoordType> NewPos;
        for (size_t i=0;i<mesh.vert.size();i++)
            NewPos.push_back(mesh.vert[i].P());

        //restore old position
        for (size_t i=0;i<mesh.vert.size();i++)
            mesh.vert[i].P()=OldPos[i];

        for (size_t i=0;i<mesh.vert.size();i++)
        {
            //if border not change
            if (mesh.vert[i].IsB())continue;
            //try new position
            mesh.vert[i].P()=NewPos[i];
            std::vector<FaceType*> FaceVec;
            std::vector<int> Indexes;
            vcg::face::VFStarVF(&mesh.vert[i],FaceVec,Indexes);

            bool IsOk=true;
            for (size_t i=0;i<FaceVec.size();i++)
            {
                size_t IndexF=vcg::tri::Index(mesh,FaceVec[i]);
                CoordType P0=FaceVec[i]->P(0);
                CoordType P1=FaceVec[i]->P(1);
                CoordType P2=FaceVec[i]->P(2);
                CoordType NewNormF=vcg::Normal(P0,P1,P2);
                NewNormF.Normalize();
                CoordType OldNormF=OldNorm[IndexF];
                if ((NewNormF*OldNormF)<0){IsOk=false;break;}
                ScalarType QFace=vcg::QualityRadii(P0,P1,P2);
                if (QFace<MinQ){IsOk=false;break;}
            }
            //restore if not ok
            if (!IsOk)
                mesh.vert[i].P()=OldPos[i];
        }
    }
    vcg::tri::UpdateSelection<MeshType>::Clear(mesh);
//    if (MinQ>0)
//        SolveFolds(mesh,MinQ);

    RestoreEdgeSel(mesh,EdgeSel);
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
                       const std::vector<size_t>  &PatchFaces,
                       const std::vector<size_t>  &PatchCorners,
                       std::vector<typename MeshType::ScalarType>  &CurvedL,
                       std::vector<typename MeshType::ScalarType>  &EuclL,
                       std::map<std::pair<size_t,size_t>,typename MeshType::ScalarType> &EdgeMap)
{
    typedef typename MeshType::CoordType CoordType;


    //then get the border sequence
    std::vector<std::vector<size_t> > BorderSeq;
    DerivePatchBorderSeq(mesh,PatchFaces,PatchCorners,BorderSeq);
    for (size_t j=0;j<BorderSeq.size();j++)
    {
        CurvedL.push_back(FieldLenght(mesh,BorderSeq[j],EdgeMap));
        CoordType Pos0=mesh.vert[BorderSeq[j][0]].P();
        CoordType Pos1=mesh.vert[BorderSeq[j].back()].P();
        EuclL.push_back((Pos1-Pos0).Norm());
    }
}

template <class MeshType>
void PatchesSideLenght(MeshType &mesh,
                       const std::vector<std::vector<size_t> > &PatchFaces,
                       const std::vector<std::vector<size_t> > &PatchCorners,
                       std::vector<std::vector<typename MeshType::ScalarType> > &CurvedL,
                       std::vector<std::vector<typename MeshType::ScalarType> > &EuclL,
                       std::map<std::pair<size_t,size_t>,typename MeshType::ScalarType> &EdgeMap)
{
    //    typedef typename MeshType::CoordType CoordType;


    //    //then get the border sequence
    //    std::vector<std::vector<std::vector<size_t> > > BorderSeq;
    //    DerivePatchBorderSeq(mesh,PatchFaces,PatchCorners,BorderSeq);

    CurvedL.clear();
    EuclL.clear();
    CurvedL.resize(PatchFaces.size());
    EuclL.resize(PatchFaces.size());

    for (size_t i=0;i<PatchFaces.size();i++)
    {
        PatchesSideLenght(mesh,PatchFaces[i],PatchCorners[i],CurvedL[i],EuclL[i],EdgeMap);
    }
    //        for (size_t j=0;j<BorderSeq[i].size();j++)
    //        {
    //            CurvedL[i].push_back(FieldLenght(mesh,BorderSeq[i][j],EdgeMap));
    //            CoordType Pos0=mesh.vert[BorderSeq[i][j][0]].P();
    //            CoordType Pos1=mesh.vert[BorderSeq[i][j].back()].P();
    //            EuclL[i].push_back((Pos1-Pos0).Norm());
    //        }
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
                                const typename MeshType::ScalarType ScaleVal,
                                const typename MeshType::ScalarType Thr)
{
    typedef typename MeshType::ScalarType ScalarType;
    std::vector<std::vector<ScalarType> > SideL,EuclL;
    PatchesSideLenght(mesh,PatchFaces,PatchCorners,SideL,EuclL,EdgeMap);
    for (size_t i=0;i<SideL.size();i++)
    {
        bool CC=IsCatmullClarkable(PatchFaces[i].size(),SideL[i],Thr);

        for (size_t j=0;j<PatchFaces[i].size();j++)
        {
            size_t IndexF=PatchFaces[i][j];
            if (CC)
                mesh.face[IndexF].C()=vcg::Color4b::Green;
            else
                mesh.face[IndexF].C()=vcg::Color4b::Red;
            //            if (CC>0)CC=1;
            //            mesh.face[IndexF].Q()=CC;
            //            if (CC>ScaleVal)CC=ScaleVal;
            //            mesh.face[IndexF].Q()=CC;
        }
    }
    //    std::pair<ScalarType, ScalarType> minmax = Stat<MeshType>::ComputePerFaceQualityMinMax(mesh);
    //    std::cout<<"Min CC: "<<minmax.first<<std::endl;
    //    std::cout<<"Max CC: "<<minmax.second<<std::endl;

    //std::cout<<"TEST CC"<<std::endl;
    //    ScalarType CC=CatmullClarkability(SideL[11]);
    //    for (size_t i=0;i<SideL[11].size();i++)
    //        std::cout<<"Value L: "<<SideL[11][i]<<std::endl;

    //    std::cout<<"Value CC: "<<CC<<std::endl;

    //std::pair<ScalarType, ScalarType> minmax(1,ScaleVal);
    //if (ScaleVal<0)ScaleVal=3;
    //    vcg::tri::UpdateColor<MeshType>::PerFaceQualityRamp(mesh,1,0);

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

//template <class MeshType>
//void GetPatchInfo(MeshType &mesh,
//                  std::vector<std::vector<size_t> > &PatchFaces,
//                  std::vector<std::vector<size_t> > &PatchCorners,
//                  std::vector<size_t> &VerticesNeeds,
//                  std::map<std::pair<size_t,size_t>,typename MeshType::ScalarType> &EdgeMap,
//                  std::vector<PatchInfo<typename MeshType::ScalarType> > &PInfo,
//                  const typename MeshType::ScalarType Thr)
//{
//    typedef typename MeshType::ScalarType ScalarType;
//    PInfo.clear();
//    PInfo.resize(PatchFaces.size());
//    std::vector<std::vector<typename MeshType::ScalarType> > CurvedL,EuclL;

//    PatchesSideLenght(mesh,PatchFaces,PatchCorners,CurvedL,EuclL,EdgeMap);

//    for (size_t i=0;i<PatchFaces.size();i++)
//    {
//        PInfo[i].NumCorners=PatchCorners[i].size();
//        PInfo[i].NumEmitters=NumEmitters(mesh,PatchFaces[i],VerticesNeeds);
//        PInfo[i].Genus=PatchGenus(mesh,PatchFaces[i]);
//        PInfo[i].CurvedL=CurvedL[i];
//        PInfo[i].CornerL=EuclL[i];

//        if ((PInfo[i].NumCorners<3)||(PInfo[i].NumCorners>6))
//        {
//            PInfo[i].CClarkability=false;
//        }
//        else
//        {
//            PInfo[i].CClarkability=IsCatmullClarkable(CurvedL[i],Thr);
//        }
//    }

//}


template <class MeshType>
void GetPatchInfo(MeshType &mesh,
                  std::vector<std::vector<size_t> > &PatchFaces,
                  std::vector<std::vector<size_t> > &PatchCorners,
                  std::vector<size_t> &VerticesNeeds,
                  std::map<std::pair<size_t,size_t>,typename MeshType::ScalarType> &EdgeMap,
                  std::vector<PatchInfo<typename MeshType::ScalarType> > &PInfo,
                  const typename MeshType::ScalarType Thr)
{
    typedef typename MeshType::ScalarType ScalarType;
    PInfo.clear();
    PInfo.resize(PatchFaces.size());

    //std::vector<std::vector<typename MeshType::ScalarType> > CurvedL,EuclL;

    //PatchesSideLenght(mesh,PatchFaces,PatchCorners,CurvedL,EuclL,EdgeMap);

    for (size_t i=0;i<PatchFaces.size();i++)
    {
        PInfo[i].NumCorners=PatchCorners[i].size();
        PInfo[i].NumEmitters=NumEmitters(mesh,PatchFaces[i],VerticesNeeds);
        PInfo[i].Genus=PatchGenus(mesh,PatchFaces[i]);
        PInfo[i].ExpectedValence=ExpectedValence(mesh,PatchFaces[i]);
        PInfo[i].NumSing=NumSingularities(mesh,PatchFaces[i]);
        //NumSing
        //PInfo[i].CurvedL=CurvedL[i];
        //PInfo[i].CornerL=EuclL[i];

        if ((PInfo[i].NumCorners<3)||
                (PInfo[i].NumCorners>6)||
                (PInfo[i].Genus!=1)||
                (PInfo[i].NumEmitters>0)
                )
        {
            PInfo[i].CClarkability=false;
        }
        else
        {
            std::vector<typename MeshType::ScalarType > CurvedL,EuclL;
            PatchesSideLenght(mesh,PatchFaces[i],PatchCorners[i],CurvedL,EuclL,EdgeMap);
            PInfo[i].CurvedL=CurvedL;
            PInfo[i].CornerL=EuclL;
            PInfo[i].CClarkability=IsCatmullClarkable(PatchFaces[i].size(),CurvedL,Thr);
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
                          const typename MeshType::ScalarType CClarkability,
                          const typename MeshType::ScalarType avgEdge,
                          bool match_sing_valence)
{
    typedef typename MeshType::ScalarType ScalarType;
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
//    size_t NonOkVal0=0;
//    size_t NonOkVal1=0;
    //at least one sing inside
    int MatchSing0=0;
    int MatchSing1=0;
    int MaxInternalSing0=1;
    int MaxInternalSing1=1;
    for (size_t i=0;i<PInf0.size();i++)
    {
        if (PInf0[i].Genus!=1)NonOKGenus0++;
        if (PInf0[i].NumEmitters>0)NonOKEmitters0++;
        if (PInf0[i].NumCorners<(int)MinSides)NonOKSize0++;
        if (PInf0[i].NumCorners>(int)MaxSides)NonOKSize0++;
        if (match_sing_valence)
        {
            MaxInternalSing0=std::max(MaxInternalSing0,PInf0[i].NumSing);
//            if (PInf0[i].NumCorners!=(int)PInf0[i].ExpectedValence)
//                NonOkVal0++;
            if ((PInf0[i].ExpectedValence!=4)&&
                (PInf0[i].NumCorners==(int)PInf0[i].ExpectedValence))
                MatchSing0++;
        }
        if (CClarkability>0)
            Sing0+=AddedSingularities(PatchFaces0[i].size(),PInf0[i].CurvedL,CClarkability*avgEdge);
        //if ((CClarkability>0)&&(PInf0[i].CClarkability>CClarkability))NonOKCC0++;
    }

    for (size_t i=0;i<PInf1.size();i++)
    {
        if (PInf1[i].Genus!=1)NonOKGenus1++;
        if (PInf1[i].NumEmitters>0)NonOKEmitters1++;
        if (PInf1[i].NumCorners<(int)MinSides)NonOKSize1++;
        if (PInf1[i].NumCorners>(int)MaxSides)NonOKSize1++;
        if (match_sing_valence)
        {
            MaxInternalSing1=std::max(MaxInternalSing1,PInf1[i].NumSing);
//            if (PInf1[i].NumCorners!=(int)PInf1[i].ExpectedValence)
//                NonOkVal1++;
            if ((PInf1[i].ExpectedValence!=4)&&
                (PInf1[i].NumCorners==(int)PInf1[i].ExpectedValence))
                MatchSing1++;
        }
        if (CClarkability>0)
            Sing1+=AddedSingularities(PatchFaces1[i].size(),PInf1[i].CurvedL,CClarkability*avgEdge);
        //if ((CClarkability>0)&&(PInf1[i].CClarkability>CClarkability))NonOKCC1++;
    }
    if (NonOKGenus1!=NonOKGenus0)
        return (NonOKGenus1<NonOKGenus0);

    if (NonOKEmitters1!=NonOKEmitters0)
        return (NonOKEmitters1<NonOKEmitters0);

    if (NonOKSize1!=NonOKSize0)
        return (NonOKSize1<NonOKSize0);

    if (MaxInternalSing1!=MaxInternalSing0)
        return (MaxInternalSing1<MaxInternalSing0);

    if (MatchSing1!=MatchSing0)
        return(MatchSing1>MatchSing0);
//    ScalarType NormNonOkVal0=((ScalarType)NonOkVal0)/((ScalarType)PInf0.size());
//    ScalarType NormNonOkVal1=((ScalarType)NonOkVal1)/((ScalarType)PInf1.size());
//    if (NormNonOkVal0!=NormNonOkVal1)return (NonOkVal1<NonOkVal0);

//    if (NonOkVal1!=NonOkVal0)
//        return (NonOkVal1<NonOkVal0);

//    //the same number of error do not remove
//    if ((NonOkVal1>0)&&(NonOkVal0>0))
//        return false;

    //if (NonOKCC1!=NonOKCC0)return (NonOKCC1<NonOKCC0);

    if ((Sing0==0)&&(Sing1==0))return true;

    return (Sing1<Sing0);
    //return true;

    //    if (Sing0!=Sing1)return (Sing1<Sing0);
    //    return true;

    //    size_t Num0=NonOkPartitions(mesh,PInf0,MinSides,MaxSides,CatmullClarkability,QThresold);
    //    size_t Num1=NonOkPartitions(mesh,PInf1,MinSides,MaxSides,CatmullClarkability,QThresold);
    //    return (Num1<=Num0);
    //    typename MeshType::ScalarType Area0=NonOkArea(mesh,PatchFaces0,PInf0,MinSides,MaxSides,QThresold);
    //    typename MeshType::ScalarType Area1=NonOkArea(mesh,PatchFaces1,PInf1,MinSides,MaxSides,QThresold);
    //    std::cout<<"Area 0:"<<Area0;
    //    std::cout<<"Area 1:"<<Area1;
    //    return (Area1<=Area0);
}

template <class MeshType>
void MeshTrace(const VertexFieldGraph<MeshType> &VFGraph,
               const CandidateTrace &CurrTrace,
               typename MeshType::ScalarType Width,
               MeshType &outMesh)
{
   typedef typename MeshType::CoordType CoordType;
   typedef typename MeshType::ScalarType ScalarType;
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
       vcg::tri::OrientedCylinder(TempMesh,P0,P1,Width,true);
       vcg::tri::Append<MeshType,MeshType>::Mesh(outMesh,TempMesh);
   }
}

template <class MeshType>
void MeshTraces(const VertexFieldGraph<MeshType> &VFGraph,
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

#endif
