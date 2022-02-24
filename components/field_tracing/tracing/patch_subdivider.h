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

#ifndef PATCH_SUBDIVIDER
#define PATCH_SUBDIVIDER

#include <vcg/space/triangle2.h>
#include <vcg/space/intersection2.h>
#include <vcg/complex/allocate.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/space/color4.h>
#include <vcg/complex/algorithms/clean.h>
#include <map>

#define EPS 0.000001

template<class MeshType>
class SelMidPointFunctor : public std::unary_function<vcg::face::Pos<typename MeshType::FaceType> ,
        typename MeshType::CoordType>
{
    typedef MeshType TriMeshType;
    typedef typename TriMeshType::CoordType    CoordType;
    typedef typename TriMeshType::VertexType   VertexType;
    typedef typename TriMeshType::FaceType     FaceType;
    typedef typename TriMeshType::ScalarType   ScalarType;

public:
    //TriMeshType *m;
    std::map<std::pair<CoordType,CoordType>, CoordType> *SplitMap;

    void operator()(VertexType &nv,
                    const vcg::face::Pos<FaceType> &ep)
    {
        VertexType *v0=ep.V();
        const VertexType *v1=ep.VFlip();

        CoordType Pos0=v0->P();
        CoordType Pos1=v1->P();

        vcg::Point2<ScalarType> t0=v0->T().P();
        vcg::Point2<ScalarType> t1=v1->T().P();
        std::pair<CoordType,CoordType> Key(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
        assert(SplitMap->count(Key)>0);
        nv.P()=(*SplitMap)[Key];

        ScalarType alpha=(Pos0-nv.P()).Norm()/(Pos0-Pos1).Norm();
        nv.T().P()=t0*(1-alpha)+t1*alpha;
        nv.SetV();
    }

    template<class FL_TYPE>
    vcg::TexCoord2<FL_TYPE,1> WedgeInterp(vcg::TexCoord2<FL_TYPE,1> &t0,
                                          vcg::TexCoord2<FL_TYPE,1> &t1)
    {
        vcg::TexCoord2<FL_TYPE,1> tmp;
        assert(t0.n()== t1.n());
        tmp.n()=t0.n();
        tmp.t()=(t0.t()+t1.t())/2;

        return tmp;
    }

    //SelMidPointFunctor(TriMeshType *_m){m=_m;}
};

template <class MeshType>
class SelEdgePred
{
    typedef MeshType TriMeshType;
    typedef typename TriMeshType::CoordType    CoordType;
    typedef typename TriMeshType::VertexType   VertexType;
    typedef typename TriMeshType::FaceType     FaceType;
    typedef typename TriMeshType::ScalarType   ScalarType;
public:

    std::map<std::pair<CoordType,CoordType>, CoordType> *SplitMap;

    bool operator()(vcg::face::Pos<typename MeshType::FaceType> ep)
    {
        CoordType Pos0=ep.f->P0(ep.z);
        CoordType Pos1=ep.f->P1(ep.z);
        std::pair<CoordType,CoordType> Key(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
        return (SplitMap->count(Key)>0);
    }
};


template <class FaceType>
bool IsInsideUV(const FaceType &face,
                const vcg::Point2<typename FaceType::ScalarType> &testUV)
{
    typedef typename FaceType::ScalarType ScalarType;
    vcg::Triangle2<ScalarType> UVTriangle(face.cV(0)->cT().P(),
                                          face.cV(1)->cT().P(),
                                          face.cV(2)->cT().P());
    return(vcg::IsInsideTrianglePoint(UVTriangle,testUV));
}


template <class MeshType>
size_t getClosestUVvert(MeshType &mesh,const vcg::Point2<typename MeshType::ScalarType> &TestPos)
{
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename vcg::Point2<ScalarType> Point2X;
    ScalarType MinD=std::numeric_limits<ScalarType>::max();
    size_t minI=0;
    for (size_t i=0;i<mesh.vert.size();i++)
    {
        Point2X UvVert=mesh.vert[i].T().P();
        ScalarType Dist=(UvVert-TestPos).Norm();
        if (Dist>MinD)continue;
        MinD=Dist;
        minI=i;
    }
    return minI;
}


template <class MeshType>
class PatchSplitter
{
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename vcg::Point2<ScalarType> Point2x;
    typedef typename vcg::Triangle2<ScalarType> Triangle2x;

    std::map<std::pair<CoordType,CoordType>, CoordType> SplitMap;
    MeshType &tri_mesh;

    void SplitCenter(MeshType &mesh,const Point2x &CenterUV)
    {
        //refine the center
        int Index=-1;
        for (size_t i=0;i<mesh.face.size();i++)
            if (IsInsideUV(mesh.face[i],CenterUV))
            {Index=i;break;}

        if (Index<0)
        {
            PatchManager<MeshType>::SetUVtoPos(mesh);
            vcg::tri::Allocator<MeshType>::AddVertex(mesh,CoordType(CenterUV.X(),CenterUV.Y(),0));
            vcg::tri::io::ExporterPLY<MeshType>::Save(mesh,"test0.ply",vcg::tri::io::Mask::IOM_FACECOLOR);
            assert(0);
        }
        assert(Index>=0);

        VertexType *v0=mesh.face[Index].V(0);
        VertexType *v1=mesh.face[Index].V(1);
        VertexType *v2=mesh.face[Index].V(2);
        size_t IndexV0=vcg::tri::Index(mesh,v0);
        size_t IndexV1=vcg::tri::Index(mesh,v1);
        size_t IndexV2=vcg::tri::Index(mesh,v2);

        vcg::Triangle2<ScalarType> UVTriangle(v0->T().P(),
                                              v1->T().P(),
                                              v2->T().P());

        CoordType bary;
        UVTriangle.InterpolationParameters(CenterUV,bary.X(),
                                           bary.Y(),bary.Z());

        //        ScalarType Eps=0.00001;

        //        //move slightly to the center if needed
        //        for (size_t i=0;i<3;i++)
        //            bary.V(i)=std::max(bary.V(i),Eps);
        //        for (size_t i=0;i<3;i++)
        //            bary.V(i)=std::min(bary.V(i),1-Eps);

        CoordType Pos3D=(v0->cP()*bary.V(0)+
                         v1->cP()*bary.V(1)+
                         v2->cP()*bary.V(2));

        vcg::Point2<typename MeshType::ScalarType> PosUV=v0->cT().P()*bary.V(0)+
                v1->cT().P()*bary.V(1)+
                v2->cT().P()*bary.V(2);

        vcg::tri::Allocator<MeshType>::AddVertex(mesh,Pos3D);
        VertexType *vnew=&mesh.vert.back();
        vnew->T().P()=PosUV;

        v0=&mesh.vert[IndexV0];
        v1=&mesh.vert[IndexV1];
        v2=&mesh.vert[IndexV2];

        mesh.face[Index].V(0)=v0;
        mesh.face[Index].V(1)=v1;
        mesh.face[Index].V(2)=vnew;

        vcg::tri::Allocator<MeshType>::AddFace(mesh,v1,v2,vnew);
        mesh.face.back().PD1()=mesh.face[Index].PD1();
        mesh.face.back().PD2()=mesh.face[Index].PD2();
        mesh.face.back().C()=mesh.face[Index].C();
        vcg::tri::Allocator<MeshType>::AddFace(mesh,v2,v0,vnew);
        mesh.face.back().PD2()=mesh.face[Index].PD2();
        mesh.face.back().PD1()=mesh.face[Index].PD1();
        mesh.face.back().PD2()=mesh.face[Index].PD2();
        mesh.face.back().C()=mesh.face[Index].C();
        mesh.face.back().C()=mesh.face[Index].C();
    }

    void UpdateSplitMap(MeshType &mesh,
                        vcg::Segment2<ScalarType> &SplitSeg,
                        std::set<std::pair<CoordType,CoordType> > &SharpSet)
    {
        Point2x IntersPos;
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                Point2x UV0= mesh.face[i].V0(j)->T().P();
                Point2x UV1= mesh.face[i].V1(j)->T().P();
                CoordType Pos0=mesh.face[i].P0(j);
                CoordType Pos1=mesh.face[i].P1(j);

                vcg::Segment2<ScalarType> Seg2(UV0,UV1);

                if (!vcg::SegmentSegmentIntersection(Seg2,SplitSeg,IntersPos))continue;

                ScalarType alpha=(UV0-IntersPos).Norm()/(UV0-UV1).Norm();
                if ((alpha<EPS)||(alpha>(1-EPS)))continue;
                CoordType Pos3D= Pos1*alpha + Pos0*(1-alpha);

                std::pair<CoordType,CoordType> key(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
                //std::cout<<"Added! "<<std::endl;
                SplitMap[key]=Pos3D;
                if (SharpSet.count(key)>0)
                {
                    std::pair<CoordType,CoordType> new_e0(std::min(Pos0,Pos3D),std::max(Pos0,Pos3D));
                    std::pair<CoordType,CoordType> new_e1(std::min(Pos1,Pos3D),std::max(Pos1,Pos3D));
                    //SharpSet.erase(key);
                    SharpSet.insert(new_e0);
                    SharpSet.insert(new_e1);
                }
            }
    }


    void GetSideSplitUV(const MeshType &mesh,
                        const std::vector<size_t> & CornersIDX,
                        const std::vector<ScalarType> & SideSplit,
                        std::vector<Point2x> &SideUV)
    {
        for (size_t i=0;i<CornersIDX.size();i++)
        {
            size_t sizeCorner=CornersIDX.size();
            size_t IndexC0=CornersIDX[i];
            size_t IndexC1=CornersIDX[(i+1)%sizeCorner];

            assert(IndexC0<mesh.vert.size());
            assert(IndexC1<mesh.vert.size());

            Point2x UV0a=mesh.vert[IndexC0].cT().P();
            Point2x UV0b=mesh.vert[IndexC1].cT().P();

            Point2x sideSplitUV=UV0a*SideSplit[i]+UV0b*(1-SideSplit[i]);
            SideUV.push_back(sideSplitUV);
        }
    }


    bool SubdivideParametrizedSubMesh(MeshType &mesh)
    {
        SelEdgePred<MeshType> SelEP;
        SelMidPointFunctor<MeshType> MidEP;//(&mesh);
        SelEP.SplitMap=&SplitMap;
        MidEP.SplitMap=&SplitMap;

        bool refined=vcg::tri::RefineE<MeshType,SelMidPointFunctor<MeshType>,SelEdgePred<MeshType> >(mesh,MidEP,SelEP);

        return refined;
    }


    void GetNewPartitionIndexes(const std::vector<bool> &MustSplit,
                                const std::vector<std::vector<size_t> > &CornersIDX,
                                std::vector<std::vector<size_t> > &NewIndex)
    {
        NewIndex.clear();
        size_t currI=0;
        for (size_t i=0;i<CornersIDX.size();i++)
        {
            if (!MustSplit[i])
            {
                NewIndex.push_back(std::vector<size_t>(1,currI));
                currI++;
            }
            else
            {
                NewIndex.resize(NewIndex.size()+1);
                for (size_t j=0;j<CornersIDX[i].size();j++)
                {
                    NewIndex.back().push_back(currI);
                    currI++;
                }
            }
        }
    }

    void SetSubPartitionIndexesOnQ(MeshType &ParamPatch,
                                   const std::vector<std::vector<Point2x> > &PartitionExtremes,
                                   const std::vector<size_t>  &NewIndex)
    {
        vcg::tri::UpdateQuality<MeshType>::FaceConstant(ParamPatch,-1);
        std::vector<Triangle2x> PartitionTris;
        for (size_t i=0;i<PartitionExtremes.size();i++)
        {
            assert(PartitionExtremes[i].size()==4);
            Point2x UV0=PartitionExtremes[i][0];
            Point2x UV1=PartitionExtremes[i][1];
            Point2x UV2=PartitionExtremes[i][2];
            Point2x UV3=PartitionExtremes[i][3];
            Triangle2x T0(UV0,UV1,UV3);
            Triangle2x T1(UV1,UV2,UV3);
            PartitionTris.push_back(T0);
            PartitionTris.push_back(T1);
        }
        assert(NewIndex.size()==PartitionExtremes.size());
        for (size_t i=0;i<ParamPatch.face.size();i++)
        {
            //find in which partition the face fall into
            Point2x FaceUV=(ParamPatch.face[i].V(0)->T().P()+
                            ParamPatch.face[i].V(1)->T().P()+
                            ParamPatch.face[i].V(2)->T().P())/3;

            int PatchIdx=-1;
            for (size_t j=0;j<PartitionTris.size();j++)
            {
                if (!vcg::IsInsideTrianglePoint(PartitionTris[j],FaceUV))continue;
                PatchIdx=j/2;
                break;
            }
            assert(PatchIdx>=0);
            assert(PatchIdx<NewIndex.size());
            ParamPatch.face[i].Q()=NewIndex[PatchIdx];
        }
    }

    void SetPartitionIndexesOnQ(std::vector<MeshType*> &ParamPatches,
                                const std::vector<bool> &MustSplit,
                                const std::vector<std::vector<std::vector<Point2x> > > &PartitionExtremes,
                                const std::vector<std::vector<size_t> > &NewIndex)
    {
        assert(ParamPatches.size()==PartitionExtremes.size());
        assert(ParamPatches.size()==MustSplit.size());
        assert(ParamPatches.size()==NewIndex.size());
        for (size_t i=0;i<ParamPatches.size();i++)
        {
            MeshType *patchMesh=ParamPatches[i];
            if (!MustSplit[i])
            {
                assert(NewIndex[i].size()==1);
                vcg::tri::UpdateQuality<MeshType>::FaceConstant((*patchMesh),NewIndex[i][0]);
            }
            else
            {
                SetSubPartitionIndexesOnQ((*patchMesh),PartitionExtremes[i],NewIndex[i]);
            }

        }
    }

    void GetPartitionCornerPos(const std::vector<MeshType*> &ParamPatches,
                               const std::vector<bool> &MustSplit,
                               const std::vector<std::vector<size_t> > &CornersIDX,
                               const std::vector<std::vector<std::vector<Point2x> > > &PartitionExtremes,
                               std::vector<std::vector<CoordType> > &CornerPos)
    {
        CornerPos.clear();
        //CornerPos.resize(NewIndex.back().back()+1);
        for (size_t i=0;i<ParamPatches.size();i++)
        {
            MeshType* currM=ParamPatches[i];
            if (!MustSplit[i])
            {
                //assert(NewIndex[i].size()==1);
                //size_t IndexP=NewIndex[i][0];
                CornerPos.resize(CornerPos.size()+1);
                //assert(PartitionExtremes[i].size()==0);
                for (size_t j=0;j<CornersIDX[i].size();j++)
                {
                    CoordType Pos=(*currM).vert[CornersIDX[i][j]].cP();
                    CornerPos.back().push_back(Pos);
                    //CornerPos[IndexP].push_back(Pos);
                }
            }else
            {
                for (size_t j=0;j<PartitionExtremes[i].size();j++)
                {
                    CornerPos.resize(CornerPos.size()+1);
                    //size_t IndexP=NewIndex[i][j];
                    //std::cout<<"IndexP:"<<IndexP<<std::endl;
                    for (size_t k=0;k<PartitionExtremes[i][j].size();k++)
                    {
                        Point2x UV=PartitionExtremes[i][j][k];
                        size_t IndexC=getClosestUVvert(*currM,UV);
                        CoordType Pos=(*currM).vert[IndexC].cP();
                        //CornerPos[IndexP].push_back(Pos);
                        CornerPos.back().push_back(Pos);
                    }
                }
            }
        }
    }

    void GetCornerIndexFromPos(const MeshType &mesh,
                               std::vector<std::vector<CoordType> > &CornerPos,
                               std::vector<std::vector<size_t> > &CornerIDx)
    {
        std::map<CoordType,size_t> CoordToIDX;
        for (size_t i=0;i<mesh.vert.size();i++)
            CoordToIDX[mesh.vert[i].P()]=i;

        CornerIDx.resize(CornerPos.size());
        for (size_t i=0;i<CornerPos.size();i++)
        {
            for (size_t j=0;j<CornerPos[i].size();j++)
            {
                CoordType pos=CornerPos[i][j];
                assert(CoordToIDX.count(pos)>0);
                size_t IndexV=CoordToIDX[pos];
                CornerIDx[i].push_back(IndexV);
                //std::cout<<IndexV<<",";
            }
            //std::cout<<";"<<std::endl;
        }
    }

    void GetSideSplitUV(const std::vector<MeshType*> ParamPatches,
                        const std::vector<std::vector<ScalarType> > &SideSplit,
                        const std::vector<std::vector<size_t> > &PatchCorners,
                        std::vector<std::vector<Point2x> > &SideUV)
    {
        SideUV.clear();
        SideUV.resize(ParamPatches.size());
        for (size_t i=0;i<ParamPatches.size();i++)
            GetSideSplitUV((*ParamPatches[i]),PatchCorners[i],SideSplit[i],SideUV[i]);
    }

    void GetCenterSubdivision(const std::vector<std::vector<Point2x> > &SideUV,
                              std::vector<Point2x > &CenterUV)
    {
        CenterUV.clear();
        for (size_t i=0;i<SideUV.size();i++)
        {
            Point2x AvgUv(0,0);
            for (size_t j=0;j<SideUV[i].size();j++)
                AvgUv+=SideUV[i][j];
            AvgUv/=SideUV[i].size();
            CenterUV.push_back(AvgUv);
        }
    }

    void GetSplittingPartition(const MeshType &mesh,
                               const std::vector<size_t> &CornersIDX,
                               const std::vector<Point2x> &SideSplit,
                               const Point2x &CenterUV,
                               std::vector<vcg::Segment2<ScalarType> > &SplitSegs,
                               std::vector<std::vector<Point2x> > &PartitionExtremes)
    {
        std::vector<Point2x> NewCorners;
        std::vector<Point2x> OldCorners;
        for (size_t i=0;i<CornersIDX.size();i++)
        {
            //            size_t sizeCorner=CornersIDX.size();
            size_t IndexC0=CornersIDX[i];
            //            size_t IndexC1=CornersIDX[(i+1)%sizeCorner];

            Point2x UV0a=mesh.vert[IndexC0].cT().P();
            //            Point2x UV0b=mesh.vert[IndexC1].T().P();

            Point2x UV0=SideSplit[i];
            Point2x UV1=CenterUV;

            //move a bit further
            Point2x UV0_mov=UV1+(UV0-UV1)*1.00001;

            //move a bit toward outside
            Point2x UV1_mov=UV0+(UV1-UV0)*0.99999;

            vcg::Segment2<ScalarType> SplitSeg(UV0_mov,UV1_mov);
            SplitSegs.push_back(SplitSeg);

            NewCorners.push_back(UV0);
            OldCorners.push_back(UV0a);
        }

        PartitionExtremes.clear();
        PartitionExtremes.resize(CornersIDX.size());
        size_t sizeCorners=PartitionExtremes.size();
        for (size_t i=0;i<PartitionExtremes.size();i++)
        {
            PartitionExtremes[i].clear();
            PartitionExtremes[i].push_back(NewCorners[i]);
            PartitionExtremes[i].push_back(OldCorners[(i+1)%sizeCorners]);
            PartitionExtremes[i].push_back(NewCorners[(i+1)%sizeCorners]);
            PartitionExtremes[i].push_back(CenterUV);
        }
    }

    void GetSplittingPartitions(const std::vector<MeshType*> ParamPatches,
                                const std::vector<std::vector<size_t> > &CornersIDX,
                                const std::vector<std::vector<Point2x> > &SideSplit,
                                const std::vector<Point2x> &CenterUV,
                                std::vector<std::vector<vcg::Segment2<ScalarType> > > &SplitSegs,
                                std::vector<std::vector<std::vector<Point2x> > > &PartitionExtremes)
    {
        PartitionExtremes.clear();
        PartitionExtremes.resize(ParamPatches.size());
        SplitSegs.clear();
        SplitSegs.resize(ParamPatches.size());
        for (size_t i=0;i<ParamPatches.size();i++)
            GetSplittingPartition((*ParamPatches[i]),CornersIDX[i],SideSplit[i],
                                  CenterUV[i],SplitSegs[i],PartitionExtremes[i]);
    }

    void CheckIndex(const std::vector<std::vector<size_t> > &NewCorners,
                    const size_t MaxPatchIndex)
    {
        std::vector<std::set<size_t> > TestIDX;
        std::vector<std::set<CoordType> > TestPos;
        TestIDX.resize(MaxPatchIndex+1);
        TestPos.resize(MaxPatchIndex+1);
        for (size_t k=0;k<tri_mesh.face.size();k++)
        {
            size_t IndexQ=tri_mesh.face[k].Q();
            assert(IndexQ<TestIDX.size());
            TestIDX[IndexQ].insert(vcg::tri::Index(tri_mesh,tri_mesh.face[k].V(0)));
            TestIDX[IndexQ].insert(vcg::tri::Index(tri_mesh,tri_mesh.face[k].V(1)));
            TestIDX[IndexQ].insert(vcg::tri::Index(tri_mesh,tri_mesh.face[k].V(2)));
            TestPos[IndexQ].insert(tri_mesh.face[k].P(0));
            TestPos[IndexQ].insert(tri_mesh.face[k].P(1));
            TestPos[IndexQ].insert(tri_mesh.face[k].P(2));
        }
        for (size_t i=0;i<NewCorners.size();i++)
        {
            //std::cout<<"There are "<<NewCorners[i].size()<<" corners"<<std::endl;
            bool found0=false;
            bool found1=false;
            for (size_t j=0;j<TestIDX.size();j++)
            {
                bool isOK=true;
                for (size_t k=0;k<NewCorners[i].size();k++)
                {
                    size_t Idx=NewCorners[i][k];
                    isOK&=(TestIDX[j].count(Idx)>0);
                }
                if (isOK)
                {
                    if (i==j)
                    {
                        found0=true;break;
                    }else
                    {

                        std::cout<<"Found Partition "<<j<<" instead of "<<i<<std::endl;
                        break;
                    }
                }
            }
            for (size_t j=0;j<TestPos.size();j++)
            {
                bool isOK=true;
                for (size_t k=0;k<NewCorners[i].size();k++)
                {
                    CoordType Pos=tri_mesh.vert[NewCorners[i][k]].P();
                    isOK&=(TestPos[j].count(Pos)>0);
                }
                if (isOK)
                {
                    if (i==j)
                    {
                        found1=true;break;
                    }else
                    {

                        std::cout<<"Found Partition "<<j<<" instead of "<<i<<std::endl;
                        break;
                    }
                }
            }
            if (!found0)
            {
                std::cout<<" 0 WARNING NOT FOUND CORNER SET "<<i<<std::endl;
                //                MeshType test;
                //                for (size_t k=0;k<CornerPos[i].size();k++)
                //                    vcg::tri::Allocator<MeshType>::AddVertex(test,CornerPos[i][k]);

                //                vcg::tri::io::ExporterPLY<MeshType>::Save(test,"WARNING.ply");
                assert(0);
            }
            if (!found1)
            {
                std::cout<<" 0 WARNING NOT FOUND CORNER SET "<<i<<std::endl;
                //                MeshType test;
                //                for (size_t k=0;k<CornerPos[i].size();k++)
                //                    vcg::tri::Allocator<MeshType>::AddVertex(test,CornerPos[i][k]);

                //                vcg::tri::io::ExporterPLY<MeshType>::Save(test,"WARNING.ply");
                assert(0);
            }
            //assert(found);
        }
    }

public:

    void Subdivide(const std::vector<std::vector<size_t> > &PatchFaces,
                   const std::vector<std::vector<size_t> > &CornersIDX,
                   const std::vector<bool> &MustSplit,
                   const std::vector<std::vector<ScalarType> > & SideSplit,
                   std::vector<std::vector<size_t> > &NewFacePaches,
                   std::vector<std::vector<size_t> > &NewCorners)
    {

        assert(MustSplit.size()==CornersIDX.size());
        assert(MustSplit.size()==SideSplit.size());
        //assert(MustSplit.size()==CenterUV.size());

        std::vector<MeshType*> ParamPatches;

        //parametrize
        std::vector<std::vector<size_t> > PatchCorners;
        PatchCorners.resize(PatchFaces.size());
        for (size_t i=0;i<PatchFaces.size();i++)
        {
            ParamPatches.push_back(new MeshType());
            PatchManager<MeshType>::ComputeParametrizedSubMesh(tri_mesh,*ParamPatches[i],
                                       PatchFaces[i],CornersIDX[i],
                                       PatchCorners[i],
                                       Arap,true,true,false);
        }

        std::vector<std::vector<Point2x> > SideUV;
        GetSideSplitUV(ParamPatches,SideSplit,PatchCorners,SideUV);

        std::vector<Point2x > CenterUV;
        GetCenterSubdivision(SideUV,CenterUV);

        std::vector<std::vector<vcg::Segment2<ScalarType> > > SplitSegs;
        std::vector<std::vector<std::vector<Point2x> > > PartitionExtremes;

        GetSplittingPartitions(ParamPatches,PatchCorners,SideUV,CenterUV,SplitSegs,PartitionExtremes);


        //split centers
        for (size_t i=0;i<ParamPatches.size();i++)
        {
            if (MustSplit[i])
            {
                SplitCenter(*ParamPatches[i],CenterUV[i]);
                ParamPatches[i]->UpdateAttributes();
            }
        }
        //then collect the split map
        bool refined=false;

        std::vector<std::pair<CoordType,CoordType> > SharpCoords;
        tri_mesh.GetSharpCoordPairs(SharpCoords);
        std::set<std::pair<CoordType,CoordType> > SharpSet(SharpCoords.begin(),SharpCoords.end());

        do {
            refined=false;
            //collect all split operations
            SplitMap.clear();

            for (size_t i=0;i<SplitSegs.size();i++)
                if (MustSplit[i]){
                    for (size_t j=0;j<SplitSegs[i].size();j++)
                        UpdateSplitMap(*ParamPatches[i],SplitSegs[i][j],SharpSet);
                }

            for (size_t i=0;i<ParamPatches.size();i++)
                    refined|=SubdivideParametrizedSubMesh(*ParamPatches[i]);
        }while (refined);

        std::vector<std::vector<size_t> > NewIndex;
        GetNewPartitionIndexes(MustSplit,PatchCorners,NewIndex);

        SetPartitionIndexesOnQ(ParamPatches,MustSplit,PartitionExtremes,NewIndex);

        std::vector<std::vector<CoordType> > CornerPos;
        GetPartitionCornerPos(ParamPatches,MustSplit,PatchCorners,PartitionExtremes,CornerPos);


        tri_mesh.Clear();
        for (size_t i=0;i<ParamPatches.size();i++)
        {
            vcg::tri::Append<MeshType,MeshType>::Mesh(tri_mesh,*ParamPatches[i]);
            delete(ParamPatches[i]);
        }

        vcg::tri::Clean<MeshType>::RemoveDuplicateVertex(tri_mesh);
        vcg::tri::Clean<MeshType>::RemoveUnreferencedVertex(tri_mesh);
        vcg::tri::Allocator<MeshType>::CompactEveryVector(tri_mesh);
        tri_mesh.UpdateAttributes();

        SharpCoords.clear();
        SharpCoords=std::vector<std::pair<CoordType,CoordType> >(SharpSet.begin(),SharpSet.end());
        tri_mesh.UpdateFromCoordPairs(SharpCoords,false);

        GetCornerIndexFromPos(tri_mesh,CornerPos,NewCorners);



        size_t MaxPatchIndex=NewIndex.back().back();
        CheckIndex(NewCorners,MaxPatchIndex);

        NewFacePaches.clear();
        NewFacePaches.resize(MaxPatchIndex+1);
        for (size_t i=0;i<tri_mesh.face.size();i++)
        {
            size_t IndexPatch=tri_mesh.face[i].Q();
            assert(IndexPatch<NewFacePaches.size());
            NewFacePaches[IndexPatch].push_back(i);
        }

        assert(NewCorners.size()==NewFacePaches.size());

        //PreProcessMesh(tri_mesh);

        //                for (size_t i=0;i<tri_mesh.face.size();i++)
        //                    tri_mesh.face[i].C()=vcg::Color4b::Scatter(MaxPatchIndex,tri_mesh.face[i].Q());

        //        NewFacePaches.clear();
        //        NewFacePaches.resize(PartitionTris.size()/2);


        //        for (size_t i=0;i<NewFacePaches.size();i++)
        //        {
        //            vcg::Color4b col=vcg::Color4b::Scatter(NewFacePaches.size(),i);
        //            for (size_t j=0;j<NewFacePaches[i].size();j++)
        //            {
        //                size_t IndexF=NewFacePaches[i][j];
        //                subdividedUV.face[IndexF].C()=col;
        //            }
        //        }
        //SetUVtoPos(subdividedUV);
        //        vcg::tri::io::ExporterPLY<MeshType>::Save(tri_mesh,"subdivided.ply",vcg::tri::io::Mask::IOM_FACECOLOR);
    }



    PatchSplitter(MeshType &_tri_mesh):tri_mesh(_tri_mesh)
    {}
};
#endif
