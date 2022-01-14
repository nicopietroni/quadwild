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

#ifndef VERT_FIELD_GRAPH
#define VERT_FIELD_GRAPH

#include <vcg/math/matrix33.h>
#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>
#include <vcg/complex/algorithms/attribute_seam.h>

#define DELTA 0.000000001

template <class MeshType>
typename MeshType::ScalarType TotArea(const MeshType &mesh)
{
    typename MeshType::ScalarType TotArea=0;
    for (size_t i=0;i<mesh.face.size();i++)
        TotArea+=vcg::DoubleArea(mesh.face[i]);
    return TotArea/2;
}

template <class FaceType>
void GetSortedPos(vcg::face::Pos<FaceType> &StartPos,
                  std::vector<vcg::face::Pos<FaceType> > &PosSeq)
{
    PosSeq.clear();
    vcg::face::Pos<FaceType> CurrPos=StartPos;
    do
    {
        PosSeq.push_back(CurrPos);
        CurrPos.FlipE();
        CurrPos.FlipF();
    }while (CurrPos!=StartPos);
}

template <class FaceType>
size_t WhichIndex(FaceType &f,typename FaceType::VertexType *v)
{
    if (f.V(0)==v)return 0;
    if (f.V(1)==v)return 1;
    if (f.V(2)==v)return 2;
    assert(0);
}

template <class MeshType>
class VertSplitter{

    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;


    static void ExtractVertex(const MeshType & srcMesh,
                              const FaceType & f,
                              int whichWedge,
                              const MeshType & dstMesh,
                              VertexType & v)
    {
        (void)srcMesh;
        (void)dstMesh;
        v.P() = f.cP(whichWedge);
        v.T() = f.cWT(whichWedge);
    }

    static bool CompareVertex(const MeshType & m,
                              const VertexType & vA,
                              const VertexType & vB)
    {
        (void)m;
        return (vA.cT() == vB.cT());
    }

public:

    static void SplitAlongEdgeSel(MeshType &mesh)
    {

        //set default WEdgeCoords
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
                mesh.face[i].WT(j).P()=vcg::Point2<ScalarType>(0,0);

        //then set discontinuities along features
        vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);
        for (size_t i=0;i<mesh.face.size();i++)
        {
            FaceType *f=&mesh.face[i];
            for (size_t j=0;j<3;j++)
            {
                if (f->IsB(j))continue;
                if (!f->IsFaceEdgeS(j))continue;
                vcg::face::Pos<FaceType> CurrPos(f,j);
                vcg::face::Pos<FaceType> StartPos=CurrPos;
                vcg::Point2<ScalarType> CurrUV(0,0);
                if (!CurrPos.V()->IsV())
                {
                    CurrPos.V()->SetV();
                    do
                    {
                        size_t VIndex=WhichIndex(*CurrPos.F(),CurrPos.V());
                        CurrPos.F()->WT(VIndex).P()=CurrUV;
                        CurrPos.FlipE();
                        CurrPos.FlipF();
                        f=CurrPos.F();
                        size_t curr_e=CurrPos.E();
                        if (f->IsFaceEdgeS(curr_e))
                            CurrUV+=vcg::Point2<ScalarType>(1,1);
                    }while (CurrPos!=StartPos);
                }
            }
        }
        vcg::tri::AttributeSeam::SplitVertex(mesh, ExtractVertex, CompareVertex);
    }
};


struct EdgeVert
{
    size_t EV0;
    size_t EV1;
    size_t CurrV;

    EdgeVert(size_t _EV0,size_t _EV1,size_t _CurrV)
    {
        EV0=std::min(_EV0,_EV1);
        EV1=std::max(_EV0,_EV1);
        CurrV=_CurrV;
    }

    inline bool operator ==(const EdgeVert &left)const
    {
        return ((EV0==left.EV0)&&
                (EV1==left.EV1)&&
                (CurrV==left.CurrV));
    }

    inline bool operator <(const EdgeVert &left)const
    {
        if ((EV0==left.EV0)&&
                (EV1==left.EV1))
            return (CurrV<left.CurrV);
        if (EV0==left.EV0)
            return (EV1<left.EV1);
        return (EV0<left.EV0);
    }
};


struct EdgeVertKeyHasher
{
  std::size_t operator()(const EdgeVert& k) const
  {
    //using std::size_t;
    //using std::hash;
    //using std::string;
    const size_t _HASH_P0 = 73856093u;
    const size_t _HASH_P1 = 19349663u;
    const size_t _HASH_P2 = 83492791u;
    return ((std::hash<size_t>()(k.EV0)*_HASH_P0)
             ^ (std::hash<size_t>()(k.EV1) *_HASH_P1)
             ^ (std::hash<size_t>()(k.CurrV) *_HASH_P2));
  }
};

template <class MeshType>
class VertexFieldGraph
{
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;

    //the input mesh to trace
    MeshType &mesh;

    //propagation steps
    size_t PropagationSteps;

    size_t TMark;

    std::vector<bool> RealBorderVert;

    bool DebugMsg;

public:

    std::map<EdgeVert,size_t> EdgeBorderDir;
    //std::unordered_map<EdgeVert,size_t,EdgeVertKeyHasher> EdgeBorderDir;

    struct NeighInfo
    {
        bool twin;
        size_t Node;
        ScalarType Angle;
        ScalarType Dist;
        bool Active;
        bool Direct;

        inline bool operator <(const NeighInfo &NInfo)const
        {
            //return (NextWeight>NInfo.NextWeight);
            return (Angle<NInfo.Angle);
        }

        NeighInfo(size_t _Node,
                  ScalarType _Angle,
                  ScalarType _Dist)
        {
            Node=_Node;
            Dist=_Dist;
            Angle=_Angle;
            Active=true;
            twin=false;
            Direct=true;
        }
    };

    static CoordType MakeTangentTO(const CoordType &Norm0,
                                   const CoordType &TestD)
    {
        CoordType TestDir=TestD;
        if(fabs(TestDir*Norm0)>0.1)
        {
            CoordType RotAxis=Norm0^TestDir;
            CoordType TargetN=TestDir^RotAxis;
            TargetN.Normalize();
            vcg::Matrix33<ScalarType> rot=vcg::RotationMatrix(TargetN,Norm0);
            TestDir=rot*TestDir;
            TestDir.Normalize();
        }
        return TestDir;
    }

    //return the direction for a vertex that is closest to an
    //input direction expressed as a CoordType
    size_t GetClosestDirTo(const size_t IndexV0,
                           const CoordType Dir)const
    {
        ScalarType Norm=-1;
        size_t BestD=4;
        for (size_t i=0;i<4;i++)
        {
            CoordType DirTest=GetDirection(IndexV0,i);
            ScalarType TestNorm=(DirTest*Dir);
            if (TestNorm<Norm)continue;
            Norm=TestNorm;
            BestD=i;
        }
        assert(BestD<4);
        return BestD;
    }

    void AddTwinsConnections(const size_t &IndexV0,const size_t &IndexV1)
    {
        size_t DirMap[4];
        DirMap[0]=FollowDirection(IndexV0,IndexV1,0);
        DirMap[1]=FollowDirection(IndexV0,IndexV1,1);
        DirMap[2]=FollowDirection(IndexV0,IndexV1,2);
        DirMap[3]=FollowDirection(IndexV0,IndexV1,3);
        for (size_t dir=0;dir<4;dir++)
        {
            size_t Node0=IndexNode(IndexV0,dir);
            size_t Node1=IndexNode(IndexV1,DirMap[dir]);
            assert(Node0!=Node1);
            assert(Node0<Nodes.size());
            assert(Node1<Nodes.size());
            //check the one that have the most outcoming directions
            if (Nodes[Node0].NonTwinNeigh()<Nodes[Node1].NonTwinNeigh())
            {
                Nodes[Node0].Neigh.push_back(NeighInfo(Node1,0,DELTA));
                Nodes[Node0].Neigh.back().twin=true;
            }
            else
            {
                Nodes[Node1].Neigh.push_back(NeighInfo(Node0,0,DELTA));
                Nodes[Node1].Neigh.back().twin=true;
            }
        }
    }

    void AddTwinsConnections(const std::vector<size_t> &VertIndexes)
    {
        for (size_t i=0;i<VertIndexes.size()-1;i++)
            for (size_t j=(i+1);j<VertIndexes.size();j++)
                AddTwinsConnections(VertIndexes[i],VertIndexes[j]);
    }

    void InitRealBorderVert()
    {
        RealBorderVert.clear();
        RealBorderVert.resize(mesh.vert.size(),false);

        MeshType swap;
        vcg::tri::Append<MeshType,MeshType>::Mesh(swap,mesh);
        vcg::tri::Clean<MeshType>::RemoveDuplicateVertex(swap);
        swap.UpdateAttributes();

        std::map<CoordType,size_t> CoordVert;
        for (size_t i=0;i<swap.vert.size();i++)
        {
            if (!swap.vert[i].IsB())continue;
            CoordType Pos=swap.vert[i].P();
            CoordVert[Pos]=i;
        }

        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (!mesh.vert[i].IsB())continue;
            CoordType Pos=mesh.vert[i].P();
            if (CoordVert.count(Pos)==0)continue;
            RealBorderVert[i]=true;
        }

    }

    void AddTwinsConnections()
    {
        std::map<CoordType,std::vector<size_t> > CoordVert;
        //        RealBorderVert.clear();
        //        RealBorderVert.resize(mesh.vert.size(),false);
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (!mesh.vert[i].IsB())continue;
            //RealBorderVert[i]=true;
            CoordType Pos=mesh.vert[i].P();
            CoordVert[Pos].push_back(i);
        }


        //set the twins ones
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            CoordType Pos=mesh.vert[i].P();
            if (CoordVert[Pos].size()<=1)continue;

            //RealBorderVert[i]=false;
            std::vector<size_t> NodesI;
            IndexNodes(i,NodesI);
            for (size_t j=0;j<NodesI.size();j++)
                Nodes[NodesI[j]].HasTwin=true;
        }

        //then add the connections
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (!mesh.vert[i].IsB())continue;
            CoordType Pos=mesh.vert[i].P();
            AddTwinsConnections(CoordVert[Pos]);
        }

        InitRealBorderVert();
    }

    bool IsBorder(const size_t &IndexN)const
    {
        size_t IndexV=NodeVertI(IndexN);
        return (mesh.vert[IndexV].IsB());
    }

    void AddConnections(const std::vector<std::vector<size_t> > &VNeigh)
    {
        for (size_t i=0;i<VNeigh.size();i++)
            for (size_t j=0;j<VNeigh[i].size();j++)
            {
                //get the two vertices
                size_t IndexV0=i;
                size_t IndexV1=VNeigh[i][j];
                //then for each direction
                for (size_t Dir0=0;Dir0<4;Dir0++)
                {
                    size_t IndexN0=IndexNode(IndexV0,Dir0);
                    ScalarType Angle,Dist;
                    size_t Dir1;
                    ConnectionData(IndexV0,Dir0,IndexV1,Angle,Dist,Dir1);
                    if (Angle>=90)continue;
                    size_t IndexN1=IndexNode(IndexV1,Dir1);
                    Nodes[IndexN0].Neigh.push_back(NeighInfo(IndexN1,Angle,Dist));
                }
            }
    }

    void RemoveDeadEnd()
    {
        bool Removed=false;
        size_t step=0;
        do
        {
            Removed=false;
            step++;
            std::vector<bool> DeadEnd(Nodes.size(),false);
            for (size_t i=0;i<Nodes.size();i++)
                DeadEnd[i]=(Nodes[i].Neigh.size()==0);

            for (size_t i=0;i<Nodes.size();i++)
            {
                if (DeadEnd[i])continue;
                bool remove_session=false;
                std::vector<bool> MustRemoved(Nodes[i].Neigh.size(),false);
                for (size_t j=0;j<Nodes[i].Neigh.size();j++)
                {
                    size_t NextV=Nodes[i].Neigh[j].Node;
                    if (IsBorder(NextV))continue;
                    if (!DeadEnd[NextV])continue;
                    MustRemoved[j]=true;
                    remove_session=true;
                }
                if (remove_session)
                {
                    std::vector<NeighInfo> Swap;
                    for (size_t j=0;j<MustRemoved.size();j++)
                    {
                        if (MustRemoved[j])continue;
                        Swap.push_back(Nodes[i].Neigh[j]);
                    }
                    Nodes[i].Neigh=Swap;
                    Removed=true;
                }
            }
        }while (Removed);

        if (DebugMsg)
            std::cout<<"permormed removed steps "<<step<<std::endl;
    }

    void GetEdgeDir(const size_t &IndexV0,
                           const size_t &IndexV1,
                           size_t &DirN0,
                           size_t &DirN1)const
    {
        CoordType Dir=Mesh().vert[IndexV1].P()-
                      Mesh().vert[IndexV0].P();
        Dir.Normalize();
        DirN0=GetClosestDirTo(IndexV0,Dir);
        DirN1=GetClosestDirTo(IndexV1,-Dir);
    }

    void GetEdgeNodes(const size_t &IndexV0,const size_t &IndexV1,
                      size_t &IndexN0,size_t &IndexN1)const
    {
        CoordType Dir=Mesh().vert[IndexV1].P()-
                      Mesh().vert[IndexV0].P();
        Dir.Normalize();
        size_t DirN0=GetClosestDirTo(IndexV0,Dir);
        size_t DirN1=GetClosestDirTo(IndexV1,Dir);
        IndexN0=IndexNode(IndexV0,DirN0);
        IndexN1=IndexNode(IndexV1,DirN1);
    }

private:

    //the node structure
    struct Node
    {
        //        size_t IndexV;
        //        size_t IndexDir;
        std::vector<NeighInfo> Neigh;
        bool Active;
        //bool IsBorder;
        size_t TMark;
        bool Selected;
        bool HasTwin;

        void AddNeigh(size_t _NextNode,
                      ScalarType _NextWeight)
        {
            Neigh.push_back(NeighInfo(_NextNode,_NextWeight));
        }

        void SortNeigh()
        {
            std::sort(Neigh.begin(),Neigh.end());
        }

        Node()
        {
            //            IndexV=0;
            //            IndexDir=0;
            Active=true;
            TMark=0;
            //IsBorder=false;
            Selected=false;
            HasTwin=false;
        }

        size_t NonTwinNeigh()
        {
            size_t num=0;
            for (size_t i=0;i<Neigh.size();i++)
                if (!Neigh[i].twin)num++;
            return num;
        }
    };


    void GetVertStar(std::vector<std::vector<size_t> > &VNeigh)
    {
        VNeigh.clear();
        VNeigh.resize(mesh.vert.size());
        for (size_t i=0;i<mesh.face.size();i++)
            for (int j=0;j<mesh.face[i].VN();j++)
            {
                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
                //add connection only with border singularities not the rest

                bool InternalSing0=((IsSingVert[IndexV0])&&(!mesh.vert[IndexV0].IsB()));
                bool InternalSing1=((IsSingVert[IndexV1])&&(!mesh.vert[IndexV1].IsB()));
                if (!remove_sign_connections)
                {
                    InternalSing0=false;
                    InternalSing1=false;
                }
                if ((InternalSing0)||(InternalSing1))continue;//no add connection with singularities
                VNeigh[IndexV0].push_back(IndexV1);
                VNeigh[IndexV1].push_back(IndexV0);
            }
    }

    void RemoveOnBorderLinks(std::vector<std::vector<size_t> > &VNeigh)
    {
        //first save all the borders one
        std::set<std::pair<size_t,size_t> > BorderEdges;
        for (size_t i=0;i<mesh.face.size();i++)
            for (int j=0;j<mesh.face[i].VN();j++)
            {
                if (!vcg::face::IsBorder(mesh.face[i],j))continue;
                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
                BorderEdges.insert(std::pair<size_t,size_t>(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1)));
            }

        for (size_t i=0;i<VNeigh.size();i++)
        {
            size_t IndexV0=i;
            std::vector<size_t> NewNeigh;
            for (size_t j=0;j<VNeigh[i].size();j++)
            {
                size_t IndexV1=VNeigh[i][j];
                std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
                if (BorderEdges.count(key)>0)continue;
                NewNeigh.push_back(IndexV1);
            }
            VNeigh[i]=NewNeigh;
        }
    }

    void PropagateNeigh(std::vector<std::vector<size_t> > &VNeigh)
    {
        std::vector<std::vector<size_t> > newNeigh;

        for (size_t i=0;i<VNeigh.size();i++)
        {
            newNeigh.push_back(VNeigh[i]);
            for (size_t j=0;j<VNeigh[i].size();j++)
            {
                size_t IndexN=VNeigh[i][j];
                if (mesh.vert[IndexN].IsB()) continue;
                for (size_t k=0;k<VNeigh[IndexN].size();k++)
                {
                    if (mesh.vert[VNeigh[IndexN][k]].IsB())continue;
                    if (VNeigh[IndexN][k]==i)continue;
                    newNeigh.back().push_back(VNeigh[IndexN][k]);
                }
            }

            std::sort(newNeigh.back().begin(),newNeigh.back().end());
            std::vector<size_t>::iterator it;
            it = std::unique (newNeigh.back().begin(),newNeigh.back().end());
            newNeigh.back().resize( std::distance(newNeigh.back().begin(),it) );
        }
        VNeigh=newNeigh;
    }



    void ConnectionData(const size_t &IndexV0,
                        const size_t &DirectionV0,
                        const size_t &IndexV1,
                        ScalarType &Angle,
                        ScalarType &Dist,
                        size_t &Direction_next)
    {
        CoordType TargetD=GetDirection(IndexV0,DirectionV0);

        CoordType Pos0=mesh.vert[IndexV0].P();
        CoordType Pos1=mesh.vert[IndexV1].P();

        CoordType TestDir=(Pos1-Pos0);
        TestDir.Normalize();
        CoordType Norm0=mesh.vert[IndexV0].N();
        Norm0.Normalize();

        //check if orthogonal
        if(fabs(TestDir*Norm0)>0.1)
        {
            CoordType RotAxis=Norm0^TestDir;
            CoordType TargetN=TestDir^RotAxis;
            TargetN.Normalize();
            vcg::Matrix33<ScalarType> rot=vcg::RotationMatrix(TargetN,Norm0);
            TestDir=rot*TestDir;
            TestDir.Normalize();
        }
        Angle=fabs(vcg::Angle(TargetD,TestDir));
        Angle=Angle* 180/M_PI;
        Dist=(Pos1-Pos0).Norm();
        Direction_next=FollowDirection(IndexV0,IndexV1,DirectionV0);
    }

public:

    //follow the direction from one vertex to another
    //it return the index of the direction in IndexV1
    size_t FollowDirection(const size_t IndexV0,
                           const size_t IndexV1,
                           const size_t DirectionV0)
    {
        CoordType N0=mesh.vert[IndexV0].cN();
        CoordType N1=mesh.vert[IndexV1].cN();

        ///find the rotation matrix that maps between normals
        vcg::Matrix33<ScalarType> rotation=vcg::RotationMatrix(N0,N1);
        CoordType DirV0=GetDirection(IndexV0,DirectionV0);
        DirV0=rotation*DirV0;
        size_t BestD=GetClosestDirTo(IndexV1,DirV0);
        return BestD;
    }

private:

    //return a direction given an index 0<=Direction<=3 for a given vertex
    CoordType GetDirection(const size_t IndexV,
                           const size_t Direction)const
    {
        assert(Direction>=0);
        assert(Direction<4);
        if (Direction==0)
            return mesh.vert[IndexV].PD1();
        if (Direction==1)
            return mesh.vert[IndexV].PD2();
        if (Direction==2)
            return (-mesh.vert[IndexV].PD1());
        //if (Direction==3)
        assert(Direction==3);
        return (-mesh.vert[IndexV].PD2());
    }

    void CheckTangentField()
    {
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            CoordType Norm=mesh.vert[i].N();
            Norm.Normalize();
            for (size_t j=0;j<4;j++)
            {
                CoordType Dir= GetDirection(i,j);
                Dir.Normalize();
                if (fabs(Dir*Norm)>0.01)
                {
                    std::cout<<"WARNING NON TANGENT FIELD"<<std::endl;
                    std::cout<<"DOT:"<<fabs(Dir*Norm)<<std::endl;
                    std::cout<<"DIR:"<<Dir.X()<<","<<Dir.Y()<<","<<Dir.Z()<<std::endl;
                    std::cout<<"NORM:"<<Norm.X()<<","<<Norm.Y()<<","<<Norm.Z()<<std::endl;
                }
            }
        }
    }

    std::vector<Node> Nodes;

    void InitDirectConnections()
    {
        std::set<std::pair<CoordType,CoordType> > EdgeSet;
        for (size_t i=0;i<Mesh().face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                CoordType P0=Mesh().face[i].P0(j);
                CoordType P1=Mesh().face[i].P1(j);
                std::pair<CoordType,CoordType> key(std::min(P0,P1),std::max(P0,P1));
                EdgeSet.insert(key);
            }
        }
        for (size_t i=0;i<NumNodes();i++)
        {
            size_t IndexN0=i;
            CoordType P0=NodePos(IndexN0);
            for (size_t j=0;j<NumNeigh(i);j++)
            {
                if (TwinNeigh(IndexN0,j))
                {
                    Nodes[IndexN0].Neigh[j].Direct=true;
                    continue;
                }
                size_t IndexN1=NodeNeigh(i,j);
                CoordType P1=NodePos(IndexN1);

                Nodes[IndexN0].Neigh[j].Direct=true;
                std::pair<CoordType,CoordType> key(std::min(P0,P1),std::max(P0,P1));
                if (EdgeSet.count(key)==0)
                    Nodes[IndexN0].Neigh[j].Direct=false;
            }
        }
    }

    std::map<std::pair<size_t,size_t>,vcg::face::Pos<FaceType> > VertPos;
    //std::unordered_map<std::pair<size_t,size_t>,vcg::face::Pos<FaceType> > VertPos;

    void InitVertPosMap()
    {
        VertPos.clear();
        for (size_t i=0;i<Mesh().face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                size_t IndexV0=vcg::tri::Index(Mesh(),Mesh().face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(Mesh(),Mesh().face[i].V1(j));
                std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),
                                             std::max(IndexV0,IndexV1));
                vcg::face::Pos<FaceType> currPos(&Mesh().face[i],j);
                VertPos[key]=currPos;
            }
    }



    void InitBorderDirMap()
    {
        //do the same for borders
        for (size_t i=0;i<Mesh().face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!vcg::face::IsBorder(Mesh().face[i],j))continue;
                size_t IndexV0=vcg::tri::Index(Mesh(),Mesh().face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(Mesh(),Mesh().face[i].V1(j));
                size_t DirFlatV0,DirFlatV1;
                GetEdgeDir(IndexV0,IndexV1,DirFlatV0,DirFlatV1);

                assert(IndexV0!=IndexV1);
                size_t MinV=std::min(IndexV0,IndexV1);
                size_t MaxV=std::max(IndexV0,IndexV1);

                EdgeVert EdgeKey0(MinV,MaxV,IndexV0);
                EdgeVert EdgeKey1(MinV,MaxV,IndexV1);

                assert(EdgeBorderDir.count(EdgeKey0)==0);
                EdgeBorderDir[EdgeKey0]=DirFlatV0;
                assert(EdgeBorderDir.count(EdgeKey1)==0);
                EdgeBorderDir[EdgeKey1]= DirFlatV1;
            }

    }

    void InitConnections()
    {
        //the neighbours of each vertex (include propagation)
        std::vector<std::vector<size_t> > VNeigh;
        GetVertStar(VNeigh);


        //remove the paths that moves along the border
        RemoveOnBorderLinks(VNeigh);

        //then propagate neighbours
        for (size_t s=0;s<PropagationSteps;s++)
            PropagateNeigh(VNeigh);

        //allocate the nodes
        Nodes.clear();
        Nodes.resize(mesh.vert.size()*4);


        AddConnections(VNeigh);

        if (DebugMsg)
            std::cout<<"removing dead end links"<<std::endl;

        RemoveDeadEnd();

        //then add twins Connection along sharp features
        if (DebugMsg)
            std::cout<<"adding twins connections"<<std::endl;
        AddTwinsConnections();

        if (DebugMsg)
            std::cout<<"removing dead end links"<<std::endl;
        RemoveDeadEnd();

        for (size_t i=0;i<Nodes.size();i++)
            std::sort(Nodes[i].Neigh.begin(),Nodes[i].Neigh.end());

        InitDirectConnections();
        if (DebugMsg)
            std::cout<<"terminated initialization connections"<<std::endl;
    }

    std::vector<ScalarType> NodeDist;
    std::vector<int> NodeFather;
    std::vector<size_t> NodeJumps;
    std::vector<size_t> NodeTwinJumps;

public:

    std::vector<size_t> SingNodes;
    std::vector<bool> IsSingVert;
    bool remove_sign_connections;

    bool HasTwin(size_t IndexN)
    {
        return(Nodes[IndexN].HasTwin);
    }

    static size_t NodeDirI(size_t IndexN)
    {
        return(IndexN%4);
    }

    static void NodeDirI(const std::vector<size_t> &IndexN,
                         std::vector<size_t> &DirI)
    {
        DirI.clear();
        for (size_t i=0;i<IndexN.size();i++)
            DirI.push_back(NodeDirI(IndexN[i]));
    }

    static size_t NodeVertI(size_t IndexN)
    {
        return(IndexN/4);
    }

    static void NodeVertI(const std::vector<size_t> &IndexN,
                          std::vector<size_t> &IndexV)
    {
        IndexV.clear();
        for (size_t i=0;i<IndexN.size();i++)
            IndexV.push_back(NodeVertI(IndexN[i]));
    }

    static void VertDir(size_t IndexN,size_t &IndexV,size_t &IndexDir)
    {
        IndexV=NodeVertI(IndexN);
        IndexDir=NodeDirI(IndexN);
    }

    static size_t IndexNode(size_t IndexV,size_t IndexDir)
    {
        assert(IndexDir<4);
        return ((IndexV*4)+IndexDir);
    }

    static void IndexNodes(size_t IndexV,std::vector<size_t> &Nodes)
    {
        Nodes.clear();
        Nodes.resize(4);
        Nodes[0]=IndexV*4;
        Nodes[1]=IndexV*4+1;
        Nodes[2]=IndexV*4+2;
        Nodes[3]=IndexV*4+3;
    }

    static size_t TangentNode(size_t NodeIndex)
    {
        size_t IndexV,IndexDir;
        VertDir(NodeIndex,IndexV,IndexDir);

        //go to tangent
        IndexDir=(IndexDir+2)%4;
        return IndexNode(IndexV,IndexDir);
    }

    static void OrthoNode(size_t NodeIndex,
                          size_t &OrthoNode0,
                          size_t &OrthoNode1)
    {
        size_t IndexV,IndexDir;
        VertDir(NodeIndex,IndexV,IndexDir);

        //go to ortho
        size_t OrthoDir0=(IndexDir+1)%4;
        size_t OrthoDir1=(IndexDir+3)%4;
        OrthoNode0=IndexNode(IndexV,OrthoDir0);
        OrthoNode1=IndexNode(IndexV,OrthoDir1);
    }

    static void OrthoNodes(std::vector<size_t> &NodeIndexes,
                           std::vector<size_t> &OrthoNodes)
    {
        for (size_t i=0;i<NodeIndexes.size();i++)
        {
            size_t OrthoNode0,OrthoNode1;
            OrthoNode(NodeIndexes[i],OrthoNode0,OrthoNode1);
            OrthoNodes.push_back(OrthoNode0);
            OrthoNodes.push_back(OrthoNode1);
        }
        std::sort(OrthoNodes.begin(),OrthoNodes.end());
        std::vector<size_t>::iterator it;
        it = std::unique (OrthoNodes.begin(),OrthoNodes.end());
        OrthoNodes.resize( std::distance(OrthoNodes.begin(),it) );
    }

    static void TangentNodes(std::vector<size_t> &NodeIndexes)
    {
        for (size_t i=0;i<NodeIndexes.size();i++)
            NodeIndexes[i]=TangentNode(NodeIndexes[i]);
    }

    CoordType NodePos(size_t IndexN)const
    {
        size_t IndexV=NodeVertI(IndexN);
        assert(IndexV<mesh.vert.size());
        CoordType Pos=mesh.vert[IndexV].P();
        return Pos;
    }

    CoordType NodeDir(size_t IndexN)const
    {
        size_t IndexV,IndexDir;
        VertDir(IndexN,IndexV,IndexDir);
        return GetDirection(IndexV,IndexDir);
    }

    size_t NumNeigh(size_t IndexN)const
    {
        return(Nodes[IndexN].Neigh.size());
    }

    ScalarType AngleNeigh(const size_t &IndexN,
                          const size_t &IndexNeigh)const
    {
        assert(IndexNeigh<Nodes[IndexN].Neigh.size());
        return (Nodes[IndexN].Neigh[IndexNeigh].Angle);
    }

    ScalarType DistNeigh(const size_t &IndexN,
                         const size_t &IndexNeigh)const
    {
        assert(IndexNeigh<Nodes[IndexN].Neigh.size());
        return (Nodes[IndexN].Neigh[IndexNeigh].Dist);
    }

    bool TwinNeigh(const size_t &IndexN,
                   const size_t &IndexNeigh)const
    {
        assert(IndexNeigh<Nodes[IndexN].Neigh.size());
        return (Nodes[IndexN].Neigh[IndexNeigh].twin);
    }

    bool AreTwin(const size_t &IndexN0,
                 const size_t &IndexN1)const
    {
        if (NodePos(IndexN0)!=NodePos(IndexN1))
            return false;
        for (size_t i=0;i<NumNeigh(IndexN0);i++)
        {
            if (!TwinNeigh(IndexN0,i))continue;
            if (NodeNeigh(IndexN0,i)==IndexN1)return true;
        }
        return false;
    }

    bool ActiveNeigh(const size_t &IndexN,
                     const size_t &IndexNeigh)const
    {
        assert(IndexNeigh<Nodes[IndexN].Neigh.size());
        return (Nodes[IndexN].Neigh[IndexNeigh].Active);
    }

    size_t NodeNeigh(const size_t &IndexN,
                     const size_t &IndexNeigh)const
    {
        assert(IndexNeigh<Nodes[IndexN].Neigh.size());
        return (Nodes[IndexN].Neigh[IndexNeigh].Node);
    }

    void GetNodeNeigh(const size_t &IndexN,
                      std::vector<size_t> &NeighNodes)const
    {
        assert(IndexN<Nodes.size());
        NeighNodes.clear();
        for (size_t i=0;i<NumNeigh(IndexN);i++)
            NeighNodes.push_back(NodeNeigh(IndexN,i));
    }

    bool DirectNeigh(const size_t &IndexN,
                     const size_t &IndexNeigh)
    {
        assert(IndexNeigh<Nodes[IndexN].Neigh.size());
        return (Nodes[IndexN].Neigh[IndexNeigh].Direct);
    }

    void NodePosDir(size_t IndexN,CoordType &Pos,CoordType &Dir)
    {
        size_t IndexV,IndexDir;
        VertDir(IndexN,IndexV,IndexDir);
        Pos=mesh.vert[IndexV].P();
        Dir=GetDirection(IndexV,IndexDir);
    }

    MeshType &Mesh()
    {
        return mesh;
    }

    MeshType &Mesh()const
    {
        return mesh;
    }

//    void Init(bool _DebugMsg=false)//std::vector<CoordType> &_Sing)
//    {
//        DebugMsg=_DebugMsg;
//        //SingPos=std::set<CoordType>(_Sing.begin(),_Sing.end());
//        //check if everything is ok
//        CheckTangentField();
//        //initialize connections
//        InitConnections();
//        //initialize pos map
//        InitVertPosMap();
//        //initialize the rest used for queries
//        NodeDist=std::vector<ScalarType>(NumNodes(),0);
//        NodeFather=std::vector<int>(NumNodes(),-1);
//        NodeJumps=std::vector<size_t> (NumNodes(),0);
//        NodeTwinJumps=std::vector<size_t> (NumNodes(),0);

//        SingNodes.clear();
//        for (size_t i=0;i<mesh.vert.size();i++)
//        {
//            int Mmatch;
//            if(vcg::tri::CrossField<MeshType>::IsSingularByCross(mesh.vert[i],Mmatch))
//            {
//                std::vector<size_t> CurrN;
//                IndexNodes(i,CurrN);
//                SingNodes.insert(SingNodes.end(),CurrN.begin(),CurrN.end());
//            }
//        }
//        if (DebugMsg)
//            std::cout<<"There are "<<SingNodes.size()<<" singular Nodes"<<std::endl;
//        RemoveSingularities();
//    }

    void Reset()
    {
        SetAllActive();
        NodeDist.clear();
        NodeFather.clear();
        NodeJumps.clear();
        NodeTwinJumps.clear();
        SingNodes.clear();
        IsSingVert.clear();
        TMark=0;
        RealBorderVert.clear();
        EdgeBorderDir.clear();
        Nodes.clear();
        VertPos.clear();
    }

    void InitGraph(bool _DebugMsg=false)//std::vector<CoordType> &_Sing)
    {

        DebugMsg=_DebugMsg;
        //SingPos=std::set<CoordType>(_Sing.begin(),_Sing.end());
        //check if everything is ok
        CheckTangentField();

        //collect singularities nodes (must be disabled)
        SingNodes.clear();
        IsSingVert=std::vector<bool>(mesh.vert.size(),false);
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            int Mmatch;
            if(vcg::tri::CrossField<MeshType>::IsSingularByCross(mesh.vert[i],Mmatch,true))
            {
                IsSingVert[i]=true;
                std::vector<size_t> CurrN;
                IndexNodes(i,CurrN);
                SingNodes.insert(SingNodes.end(),CurrN.begin(),CurrN.end());
            }
        }
        if (DebugMsg)
            std::cout<<"There are "<<SingNodes.size()<<" singular Nodes"<<std::endl;
        //RemoveSingularities();


        //initialize connections
        InitConnections();

        //initialize pos map
        InitVertPosMap();

        //initialize Border direction map
        InitBorderDirMap();

        //initialize the rest used for queries
        NodeDist=std::vector<ScalarType>(NumNodes(),0);
        NodeFather=std::vector<int>(NumNodes(),-1);
        NodeJumps=std::vector<size_t> (NumNodes(),0);
        NodeTwinJumps=std::vector<size_t> (NumNodes(),0);
    }


    int &Father(const size_t &IndexNode)
    {
        assert(IndexNode<NodeFather.size());
        return (NodeFather[IndexNode]);
    }

    size_t &Jumps(const size_t &IndexNode)
    {
        assert(IndexNode<NodeJumps.size());
        return (NodeJumps[IndexNode]);
    }

//    ScalarType FieldL(const size_t &IndexNode0,
//                      const size_t &IndexNode1)
//    {
//        CoordType Dir0=NodeDir(IndexNode0);
//        CoordType Dir1=NodeDir(IndexNode1);
//        CoordType AvgDir=(Dir0+Dir1);
//        AvgDir.Normalize();
//        CoordType Pos0=NodePos(IndexNode0);
//        CoordType Pos1=NodePos(IndexNode1);
//        ScalarType L=(Pos1-Pos0)*AvgDir;
//        return L;
//        //        assert(IndexNode<NodeJumps.size());
//        //        return (NodeJumps[IndexNode]);
//    }

    size_t &TwinJumps(const size_t &IndexNode)
    {
        assert(IndexNode<NodeTwinJumps.size());
        return (NodeTwinJumps[IndexNode]);
    }

    ScalarType &Distance(const size_t &IndexNode)
    {
        assert(IndexNode<NodeDist.size());
        return (NodeDist[IndexNode]);
    }

    ScalarType Distance(const size_t &IndexNode)const
    {
        assert(IndexNode<NodeDist.size());
        return (NodeDist[IndexNode]);
    }

    std::vector<ScalarType> &Distances()
    {
        return (NodeDist);
    }

    bool IsMarked(const size_t &IndexNode)const
    {
        assert(IndexNode<Nodes.size());
        //assert(Nodes[IndexNode].Active);
        assert(Nodes[IndexNode].TMark<=TMark);
        return (Nodes[IndexNode].TMark==TMark);
    }

    bool IsSelected(const size_t &IndexNode)const
    {
        assert(IndexNode<Nodes.size());
        //assert(Nodes[IndexNode].Active);
        return (Nodes[IndexNode].Selected);
    }

    void ClearSelection()
    {
        for (size_t i=0;i<NumNodes();i++)
            Nodes[i].Selected=false;
    }

    void Select(const size_t &IndexNode)
    {
        assert(IndexNode<Nodes.size());
        //assert(Nodes[IndexNode].Active);
        Nodes[IndexNode].Selected=true;
    }

    size_t Select(const std::vector<bool> &ToSelect)
    {
        assert(ToSelect.size()==Nodes.size());
        size_t num=0;
        for (size_t i=0;i<ToSelect.size();i++)
        {
            Nodes[i].Selected=ToSelect[i];
            if (Nodes[i].Selected)num++;
        }
        return num;
    }


    void DeSelect(const size_t &IndexNode)
    {
        assert(IndexNode<Nodes.size());
        //assert(Nodes[IndexNode].Active);
        Nodes[IndexNode].Selected=false;
    }

    size_t NumSelected()
    {
        size_t NumSel=0;
        for (size_t i=0;i<NumNodes();i++)
            if (IsSelected(i))NumSel++;
        return NumSel;
    }

    void Mark(const size_t &IndexNode)
    {
        assert(IndexNode<Nodes.size());
        //assert(Nodes[IndexNode].Active);
        assert(Nodes[IndexNode].TMark<=TMark);
        Nodes[IndexNode].TMark=TMark;
    }

    void MarkAll()
    {
        for (size_t i=0;i<NumNodes();i++)
            Nodes[i].TMark=TMark;
    }

    size_t NumNodes()const
    {
        return Nodes.size();
    }

    void UnMarkAll()
    {
        TMark++;
    }

    //    bool IsRealBorderNode(const size_t &IndexNode)
    //    {
    //        if (IsBorder(IndexNode))return false;
    //        return(RealBorderVert[IndexNode]);
    //    }

    bool IsRealBorderVert(const size_t &IndexVert)
    {
        return(RealBorderVert[IndexVert]);
        //        if (!mesh.vert[IndexVert].IsB())return false;
        //        if (HasTwin(IndexNode))return true;
        //        return falsel
    }

    bool IsActive(const size_t &IndexNode)const
    {
        assert(IndexNode<Nodes.size());
        return(Nodes[IndexNode].Active);
    }

    void IsActiveNodes(std::vector<bool> &ActiveFlag)const
    {
        ActiveFlag.clear();
        for (size_t i=0;i<NumNodes();i++)
            ActiveFlag.push_back(IsActive(i));
    }

    void SetDisabledNodes(std::vector<bool> &DisabledFlag)
    {
        DisabledFlag.clear();
        for (size_t i=0;i<NumNodes();i++)
            SetActive(i,!DisabledFlag[i]);
    }

    void SetActiveNodes(std::vector<bool> &ActiveFlag)
    {
        ActiveFlag.clear();
        for (size_t i=0;i<NumNodes();i++)
            SetActive(i,ActiveFlag[i]);
    }

    //    bool IsSingularPos(const CoordType &pos)const
    //    {
    //        return(SingPos.count(pos)>0);
    //    }

    bool FirstValid(const size_t &IndexNode,size_t &NextNode)const
    {
        assert(IndexNode<Nodes.size());
        assert(IsActive(IndexNode));
        for (size_t i=0;i<NumNeigh(IndexNode);i++)
        {
            if (!ActiveNeigh(IndexNode,i))continue;
            NextNode=NodeNeigh(IndexNode,i);
            assert(NextNode>=0);
            assert(NextNode<Nodes.size());
            if (!IsActive(NextNode))continue;
            return true;
        }
        return false;
    }

//    void DeactivateInternalSingularities()
//    {
//        for (size_t i=0;i<SingNodes.size();i++)
//            SetActive(SingNodes[i],false);
//    }

//    void ActivateSingularities()
//    {
//        for (size_t i=0;i<SingNodes.size();i++)
//            SetActive(SingNodes[i],true);
//    }

    void RemoveConnections(const size_t &IndexNode)
    {
        Nodes[IndexNode].Neigh.clear();
        for (size_t i=0;i<Nodes.size();i++)
        {
            std::vector<NeighInfo> SwapNeigh;
            for (size_t j=0;j<Nodes[i].Neigh.size();j++)
            {
                if (Nodes[i].Neigh[j].Node==IndexNode)continue;
                SwapNeigh.push_back(Nodes[i].Neigh[j]);
            }
            if (SwapNeigh.size()==Nodes[i].Neigh.size())continue;
            Nodes[i].Neigh=SwapNeigh;
        }
    }

    void RemoveConnections(const std::set<size_t> &IndexNodes)
    {
        std::set<size_t>::iterator IteSet;
        for (IteSet=IndexNodes.begin();IteSet!=IndexNodes.end();IteSet++)
            Nodes[(*IteSet)].Neigh.clear();

        for (size_t i=0;i<Nodes.size();i++)
        {
            std::vector<NeighInfo> SwapNeigh;
            for (size_t j=0;j<Nodes[i].Neigh.size();j++)
            {
                if (IndexNodes.count(Nodes[i].Neigh[j].Node>0))continue;
                SwapNeigh.push_back(Nodes[i].Neigh[j]);
            }
            if (SwapNeigh.size()==Nodes[i].Neigh.size())continue;
            Nodes[i].Neigh=SwapNeigh;
        }
    }

//    void RemoveSingularities()
//    {
//        for (size_t i=0;i<SingNodes.size();i++)
//            RemoveConnections(SingNodes[i]);
//    }

    void SetActive(const size_t &IndexNode,bool ActiveVal)
    {
        assert(IndexNode<Nodes.size());
        Nodes[IndexNode].Active=ActiveVal;
    }

    void SetAllActive()
    {
        for (size_t i=0;i<Nodes.size();i++)
            Nodes[i].Active=true;
    }

    vcg::face::Pos<FaceType> GetNodesPos(const size_t &IndexN0,const size_t &IndexN1)const
    {
        size_t IndexV0=NodeVertI(IndexN0);
        size_t IndexV1=NodeVertI(IndexN1);
        std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),
                                     std::max(IndexV0,IndexV1));
        assert(VertPos.count(key)>0);
        return (VertPos.at(key));
        //return(VertPos[key]);
    }

    void GetNodesPos(const std::vector<size_t> &PathNodes,bool IsLoop,
                     std::vector<vcg::face::Pos<FaceType> > &PathPos)const
    {

        PathPos.clear();
        if (PathNodes.size()==0)return;
        assert(PathNodes.size()>=2);
        //retrieve the indexes
        size_t Limit=PathNodes.size()-1;
        if (IsLoop)Limit++;
        for (size_t j=0;j<Limit;j++)
        {
            size_t IndexN0=PathNodes[j];
            size_t IndexN1=PathNodes[(j+1)%PathNodes.size()];
            CoordType Pos0=NodePos(IndexN0);
            CoordType Pos1=NodePos(IndexN1);
            //check if twin
            if (Pos0==Pos1)continue;
            vcg::face::Pos<FaceType> CurrPos=GetNodesPos(PathNodes[j],PathNodes[(j+1)%PathNodes.size()]);
            PathPos.push_back(CurrPos);
        }
    }

    //*** CONSTRUCTORS ***
    //draw a trace
    VertexFieldGraph(MeshType &_mesh,size_t _PropagationSteps=1):mesh(_mesh)
    {
        PropagationSteps=_PropagationSteps;
        assert(PropagationSteps>=1);
        TMark=0;
        remove_sign_connections=true;
    }
};



template <class MeshType>
void SplitAdjacentSingularities(MeshType &mesh)
{
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef std::pair<CoordType,CoordType> CoordPair;

    // Basic subdivision class
    struct SplitLev : public   std::unary_function<vcg::face::Pos<FaceType> ,CoordType >
    {
        std::map<CoordPair,CoordType> *SplitOps;

        void operator()(VertexType &nv,vcg::face::Pos<FaceType>  ep)
        {
            VertexType* v0=ep.f->V0(ep.z);
            VertexType* v1=ep.f->V1(ep.z);

            assert(v0!=v1);

            CoordType Pos0=v0->P();
            CoordType Pos1=v1->P();

            CoordPair CoordK(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
            assert(SplitOps->count(CoordK)>0);
            nv.P()=(*SplitOps)[CoordK];
        }

        vcg::TexCoord2<ScalarType> WedgeInterp(vcg::TexCoord2<ScalarType> &t0,
                                               vcg::TexCoord2<ScalarType> &t1)
        {
            (void)t0;
            (void)t1;
            return (vcg::TexCoord2<ScalarType>(0,0));
        }

        SplitLev(std::map<CoordPair,CoordType> *_SplitOps){SplitOps=_SplitOps;}
    };

    class EdgePred
    {

        std::map<CoordPair,CoordType> *SplitOps;

    public:

        bool operator()(vcg::face::Pos<FaceType> ep) const
        {
            VertexType* v0=ep.f->V0(ep.z);
            VertexType* v1=ep.f->V1(ep.z);

            assert(v0!=v1);

            CoordType Pos0=v0->P();
            CoordType Pos1=v1->P();

            CoordPair CoordK(std::min(Pos0,Pos1),std::max(Pos0,Pos1));

            return (SplitOps->count(CoordK)>0);
        }

        EdgePred(std::map<CoordPair,CoordType> *_SplitOps){SplitOps=_SplitOps;}
    };


    //InitFeatureCoordsTable();
    vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);
    for (size_t i=0;i<mesh.vert.size();i++)
    {
        int Mmatch;
        if(!vcg::tri::CrossField<MeshType>::IsSingularByCross(mesh.vert[i],Mmatch,true))continue;
        mesh.vert[i].SetS();
    }

    //then also split single sharp edges
    std::vector<size_t> VertCreases(mesh.vert.size(),0);
    for (size_t i=0;i<mesh.face.size();i++)
    {
        for (size_t j=0;j<3;j++)
        {
            size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
            size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
            if (!mesh.face[i].IsFaceEdgeS(j))continue;
            if (vcg::face::IsBorder(mesh.face[i],j))continue;
            if (mesh.face[i].FFp(j)>&mesh.face[i])continue;
            VertCreases[IndexV0]++;
            VertCreases[IndexV1]++;
        }
    }

    std::map<CoordPair,CoordType> ToBeSplitted;
    std::vector<CoordPair> Creases;
    for (size_t i=0;i<mesh.face.size();i++)
    {
        for (size_t j=0;j<3;j++)
        {
            CoordType P0=mesh.face[i].P0(j);
            CoordType P1=mesh.face[i].P1(j);
            std::pair<CoordType,CoordType> key(std::min(P0,P1),std::max(P0,P1));
            if (mesh.face[i].IsFaceEdgeS(j))
                Creases.push_back(key);

            bool SplitForSing=((mesh.face[i].V0(j)->IsS())||(mesh.face[i].V1(j)->IsS()));
            size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
            size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
            bool SplitForSingleE=((VertCreases[IndexV0]==1)&&(VertCreases[IndexV1]==1));
            if (SplitForSing || SplitForSingleE)
            {
                CoordType Avg=(P0+P1)/2;
                ToBeSplitted[key]=Avg;
                if (mesh.face[i].IsFaceEdgeS(j))
                {
                    std::pair<CoordType,CoordType> Key0(std::min(P0,Avg),std::max(P0,Avg));
                    std::pair<CoordType,CoordType> Key1(std::min(P1,Avg),std::max(P1,Avg));
                    Creases.push_back(Key0);
                    Creases.push_back(Key1);
                }
            }
        }
    }

    SplitLev splMd(&ToBeSplitted);
    EdgePred eP(&ToBeSplitted);

    //do the final split
    vcg::tri::RefineE<MeshType,SplitLev,EdgePred>(mesh,splMd,eP);
    mesh.UpdateFromCoordPairs(Creases,false);

//    for (size_t i=0;i<mesh.face.size();i++)
//    {
//        for (size_t j=0;j<3;j++)
//        {
//            CoordType P0=mesh.face[i].P0(j);
//            CoordType P1=mesh.face[i].P1(j);
//            std::pair<CoordType,CoordType> key(std::min(P0,P1),std::max(P0,P1));
//            if (Creases.count(key)==0);
//            mesh.face[i].SetFaceEdgeS(j);
//        }
//    }

    mesh.UpdateAttributes();
    //    SetFeatureFromTable();
    //    return done;
}

template <class MeshType>
void PreProcessMesh(MeshType &mesh,bool DebugMsg=true)
{


    mesh.SelectSharpFeatures();


    SplitAdjacentSingularities(mesh);


    //split along marked sharp features
    if (DebugMsg)
        std::cout<<"splitting along sharp features"<<std::endl;
    VertSplitter<MeshType>::SplitAlongEdgeSel(mesh);
    if (DebugMsg)
        std::cout<<"done"<<std::endl;

    if (DebugMsg)
        std::cout<<"updating attributes"<<std::endl;

    mesh.UpdateAttributes();

    if (DebugMsg)
        std::cout<<"compact vectors"<<std::endl;

    vcg::tri::Allocator<MeshType>::CompactEveryVector(mesh);


    if (DebugMsg)
        std::cout<<"counting non manifold Vert"<<std::endl;

    size_t Test1=vcg::tri::Clean<MeshType>::CountNonManifoldVertexFF(mesh);
    if (Test1>0)
    {
        std::cout<<"WARNING NON MANIFOLD VERTEX SPLIT! "<<std::endl;
        vcg::tri::Clean<MeshType>::SplitNonManifoldVertex(mesh,std::numeric_limits<typename MeshType::ScalarType>::epsilon()*100);
    }
    if (DebugMsg)
        std::cout<<"splitted"<<std::endl;
    //then update attributes
    mesh.UpdateAttributes();

    size_t Test2=vcg::tri::Clean<MeshType>::CountNonManifoldVertexFF(mesh);
    if (DebugMsg)
        std::cout<<"Non Manif "<<Test2<<std::endl;
    //vcg::tri::io::ExporterPLY<MeshType>::Save(mesh,"test0.ply");
    //then reupdate the vert cross field
    vcg::tri::CrossField<MeshType>::UpdateSingularByCross(mesh,true);
    vcg::tri::CrossField<MeshType>::SetVertCrossVectorFromFace(mesh);

    if (DebugMsg)
        std::cout<<"setting rest pos"<<std::endl;

    for (size_t i=0;i<mesh.vert.size();i++)
        mesh.vert[i].RPos=mesh.vert[i].P();
    for (size_t i=0;i<mesh.face.size();i++)
        mesh.face[i].FullTraced=false;

    if (DebugMsg)
        std::cout<<"done"<<std::endl;

    mesh.InitSingVert();
}

#endif
