#ifndef QR_FIELD_TRACER_H
#define QR_FIELD_TRACER_H

#include <vcg/math/matrix33.h>
#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>
#include <vcg/complex/algorithms/isotropic_remeshing.h>
#include <vcg/complex/algorithms/geodesic.h>
#include <wrap/io_trimesh/export.h>

namespace QuadRetopology {
namespace internal {

template <class MeshType>
class VertexFieldTracer
{
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;

    //the input mesh to trace
    MeshType &mesh;

    //neighbours only connected by original edges
    std::vector<std::vector<size_t> > VNeigh1;
    //propagation steps
    size_t PropagationSteps;

    struct NeighInfo
    {
        size_t NextV;
        size_t NextDir;
        ScalarType NextWeight;

        inline bool operator <(const NeighInfo &NInfo)const
        {
            return (NextWeight>NInfo.NextWeight);
        }

        NeighInfo(size_t &_NextV,
                  size_t &_NextDir,
                  ScalarType &_NextWeight)
        {
            NextV=_NextV;
            NextDir=_NextDir;
            NextWeight=_NextWeight;
        }
    };


public:

    //per vertex per directions neighbors
    std::vector<std::vector<std::vector<NeighInfo> > > VertNeigh;

    void ConnectionData(const size_t &IndexV0,
                        const size_t &DirectionV0,
                        const size_t &IndexV1,
                        ScalarType &Weight,
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
        Weight=(TargetD*TestDir);
        Direction_next=FollowDirection(IndexV0,IndexV1,DirectionV0);
    }

    void GetEdgeNeightbours(std::vector<std::vector<size_t> > &VNeigh)
    {
        VNeigh.clear();
        VNeigh.resize(mesh.vert.size());
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<mesh.face[i].VN();j++)
            {
                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
                VNeigh[IndexV0].push_back(IndexV1);
                VNeigh[IndexV1].push_back(IndexV0);
            }
    }

    void RemoveOnBorderPaths(std::vector<std::vector<size_t> > &VNeigh)
    {
        //first save all the borders one
        std::set<std::pair<size_t,size_t> > BorderEdges;
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<mesh.face[i].VN();j++)
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

    void RemoveNodes(std::vector<std::vector<size_t> > &VNeigh,
                     const std::set<size_t> &ExcludeNodes)
    {
        for (size_t i=0;i<VNeigh.size();i++)
        {
            if (ExcludeNodes.count(i)>0)
            {
                VNeigh[i].clear();
                continue;
            }
            std::vector<size_t> NewNeigh;
            for (size_t j=0;j<VNeigh[i].size();j++)
            {
                size_t IndexVN=VNeigh[i][j];
                if (ExcludeNodes.count(IndexVN)>0)continue;
                NewNeigh.push_back(IndexVN);
            }
            VNeigh[i]=NewNeigh;
        }
    }

    void InitConnections(const std::set<size_t> &ExcludeNodes)
    {
        //the neighbours of each vertex (include propagation)
        std::vector<std::vector<size_t> > VNeigh;
        GetEdgeNeightbours(VNeigh);
        VNeigh1=VNeigh;

        RemoveOnBorderPaths(VNeigh);

        for (size_t s=0;s<PropagationSteps;s++)
            PropagateNeigh(VNeigh);

        RemoveNodes(VNeigh,ExcludeNodes);

        VertNeigh.clear();
        VertNeigh.resize(mesh.vert.size());

//        std::cout<<"initializing connections"<<std::endl;
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            //then set per direction
            VertNeigh[i].resize(4);

            for (size_t j=0;j<VNeigh[i].size();j++)
            {
                size_t IndexV0=i;
                size_t IndexV1=VNeigh[i][j];
                for (size_t k=0;k<4;k++)
                {
                    ScalarType Weight;
                    size_t NextDir;
                    ConnectionData(IndexV0,k,IndexV1,Weight,NextDir);

                    if (Weight<=0)continue;

                    VertNeigh[i][k].push_back(NeighInfo(IndexV1,NextDir,Weight));
                }
            }
        }


//        std::cout<<"removing unconnected nodes"<<std::endl;

        bool Removed=false;
        size_t step=0;
        do
        {
            Removed=false;
            step++;
            for (size_t i=0;i<VertNeigh.size();i++)
            {
                for (size_t k=0;k<4;k++)
                {
                    std::vector<bool> MustRemoved(VertNeigh[i][k].size(),false);
                    bool remove_session=false;
                    for (size_t j=0;j<VertNeigh[i][k].size();j++)
                    {
                        size_t NextV=VertNeigh[i][k][j].NextV;
                        size_t NextD=VertNeigh[i][k][j].NextDir;
                        if (mesh.vert[NextV].IsB())continue;

                        bool DeadEnd=(VertNeigh[NextV][NextD].size()==0);
                        if (!DeadEnd)continue;
                        MustRemoved[j]=true;
                        remove_session=true;
                        Removed=true;
                    }
                    if (remove_session)
                    {
                        std::vector<NeighInfo> Swap;
                        for (size_t j=0;j<MustRemoved.size();j++)
                        {
                            if (MustRemoved[j])continue;
                            Swap.push_back(VertNeigh[i][k][j]);
                        }
                        VertNeigh[i][k]=Swap;
                    }
                }
            }
        }while (Removed);

//        std::cout<<"permormed removed steps "<<step<<std::endl;
        for (size_t i=0;i<VertNeigh.size();i++)
            for (size_t j=0;j<VertNeigh[i].size();j++)
            {
                std::sort(VertNeigh[i][j].begin(),VertNeigh[i][j].end());
            }

//        std::cout<<"terminated initialization connections"<<std::endl;
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
                    std::cout<<"WARNING NON TANGENT FIELD 1"<<std::endl;
            }
        }
    }

    //*** TRACING FUNCTIONS ***
    //return a direction given an index 0<=Direction<=3 for a given vertex
    CoordType GetDirection(const size_t IndexV,
                           const size_t Direction)
    {
        assert(Direction>=0);
        assert(Direction<4);
        if (Direction==0)
            return mesh.vert[IndexV].PD1();
        if (Direction==1)
            return mesh.vert[IndexV].PD2();
        if (Direction==2)
            return (-mesh.vert[IndexV].PD1());
        if (Direction==3)
            return (-mesh.vert[IndexV].PD2());
    }

    //return the direction for a vertex that is closest to an
    //input direction expressed as a CoordType
    size_t GetClosestDirTo(const size_t IndexV0,
                           const CoordType Dir)
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

    //follow the direction from one vertex to another
    //it return the index of the direction in IndexV1
    bool TraceNext(const size_t IndexV0,
                   const size_t DirectionV0,
                   size_t &IndexV_next,
                   size_t &Direction_next)
    {
        assert(IndexV0>=0);
        assert(IndexV0<mesh.vert.size());
        if (VertNeigh[IndexV0][DirectionV0].size()==0)return false;
        IndexV_next=VertNeigh[IndexV0][DirectionV0][0].NextV;
        Direction_next=VertNeigh[IndexV0][DirectionV0][0].NextDir;
        return true;
    }

    //trace from one Index of Vertex toward a direction
    //it can be also tuned to stop when it encounter a selected vertex
    bool TraceFrom(const size_t IndexV0,
                   const size_t DirectionV0,
                   std::vector<size_t> &IndexV,
                   std::vector<size_t> &IndexDir,
                   bool StopAtSel=false)
    {
        size_t IndexVCurr=IndexV0;
        size_t DirectionVCurr=DirectionV0;

        IndexV.clear();
        IndexDir.clear();
        IndexV.push_back(IndexV0);
        IndexDir.push_back(DirectionV0);
        std::set<std::pair<size_t,size_t> >Traced;
        bool has_terminated=false;
        do {
            Traced.insert(std::pair<size_t,size_t>(IndexVCurr,DirectionVCurr%2));
            size_t IndexVNext;
            size_t DirectionVNext;
            bool traced=TraceNext(IndexVCurr,DirectionVCurr,IndexVNext,DirectionVNext);
            if (!traced)
            {
//                std::cout<<"stopped tracing"<<std::endl;
                return false;
            }
            if (Traced.count(std::pair<size_t,size_t>(IndexVNext,DirectionVNext%2))>0)
            {
//                std::cout<<"tangent self intersecting"<<std::endl;
                return false;
            }
            //std::cout<<"next step"<<std::endl;

            assert(IndexVNext!=IndexVCurr);
            IndexVCurr=IndexVNext;
            DirectionVCurr=DirectionVNext;

            IndexV.push_back(IndexVCurr);
            IndexDir.push_back(DirectionVCurr);

            has_terminated=mesh.vert[IndexVCurr].IsB();
            if (StopAtSel)
                has_terminated|=mesh.vert[IndexVCurr].IsS();

        }while(!has_terminated);
        return true;
    }

    //    //*** INITIALIZATION FUNCTIONS ***
    //    //initialize neighbours
    //    void InitVertNeigh()
    //    {
    //        VNeigh.clear();
    //        VNeigh.resize(mesh.vert.size());
    //        for (size_t i=0;i<mesh.face.size();i++)
    //            for (size_t j=0;j<mesh.face[i].VN();j++)
    //            {
    //                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
    //                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
    //                VNeigh[IndexV0].push_back(IndexV1);
    //                VNeigh[IndexV1].push_back(IndexV0);
    //            }
    //        VNeigh1=VNeigh;
    //    }

    //propagate neighbours
    //    void PropagateNeigh()
    //    {
    //        std::vector<std::vector<size_t> > newNeigh;

    //        for (size_t i=0;i<VNeigh.size();i++)
    //        {
    //            newNeigh.push_back(VNeigh[i]);
    //            for (size_t j=0;j<VNeigh[i].size();j++)
    //            {
    //                size_t IndexN=VNeigh[i][j];
    //                if (mesh.vert[IndexN].IsB()) continue;
    //                for (size_t k=0;k<VNeigh[IndexN].size();k++)
    //                {
    //                    if (mesh.vert[VNeigh[IndexN][k]].IsB())continue;
    //                    if (VNeigh[IndexN][k]==i)continue;
    //                    newNeigh.back().push_back(VNeigh[IndexN][k]);
    //                }
    //            }

    //            std::sort(newNeigh.back().begin(),newNeigh.back().end());
    //            std::vector<size_t>::iterator it;
    //            it = std::unique (newNeigh.back().begin(),newNeigh.back().end());   // 10 20 30 20 10 ?  ?  ?  ?
    //            newNeigh.back().resize( std::distance(newNeigh.back().begin(),it) ); // 10 20 30 20 10
    //        }
    //        VNeigh=newNeigh;
    //    }

    //    //remove the neighbours that pass along the border
    //    void RemoveOnBorderPaths()
    //    {
    //        //first save all the borders one
    //        std::set<std::pair<size_t,size_t> > BorderEdges;
    //        for (size_t i=0;i<mesh.face.size();i++)
    //            for (size_t j=0;j<mesh.face[i].VN();j++)
    //            {
    //                if (!vcg::face::IsBorder(mesh.face[i],j))continue;
    //                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
    //                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
    //                BorderEdges.insert(std::pair<size_t,size_t>(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1)));
    //            }

    //        for (size_t i=0;i<VNeigh.size();i++)
    //        {
    //            size_t IndexV0=i;
    //            std::vector<size_t> NewNeigh;
    //            for (size_t j=0;j<VNeigh[i].size();j++)
    //            {
    //                size_t IndexV1=VNeigh[i][j];
    //                std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
    //                if (BorderEdges.count(key)>0)continue;
    //                NewNeigh.push_back(IndexV1);
    //            }
    //            VNeigh[i]=NewNeigh;
    //        }
    //    }

//    void RemoveNodes(const std::set<size_t> &ExcludeNodes)
//    {
//        for (size_t i=0;i<VNeigh.size();i++)
//        {
//            if (ExcludeNodes.count(i)>0)
//            {
//                VNeigh[i].clear();
//                continue;
//            }
//            std::vector<size_t> NewNeigh;
//            for (size_t j=0;j<VNeigh[i].size();j++)
//            {
//                size_t IndexVN=VNeigh[i][j];
//                if (ExcludeNodes.count(IndexVN)>0)continue;
//                NewNeigh.push_back(IndexVN);
//            }
//            VNeigh[i]=NewNeigh;
//        }
//    }

    void Init(const std::set<size_t> &ExcludeNodes)
    {
        //        InitVertNeigh();
        //        for (size_t s=0;s<PropagationSteps;s++)
        //            PropagateNeigh();

        assert(0);

        InitConnections(ExcludeNodes);

        //RemoveOnBorderPaths();

        //RemoveNodes(ExcludeNodes);

        CheckTangentField();

        //InitConnectionWeights();
    }

    //*** LENGHT FUNCTIONS ***
    //return the lenght of an edge, considering the field
    void EdgeLengh(const size_t IndexV0,
                   const size_t IndexV1,
                   const size_t FieldDirIndex0,
                   CoordType &Dir,
                   ScalarType &Len)
    {
        CoordType edgedir=(mesh.vert[IndexV0].P()-
                           mesh.vert[IndexV1].P());

        Dir=GetDirection(IndexV0,FieldDirIndex0);
        Len=fabs(edgedir*Dir);
    }

    //return the lenght of an trace, considering the field
    ScalarType TraceLenght(std::vector<size_t> &CurrTrace,
                           std::vector<size_t> &DirTrace)
    {
        ScalarType Len=0;
        for (size_t i=0;i<CurrTrace.size()-1;i++)
        {
            ScalarType currL;
            CoordType Dir;
            EdgeLengh(CurrTrace[i],CurrTrace[i+1],DirTrace[i],Dir,currL);
            Len+=currL;
        }
        return Len;
    }

    bool CollideSubSequence(const std::vector<size_t > &VIndex0,
                            const std::vector<size_t > &VIndex1)
    {
        assert(VIndex0.size()==3);
        assert(VIndex1.size()==3);
        if (VIndex0[1]!=VIndex1[1])return false;

        //get a first face
        VertexType *vMiddle=&mesh.vert[VIndex0[1]];
        vcg::face::VFIterator<FaceType> VfI(vMiddle);
        FaceType *F0 = VfI.F();
        size_t IndexV= VfI.I();
        assert(F0->V(IndexV)==vMiddle);

        //initialize a pos
        vcg::face::Pos<FaceType> InitPos(F0,IndexV);
        std::vector<vcg::face::Pos<FaceType> > PosStar;
        vcg::face::VFOrderedStarFF(InitPos,PosStar);

        int IndexPartition=-1;
        for (size_t i=0;i<PosStar.size();i++)
        {
            assert(PosStar[i].V()==vMiddle);
            VertexType *vOpp=PosStar[i].VFlip();
            size_t IndexVOpp=vcg::tri::Index(mesh,vOpp);

            if ((IndexVOpp==VIndex0[0])||(IndexVOpp==VIndex0[2]))
            {
                if (IndexPartition==-1)
                    IndexPartition=0;
                else
                    if (IndexPartition==0)
                        return true;
                    else
                        if (IndexPartition==1)
                            IndexPartition=0;

            }

            if ((IndexVOpp==VIndex1[0])||(IndexVOpp==VIndex1[2]))
            {
                if (IndexPartition==-1)
                    IndexPartition=1;
                else
                    if (IndexPartition==1)
                        return true;
                    else
                        if (IndexPartition==0)
                            IndexPartition=1;
            }
        }
        return false;
        //vcg::face::Pos<FaceType> currPos(F0,IndexV);
    }

    //*** COLLIDING FUNCTIONS ***
    //return true if two traces collides
    bool CollideTraces(const std::vector<size_t > &TraceV0,
                       const std::vector<size_t > &TraceDir0,
                       const std::vector<size_t > &TraceV1,
                       const std::vector<size_t > &TraceDir1)
    {
        assert(TraceV0.size()==TraceDir0.size());
        assert(TraceV1.size()==TraceDir1.size());

        //not consider last one, they can intersect there!
        std::set<std::pair<size_t,size_t> > PathNodes0;
        for (size_t i=1;i<TraceV0.size()-1;i++)
            PathNodes0.insert(std::pair<size_t,size_t>(TraceV0[i],TraceDir0[i]%2));

        for (size_t i=1;i<TraceV1.size()-1;i++)
        {
            std::pair<size_t,size_t> key(TraceV1[i],TraceDir1[i]%2);
            if (PathNodes0.count(key)>0)return true;
        }

        std::set<std::pair<size_t,size_t> > PathEdges0;
        for (size_t i=0;i<TraceV0.size()-1;i++)
            PathEdges0.insert(std::pair<size_t,size_t>(std::min(TraceV0[i],TraceV0[i+1]),
                              std::max(TraceV0[i],TraceV0[i+1])));

        for (size_t i=0;i<TraceV1.size()-1;i++)
        {
            std::pair<size_t,size_t> key=std::pair<size_t,size_t>(std::min(TraceV1[i],TraceV1[i+1]),
                    std::max(TraceV1[i],TraceV1[i+1]));
            if (PathEdges0.count(key)>0)return true;
        }


        for (size_t i=1;i<TraceV0.size()-1;i++)
        {
            std::vector<size_t> VIndex0;
            VIndex0.push_back(TraceV0[i-1]);
            VIndex0.push_back(TraceV0[i]);
            VIndex0.push_back(TraceV0[i+1]);
            for (size_t j=1;j<TraceV1.size()-1;j++)
            {
                std::vector<size_t> VIndex1;
                VIndex1.push_back(TraceV1[j-1]);
                VIndex1.push_back(TraceV1[j]);
                VIndex1.push_back(TraceV1[j+1]);
                if (CollideSubSequence(VIndex0,VIndex1))return true;
            }
        }

        return false;
    }

    //*** EXPANDING SEQUENCE FUNCTIONS ***
    //return the lenght of an edge, considering the field
    struct HeapEntry
    {
        ScalarType Dist;
        size_t Vert;
        bool operator < (const HeapEntry &He) const
        {
            if( Dist != He.Dist)
                return Dist > He.Dist;
            return Vert<He.Vert;
        }

        HeapEntry(size_t _Vert,
                  ScalarType _Dist)
        {
            Dist=_Dist;
            Vert=_Vert;
        }
    };

    //return the subsequence that connect two nodes, used to expand paths
    bool GetSubSequence(const size_t IndexV0,
                        const size_t IndexV1,
                        std::vector<size_t> &IndexV)
    {

        std::map<size_t,ScalarType> VertDist;
        std::map<size_t,size_t> VertJump;
        std::map<size_t,size_t> Father;
        std::vector<HeapEntry> Heap;

        VertDist[IndexV0]=0;
        VertJump[IndexV0]=0;
        Father[IndexV0]=IndexV0;

        Heap.push_back(HeapEntry(IndexV0,0));
        std::make_heap(Heap.begin(),Heap.end());

        do
        {
            pop_heap(Heap.begin(),Heap.end());
            size_t CurrV=(Heap.back()).Vert;
            ScalarType CurrDist=(Heap.back()).Dist;
            Heap.pop_back();

            if (CurrV==IndexV1)
            {
                IndexV.clear();
                //retrieve the sequence
                do {
                    IndexV.push_back(CurrV);
                    assert(Father.count(CurrV)>0);
                    CurrV=Father[CurrV];
                }while (CurrV!=Father[CurrV]);
                IndexV.push_back(CurrV);
                std::reverse(IndexV.begin(),IndexV.end());
                //std::cout<<"found at "<<VertJump[CurrV]<<std::endl;
                return true;
            }
            int CurrJump=VertJump[CurrV];
            //std::cout<<"test at "<<CurrJump<<std::endl;
            if (CurrJump==(pow(PropagationSteps,2)+1))continue;//maximum number of jumps reached
            //assert(VertDist.count(CurrV)>0);
            //ScalarType CurrDist=VertDist[CurrV];

            for (size_t i=0;i<VNeigh1[CurrV].size();i++)
            {
                size_t NextV=VNeigh1[CurrV][i];

                //internal vertices cannot be on borders
                if ((NextV!=IndexV1)&&(mesh.vert[NextV].IsB()))continue;

                assert(CurrV!=NextV);
                ScalarType EdgeLen=(mesh.vert[CurrV].P()-
                                    mesh.vert[NextV].P()).Norm();
                ScalarType NextDist=CurrDist+EdgeLen;
                if ((VertDist.count(NextV)==0)||(VertDist[NextV]>NextDist))
                {
                    VertDist[NextV]=NextDist;
                    //std::cout<<"next dist "<<NextDist<<std::endl;
                    Father[NextV]=CurrV;
                    VertJump[NextV]=CurrJump+1;
                    Heap.push_back(HeapEntry(NextV,NextDist));
                    push_heap(Heap.begin(),Heap.end());
                }
            }
        }while (!Heap.empty());
        return false;
    }

    //return the subsequence that connect two nodes, used to expand paths
    //also return the directions
    bool GetSubSequence(const size_t IndexV0,
                        const size_t DirectionV0,
                        const size_t IndexV1,
                        const size_t DirectionV1,
                        std::vector<size_t> &IndexV,
                        std::vector<size_t> &DirV)
    {

        bool found=GetSubSequence(IndexV0,IndexV1,IndexV);

        if (!found)return false;
        DirV.push_back(DirectionV0);

        if (IndexV.size()>2)
            for (size_t i=0;i<IndexV.size()-2;i++)
                DirV.push_back(FollowDirection(IndexV0,IndexV[i],DirectionV0));

        DirV.push_back(DirectionV1);
        //        std::cout<<"size 0 "<<IndexV.size()<<std::endl;
        //        std::cout<<"size 1 "<<DirV.size()<<std::endl;
        return true;
    }

    //used to expand a path
    bool ExpandPath(std::vector<size_t > &TraceVert,
                    std::vector<size_t > &TraceDir)
    {
        std::vector<size_t > SwapTraceVert;
        std::vector<size_t > SwapTraceDir;
        for (size_t i=0;i<TraceVert.size()-1;i++)
        {
            size_t IndexV0=TraceVert[i];
            size_t DirV0=TraceDir[i];
            size_t IndexV1=TraceVert[i+1];
            size_t DirV1=TraceDir[i+1];
            std::vector<size_t> IndexV;
            std::vector<size_t> DirV;
            bool found=GetSubSequence(IndexV0,DirV0,IndexV1,DirV1,IndexV,DirV);
            if(!found)return false;
            assert(IndexV.size()==DirV.size());
            assert(IndexV.size()>=2);
            SwapTraceVert.insert(SwapTraceVert.end(),IndexV.begin(),IndexV.end()-1);
            SwapTraceDir.insert(SwapTraceDir.end(),DirV.begin(),DirV.end()-1);
        }
        SwapTraceVert.push_back(TraceVert.back());
        SwapTraceDir.push_back(TraceDir.back());
        TraceVert=SwapTraceVert;
        TraceDir=SwapTraceDir;
        return true;
    }

#ifdef DRAWTRACE
    //*** DRAWING FUNCTIONS ***
    //draw a trace
    void DrawTrace(std::vector<size_t> &CurrTrace,
                   vcg::Color4b TraceCol)
    {
        vcg::glColor(TraceCol);

        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.99998);
        glLineWidth(12);
        glBegin(GL_LINE_STRIP);
        for (size_t i=0;i<CurrTrace.size();++i)
            vcg::glVertex(mesh.vert[CurrTrace[i]].P());
        glEnd();

        glPopAttrib();
    }
#endif

    //*** CONSTRUCTORS ***
    //draw a trace
    VertexFieldTracer(MeshType &_mesh,size_t _PropagationSteps=1):mesh(_mesh)
    {
        PropagationSteps=_PropagationSteps;
        assert(PropagationSteps>=1);
    }
};

}
}

#endif
