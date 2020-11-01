#ifndef PATCH_TRACER
#define PATCH_TRACER


#include "vert_field_graph.h"
#include "graph_query.h"
#include <vcg/math/matrix33.h>
#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>
#include <vcg/complex/algorithms/point_sampling.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/algorithms/closest.h>
#include <vcg/complex/algorithms/symmetry.h>
#include <wrap/io_trimesh/export_ply.h>
#include <wrap/io_trimesh/export_obj.h>
#include "vertex_emitter.h"
#include "vertex_classifier.h"
#include <tracing/candidate_path.h>
#include <tracing/patch_manager.h>
#include <tracing/patch_subdivider.h>

#define MAX_SAMPLES 1000
#define MAX_NARROW_CONST 0.05
#define NARROW_NEED 1
#define MAX_TWIN_DIKSTRA 1
#define MIN_ADMITTIBLE 3
#define MAX_ADMITTIBLE 6

//enum TraceType{TraceDirect,DijkstraReceivers,TraceLoop};
enum PatchType{LowCorners,HighCorners,NonDisk,
               HasEmitter,MaxCClarkability,NonMatchValence,
               IsOK};//MoreSing,IsOK};


template <class MeshType>
void RetrievePatchFromSelEdges(MeshType &mesh,
                               const size_t &IndexF,
                               std::vector<size_t> &partition)
{
    partition.clear();

    std::vector<size_t> stack;
    std::vector<bool> explored(mesh.face.size(),false);

    stack.push_back(IndexF);
    explored[IndexF]=true;
    do
    {
        size_t currF=stack.back();
        stack.pop_back();

        partition.push_back(currF);
        for (size_t j=0;j<3;j++)
        {
            if (mesh.face[currF].IsFaceEdgeS(j))continue;

            int NextFIndex=vcg::tri::Index(mesh,mesh.face[currF].FFp(j));

            if (explored[NextFIndex])continue;

            explored[NextFIndex]=true;
            stack.push_back(NextFIndex);
        }
    }while (!stack.empty());
}

template <class MeshType>
void RetrievePatchesFromSelEdges(MeshType &mesh,
                                 std::vector<size_t> &StartF,
                                 std::vector<std::vector<size_t> > &Partitions)
{
    Partitions.clear();
    vcg::tri::UpdateFlags<MeshType>::FaceClearV(mesh);
    for (size_t i=0;i<StartF.size();i++)
    {
        size_t IndexF=StartF[i];
        if (mesh.face[IndexF].IsV())continue;

        std::vector<size_t> partition;
        RetrievePatchFromSelEdges(mesh,IndexF,partition);

        for (size_t j=0;j<partition.size();j++)
            mesh.face[partition[j]].SetV();

        Partitions.push_back(partition);
    }
}

template <class MeshType>
void ColorMeshByPartitions(MeshType &mesh,const std::vector<std::vector<size_t> > &Partitions)
{
    for (size_t i=0;i<Partitions.size();i++)
    {
        vcg::Color4b CurrCol=vcg::Color4b::Scatter(Partitions.size(),i);
        for (size_t j=0;j<Partitions[i].size();j++)
            mesh.face[Partitions[i][j]].C()=CurrCol;
    }
}


std::vector<size_t> FilterFromSet(const std::vector<size_t> &ToFilter,
                                  const std::set<size_t> &FilterSet)
{
    std::vector<size_t> Filtered;
    for (size_t i=0;i<ToFilter.size();i++)
    {
        if(FilterSet.count(ToFilter[i])==0)continue;
        Filtered.push_back(ToFilter[i]);
    }
    return Filtered;
}

template <class MeshType>
bool SelectMeshPatchBorders(MeshType &mesh,
                            const std::vector<CandidateTrace> &TraceSet)
{
    //std::set<std::pair<size_t,size_t> > BorderPatches;
    vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(mesh);
    //first add borders
    for (size_t i=0;i<mesh.face.size();i++)
        for (size_t j=0;j<3;j++)
        {
            if (!mesh.face[i].IsB(j))continue;
            mesh.face[i].SetFaceEdgeS(j);
        }

    std::set<std::pair<size_t,size_t> > BorderEdges;
    for (size_t i=0;i<TraceSet.size();i++)
    {
        if (TraceSet[i].PathNodes.size()==0)continue;

        size_t Limit=TraceSet[i].PathNodes.size()-1;
        if (TraceSet[i].IsLoop)Limit++;
        for (size_t j=0;j<Limit;j++)
        {
            size_t IndexN0=TraceSet[i].PathNodes[j];
            size_t IndexN1=TraceSet[i].PathNodes[(j+1)%TraceSet[i].PathNodes.size()];
            size_t IndexV0=VertexFieldGraph<MeshType>::NodeVertI(IndexN0);
            size_t IndexV1=VertexFieldGraph<MeshType>::NodeVertI(IndexN1);
            BorderEdges.insert(std::pair<size_t,size_t>(std::min(IndexV0,IndexV1),
                                                        std::max(IndexV0,IndexV1)));
        }
    }

    std::vector<size_t> NumSel(mesh.vert.size(),0);
    for (size_t i=0;i<mesh.face.size();i++)
        for (size_t j=0;j<3;j++)
        {
            size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
            size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
            std::pair<size_t,size_t> Key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
            if (BorderEdges.count(Key)==0)continue;
            mesh.face[i].SetFaceEdgeS(j);
            NumSel[IndexV0]++;
            NumSel[IndexV1]++;
        }

    for (size_t i=0;i<NumSel.size();i++)
    {
        if (mesh.vert[i].IsB())continue;
        assert(NumSel[i]%2==0);
        if (NumSel[i]==2)return false;
    }
    return true;
}

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

//this is a function to remove some trace that is collinear, it should never happen
//but sometime it happens for a few cases, still should check the cause
template <class MeshType>
void RemoveCollinearTraces(std::vector<CandidateTrace> &TraceSet)
{
    
    std::vector<CandidateTrace> SwapTraceSet;
    std::set<std::pair<size_t,size_t> > PathEdges;
    //for each edge set the direction per vertex
    for (size_t i=0;i<TraceSet.size();i++)
    {
        if (TraceSet[i].PathNodes.size()==0)continue;
        size_t Limit=TraceSet[i].PathNodes.size()-1;
        if (TraceSet[i].IsLoop)
            Limit++;
        bool isOK=true;
        for (size_t j=0;j<Limit;j++)
        {
            size_t IndexN0=TraceSet[i].PathNodes[j];
            size_t IndexN1=TraceSet[i].PathNodes[(j+1)%TraceSet[i].PathNodes.size()];
            size_t IndexV0=VertexFieldGraph<MeshType>::NodeVertI(IndexN0);
            size_t IndexV1=VertexFieldGraph<MeshType>::NodeVertI(IndexN1);

            assert(IndexV0!=IndexV1);

            std::pair<size_t,size_t> EdgeKey(std::min(IndexV0,IndexV1),
                                             std::max(IndexV0,IndexV1));

            if (PathEdges.count(EdgeKey)>0)
            {
                std::cout<<"WARNING DOUBLE EDGE"<<std::endl;
                isOK=false;
                break;
            }

            PathEdges.insert(EdgeKey);
        }
        if (isOK)
            SwapTraceSet.push_back(TraceSet[i]);
    }
    TraceSet=SwapTraceSet;
}

template <class MeshType>
void GetEdgeDirVertMap(const VertexFieldGraph<MeshType> &VFGraph,
                       const std::vector<CandidateTrace> &TraceSet,
                       std::map<EdgeVert,size_t> &EdgeDirVert)
{
    EdgeDirVert.clear();

    //for each edge set the direction per vertex
    for (size_t i=0;i<TraceSet.size();i++)
    {
        if (TraceSet[i].PathNodes.size()==0)continue;
        size_t Limit=TraceSet[i].PathNodes.size()-1;
        if (TraceSet[i].IsLoop)
            Limit++;
        for (size_t j=0;j<Limit;j++)
        {
            size_t IndexN0=TraceSet[i].PathNodes[j];
            size_t IndexN1=TraceSet[i].PathNodes[(j+1)%TraceSet[i].PathNodes.size()];
            size_t IndexV0=VertexFieldGraph<MeshType>::NodeVertI(IndexN0);
            size_t IndexV1=VertexFieldGraph<MeshType>::NodeVertI(IndexN1);

            assert(IndexV0!=IndexV1);
            size_t MinV=std::min(IndexV0,IndexV1);
            size_t MaxV=std::max(IndexV0,IndexV1);

            size_t DirV0=VertexFieldGraph<MeshType>::NodeDirI(IndexN0);
            size_t DirV1=VertexFieldGraph<MeshType>::NodeDirI(IndexN1);
            EdgeVert EdgeKey0(MinV,MaxV,IndexV0);
            EdgeVert EdgeKey1(MinV,MaxV,IndexV1);

            if (EdgeDirVert.count(EdgeKey0)>0)
            {
                std::cout<<"WARNING DOUBLE EDGE"<<std::endl;
                MeshType traceMesh;
                std::vector<bool> Selected(TraceSet.size(),false);
                Selected[i]=true;
                MeshTraces(VFGraph,TraceSet,Selected,traceMesh);
                vcg::tri::io::ExporterPLY<MeshType>::Save(VFGraph.Mesh(),"double_direction_domain.ply",vcg::tri::io::Mask::IOM_FACECOLOR);
                vcg::tri::io::ExporterPLY<MeshType>::Save(traceMesh,"double_direction_error.ply",vcg::tri::io::Mask::IOM_FACECOLOR);
                assert(0);
            }

            if (EdgeDirVert.count(EdgeKey1)>0)
            {
                std::cout<<"WARNING DOUBLE EDGE"<<std::endl;
                MeshType traceMesh;
                std::vector<bool> Selected(TraceSet.size(),false);
                Selected[i]=true;
                MeshTraces(VFGraph,TraceSet,Selected,traceMesh);
                vcg::tri::io::ExporterPLY<MeshType>::Save(VFGraph.Mesh(),"double_direction_domain.ply",vcg::tri::io::Mask::IOM_FACECOLOR);
                vcg::tri::io::ExporterPLY<MeshType>::Save(traceMesh,"double_direction_error.ply",vcg::tri::io::Mask::IOM_FACECOLOR);
                assert(0);
            }
            assert(EdgeDirVert.count(EdgeKey0)==0);
            assert(EdgeDirVert.count(EdgeKey1)==0);

            EdgeDirVert[EdgeKey0]=DirV0;
            //            if (TraceSet[i].IsLoop)
            //                EdgeDirVert[EdgeKey0].push_back((DirV0+2)%4);
            //            else
            //                if (j>0)//add the inverse direction in case it is not starting
            //                    EdgeDirVert[EdgeKey0].push_back((DirV0+2)%4);

            EdgeDirVert[EdgeKey1]=((DirV1+2)%4);//put the inverse cause look internally the interval
            //            if (TraceSet[i].IsLoop)
            //                EdgeDirVert[EdgeKey1].push_back((DirV1+2)%4);
            //            else
            //                if (j!=(Limit-1))//add the inverse direction in case it is not ending
            //                    EdgeDirVert[EdgeKey1].push_back((DirV1+2)%4);
        }
    }

    //do the same for borders
    for (size_t i=0;i<VFGraph.Mesh().face.size();i++)
        for (size_t j=0;j<3;j++)
        {
            if (!vcg::face::IsBorder(VFGraph.Mesh().face[i],j))continue;
            size_t IndexV0=vcg::tri::Index(VFGraph.Mesh(),VFGraph.Mesh().face[i].V0(j));
            size_t IndexV1=vcg::tri::Index(VFGraph.Mesh(),VFGraph.Mesh().face[i].V1(j));
            size_t DirFlatV0,DirFlatV1;
            VertexFieldQuery<MeshType>::GetEdgeDir(VFGraph,IndexV0,IndexV1,DirFlatV0,DirFlatV1);

            assert(IndexV0!=IndexV1);
            size_t MinV=std::min(IndexV0,IndexV1);
            size_t MaxV=std::max(IndexV0,IndexV1);

            EdgeVert EdgeKey0(MinV,MaxV,IndexV0);
            EdgeVert EdgeKey1(MinV,MaxV,IndexV1);

            assert(EdgeDirVert.count(EdgeKey0)==0);
            EdgeDirVert[EdgeKey0]=DirFlatV0;
            assert(EdgeDirVert.count(EdgeKey1)==0);
            EdgeDirVert[EdgeKey1]= DirFlatV1;
        }
}

template <class MeshType>
void FindPerVertDirs(const VertexFieldGraph<MeshType> &VFGraph,
                     const std::vector<size_t> &Partition,
                     std::map<EdgeVert,size_t> &EdgeDirVert,
                     std::vector<std::vector<size_t> > &DirVert)
{
    typedef typename MeshType::CoordType CoordType;

    DirVert.clear();
    DirVert.resize(VFGraph.Mesh().vert.size());
    std::set<std::pair<CoordType,CoordType> > AddedEdges;
    for (size_t i=0;i<Partition.size();i++)
    {
        size_t IndexF=Partition[i];
        for (size_t e=0;e<3;e++)
        {
            size_t IndexV0=vcg::tri::Index(VFGraph.Mesh(),
                                           VFGraph.Mesh().face[IndexF].V0(e));

            size_t IndexV1=vcg::tri::Index(VFGraph.Mesh(),
                                           VFGraph.Mesh().face[IndexF].V1(e));

            CoordType Pos0=VFGraph.Mesh().vert[IndexV0].P();
            CoordType Pos1=VFGraph.Mesh().vert[IndexV1].P();
            std::pair<CoordType,CoordType> KeyE(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
            if (AddedEdges.count(KeyE)>0)continue;
            AddedEdges.insert(KeyE);


            size_t MinV=std::min(IndexV0,IndexV1);
            size_t MaxV=std::max(IndexV0,IndexV1);

            EdgeVert EdgeKey0(MinV,MaxV,IndexV0);
            if (EdgeDirVert.count(EdgeKey0)>0)
            {
                size_t EdgeDir=EdgeDirVert[EdgeKey0];
                DirVert[IndexV0].push_back(EdgeDir);
            }


            EdgeVert EdgeKey1(MinV,MaxV,IndexV1);
            if (EdgeDirVert.count(EdgeKey1)>0)
            {
                size_t EdgeDir=EdgeDirVert[EdgeKey1];
                DirVert[IndexV1].push_back(EdgeDir);
            }
        }
    }
}

bool HasOrthogonalCross(const std::vector<size_t> &Directions)
{
    if (Directions.size()<2)return false;
    for (size_t i=0;i<Directions.size()-1;i++)
        for (size_t j=(i+1);j<Directions.size();j++)
        {
            size_t Dir0=Directions[i];
            size_t Dir1=Directions[j];
            if ((Dir0 % 2)!=(Dir1 % 2))return true;
        }
    return false;
}

bool HasNarrowCross(const std::vector<size_t> &Directions)
{
    if (Directions.size()<2)return false;
    for (size_t i=0;i<Directions.size()-1;i++)
        for (size_t j=(i+1);j<Directions.size();j++)
        {
            size_t Dir0=Directions[i];
            size_t Dir1=Directions[j];
            if (Dir0 == Dir1)return true;
        }
    return false;
}

template <class MeshType>
void FindPartitionsCorners(const VertexFieldGraph<MeshType> &VFGraph,
                           const std::vector<TypeVert> &VertType,
                           const std::vector<CandidateTrace> &TraceSet,
                           const std::vector<std::vector<size_t> > &Partitions,
                           std::vector<std::vector<size_t> > &PartitionCorners)
{
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::ScalarType ScalarType;

    PartitionCorners.clear();
    PartitionCorners.resize(Partitions.size());

    //for each edge set the direction per vertex
    std::map<EdgeVert,size_t> EdgeDirVert;
    GetEdgeDirVertMap(VFGraph,TraceSet,EdgeDirVert);

    //    //then find per partition per edge angle
    //    std::map<std::pair<size_t,size_t>,ScalarType> PartitionVertAngle;
    //    FindTraceAngles(VFGraph,Partitions,PartitionVertAngle);

    //first initialize the quality of each face with the partition
    for (size_t i=0;i<Partitions.size();i++)
    {
        //then add the corners
        for (size_t j=0;j<Partitions[i].size();j++)
        {
            size_t IndexF=Partitions[i][j];
            VFGraph.Mesh().face[IndexF].Q()=i;
        }
    }

    //then store per vertex
    //then go over all partitions
    for (size_t i=0;i<Partitions.size();i++)
    {
        //map for each vertex of the partitions the number of directions
        //that have been sampled
        std::vector<std::vector<size_t> > DirVert;
        FindPerVertDirs(VFGraph,Partitions[i],EdgeDirVert,DirVert);

        //then add the corners
        for (size_t j=0;j<Partitions[i].size();j++)
        {
            size_t IndexF=Partitions[i][j];
            for (size_t e=0;e<3;e++)
            {
                size_t IndexV=vcg::tri::Index(VFGraph.Mesh(),VFGraph.Mesh().face[IndexF].V(e));
                if (VertType[IndexV]==TVConvex)//this is convex then is a corner for sure
                    PartitionCorners[i].push_back(IndexV);
                else
                {
                    bool PossibleCorner=(HasOrthogonalCross(DirVert[IndexV])||
                                         HasNarrowCross(DirVert[IndexV]));
                    if (!PossibleCorner)continue;

                    //check incomplete concave
                    if (VertType[IndexV]==TVConcave)
                    {
                        std::vector<FaceType*> StarF;
                        std::vector<int> LocalIndexV;
                        vcg::face::VFStarVF(&VFGraph.Mesh().vert[IndexV],StarF,LocalIndexV);
                        ScalarType Angle=0;
                        for (size_t j=0;j<StarF.size();j++)
                        {
                            FaceType* f=StarF[j];
                            int currV=LocalIndexV[j];
                            if (f->Q()!=i)continue;
                            Angle+=vcg::face::WedgeAngleRad(*f,currV);
                        }
                        if (Angle>(M_PI))continue;
                    }
                    PartitionCorners[i].push_back(IndexV);
                }
            }
        }
        std::sort(PartitionCorners[i].begin(),PartitionCorners[i].end());
        std::vector<size_t>::iterator it;
        it = std::unique (PartitionCorners[i].begin(),PartitionCorners[i].end());
        PartitionCorners[i].resize( std::distance(PartitionCorners[i].begin(),it) );
    }


    //std::cout<<"3"<<std::endl;
}

bool MergeContiguousPathStep(std::vector<CandidateTrace> &TraceSet)
{
    for (size_t i=0;i<TraceSet.size();i++)
        for (size_t j=0;j<TraceSet.size();j++)
        {
            if (i==j)continue;
            if (TraceSet[i].IsLoop)continue;
            if (TraceSet[j].IsLoop)continue;
            if(TraceSet[i].PathNodes.size()==0)continue;
            if(TraceSet[j].PathNodes.size()==0)continue;
            size_t Node0_0=TraceSet[i].PathNodes[0];
            size_t Node0_1=TraceSet[i].PathNodes.back();
            size_t Node1_0=TraceSet[j].PathNodes[0];
            size_t Node1_1=TraceSet[j].PathNodes.back();
            if (Node0_1==Node1_0)
            {
                //remove last one
                TraceSet[i].PathNodes.pop_back();
                //add the rest
                TraceSet[i].PathNodes.insert(TraceSet[i].PathNodes.end(),
                                             TraceSet[j].PathNodes.begin(),
                                             TraceSet[j].PathNodes.end());
                TraceSet[i].ToType=TraceSet[j].ToType;

                TraceSet[j].PathNodes.clear();
                return true;
            }
            if (Node1_1==Node0_0)
            {
                //remove last one
                TraceSet[j].PathNodes.pop_back();
                //add the rest
                TraceSet[j].PathNodes.insert(TraceSet[j].PathNodes.end(),
                                             TraceSet[i].PathNodes.begin(),
                                             TraceSet[i].PathNodes.end());

                TraceSet[j].ToType=TraceSet[i].ToType;

                TraceSet[i].PathNodes.clear();
                return true;
            }
        }
    return false;
}

void MergeContiguousPaths(std::vector<CandidateTrace> &TraceSet)
{
    while(MergeContiguousPathStep(TraceSet)){};
}



void GetIsLoop(const std::vector<CandidateTrace> &TraceSet,
               std::vector<bool> &CurrChosenIsLoop)
{
    CurrChosenIsLoop.clear();
    for (size_t i=0;i<TraceSet.size();i++)
        CurrChosenIsLoop.push_back(TraceSet[i].IsLoop);
}

template <class MeshType>
bool SplitSomeNonOKPartition(const VertexFieldGraph<MeshType> &VFGraph,
                             const CandidateTrace &TestTrace,
                             const std::vector<std::vector<size_t> > &Partitions,
                             const std::vector<int > &FacePartitions,
                             const std::vector<PatchType> &PartitionType)
{
    typedef typename MeshType::FaceType FaceType;

    bool IsCurrLoop=TestTrace.IsLoop;
    size_t TestL=TestTrace.PathNodes.size()-1;
    if (IsCurrLoop)TestL++;
    //        std::cout<<"There are "<<Partitions.size()<<" partitions "<<std::endl;
    //        if (PartitionType[0]==IsOK)
    //            std::cout<<"Is OK "<<std::endl;

    for (size_t i=0;i<TestL;i++)
    {
        size_t IndexN0=TestTrace.PathNodes[i];
        size_t IndexN1=TestTrace.PathNodes[(i+1)%TestTrace.PathNodes.size()];
        if (VFGraph.AreTwin(IndexN0,IndexN1))continue;

        vcg::face::Pos<FaceType> Pos=VFGraph.GetNodesPos(IndexN0,IndexN1);

        //border, cannot split any partition
        if (Pos.IsBorder())continue;
        FaceType *f0=Pos.F();
        FaceType *f1=Pos.FFlip();
        assert(f0!=f1);
        size_t IndexF0=vcg::tri::Index(VFGraph.Mesh(),f0);
        size_t IndexF1=vcg::tri::Index(VFGraph.Mesh(),f1);
        assert((IndexF0>=0)&&(IndexF0<FacePartitions.size()));
        assert((IndexF1>=0)&&(IndexF1<FacePartitions.size()));
        size_t Partition0=FacePartitions[IndexF0];
        size_t Partition1=FacePartitions[IndexF1];

        assert(Partition0>=0);
        assert(Partition0<Partitions.size());
        assert(Partition0<PartitionType.size());

        assert(Partition1>=0);
        assert(Partition1<Partitions.size());
        assert(Partition1<PartitionType.size());

        //passing across two partitions
        if (Partition0!=Partition1)continue;

        if (PartitionType[Partition0]!=IsOK)return true;
    }
    return false;
}

template <class MeshType>
void ColorMeshByExpValence(MeshType &mesh,
                           const std::vector<std::vector<size_t> > &PartitionCorners,
                           const std::vector<std::vector<size_t> > &Partitions)
{
    //vcg::tri::UpdateSelection<MeshType>::FaceClear(Mesh());
    for (size_t i=0;i<Partitions.size();i++)
    {
        int ExpVal=ExpectedValence(mesh,Partitions[i]);
        vcg::Color4b CurrCol;
        if (ExpVal==-1)
            CurrCol=vcg::Color4b::Gray;
        else
        {
            if (PartitionCorners[i].size()==ExpVal)
                CurrCol=vcg::Color4b::Green;
            if (PartitionCorners[i].size()!=ExpVal)
                CurrCol=vcg::Color4b::Red;
        }
        for (size_t j=0;j<Partitions[i].size();j++)
            mesh.face[Partitions[i][j]].C()=CurrCol;
    }
}

template <class MeshType>
void ColorMeshByValence(MeshType &mesh,
                        const std::vector<std::vector<size_t> > &PartitionCorners,
                        const std::vector<std::vector<size_t> > &Partitions,
                        size_t MinVal,size_t MaxVal)
{
    //vcg::tri::UpdateSelection<MeshType>::FaceClear(Mesh());
    for (size_t i=0;i<PartitionCorners.size();i++)
    {
        vcg::Color4b CurrCol=vcg::Color4b::Gray;

        if (PartitionCorners[i].size()==MinVal)
            CurrCol=vcg::Color4b::Yellow;
        if (PartitionCorners[i].size()==5)
            CurrCol=vcg::Color4b::Blue;
        if (PartitionCorners[i].size()==MaxVal)
            CurrCol=vcg::Color4b::Cyan;
        if ((PartitionCorners[i].size()<MinVal)||
                (PartitionCorners[i].size()>MaxVal))
            CurrCol=vcg::Color4b::Red;


        for (size_t j=0;j<Partitions[i].size();j++)
            mesh.face[Partitions[i][j]].C()=CurrCol;
    }
}

template <class MeshType>
void ColorMeshByTopology(MeshType &mesh,
                         const std::vector<std::vector<size_t> > &Partitions,
                         const std::vector<PatchType> &PartitionType)
{
    for (size_t i=0;i<Partitions.size();i++)
    {
        vcg::Color4b CurrCol=vcg::Color4b::Gray;

        switch (PartitionType[i])
        {
        case LowCorners:CurrCol=vcg::Color4b::Blue;break;
        case HighCorners:CurrCol=vcg::Color4b::Magenta;break;
        case NonDisk:CurrCol=vcg::Color4b::Red;break;
        case HasEmitter:CurrCol=vcg::Color4b::Yellow;break;
        case MaxCClarkability:CurrCol=vcg::Color4b::LightGreen;break;
        case NonMatchValence:CurrCol=vcg::Color4b::Magenta;break;
        default: CurrCol=vcg::Color4b::Gray;break;
        }

        for (size_t j=0;j<Partitions[i].size();j++)
            mesh.face[Partitions[i][j]].C()=CurrCol;
    }
}

template <class MeshType>
bool TraceDirectPath(VertexFieldGraph<MeshType> &VFGraph,
                     const size_t StartingNode,
                     std::vector<size_t> &PathN)
{
    //deselect source seletion in case is also a receiver
    bool restoreSel=false;
    if (VFGraph.IsSelected(StartingNode))
    {
        VFGraph.DeSelect(StartingNode);
        restoreSel=true;
    }
    bool Traced=VertexFieldQuery<MeshType>::TraceToSelected(VFGraph,StartingNode,PathN);
    if (restoreSel)
        VFGraph.Select(StartingNode);

    if (!Traced)return false;

    bool SelfInt=VertexFieldQuery<MeshType>::SelfIntersect(VFGraph,PathN,false);
    if (SelfInt)return false;

    return true;
}

template <class MeshType>
bool TraceDijkstraPath(VertexFieldGraph<MeshType> &VFGraph,
                       const size_t StartingNode,
                       const typename MeshType::ScalarType &Drift,
                       const typename MeshType::ScalarType &MaxWeight,
                       std::vector<size_t> &PathN)
{
    typename VertexFieldQuery<MeshType>::ShortParam SParam;
    SParam.StartNode.push_back(StartingNode);
    SParam.MaxTwin=MAX_TWIN_DIKSTRA;
    SParam.MaxWeight=MaxWeight;
    SParam.OnlyDirect=false;
    SParam.DriftPenalty=Drift;

    //deselect source selection in case is also a receiver
    bool restoreSel=false;
    if (VFGraph.IsSelected(StartingNode))
    {
        VFGraph.DeSelect(StartingNode);
        restoreSel=true;
    }
    bool Traced=VertexFieldQuery<MeshType>::ShortestPath(VFGraph,SParam,PathN);
    if (restoreSel)
        VFGraph.Select(StartingNode);
    if (!Traced)return false;
    bool SelfInt=VertexFieldQuery<MeshType>::SelfIntersect(VFGraph,PathN,false);
    if (SelfInt)return false;
    return true;
}

template <class MeshType>
bool TraceLoopPath(VertexFieldGraph<MeshType> &VFGraph,
                   const size_t StartingNode,
                   const typename MeshType::ScalarType &Drift,
                   std::vector<size_t> &PathN)
{
    bool ClosedLoop=VertexFieldQuery<MeshType>::FindLoop(VFGraph,StartingNode,PathN,Drift);
    if (!ClosedLoop)return false;
    bool SelfInt=VertexFieldQuery<MeshType>::SelfIntersect(VFGraph,PathN,true);
    if (SelfInt)return false;
    return true;
}

template <class MeshType>
void GetFarthestCoupleVertx(MeshType &mesh,
                            int &MinIdx,
                            int &MaxIdx)
{
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;

    ScalarType maxD=0;
    for (size_t i=0;i<mesh.vert.size()-1;i++)
        for (size_t j=(i+1);j<mesh.vert.size();j++)
        {
            CoordType Pos0=mesh.vert[i].P();
            CoordType Pos1=mesh.vert[j].P();
            ScalarType testD=(Pos0-Pos1).Norm();
            if (testD<=maxD)continue;
            maxD=testD;
            MinIdx=i;
            MaxIdx=j;
        }
}


template <class MeshType>
class PatchTracer
{
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;

    VertexFieldGraph<MeshType> &VFGraph;

public:
    //vertices types and ortho directions, used to select sources
    ScalarType Drift;

    bool split_on_removal;
    bool DebugMsg;
    //bool avoid_increase_valence;
    //bool avoid_collapse_irregular;
    bool away_from_singular;
    //    ScalarType max_lenght_distortion;
    //    ScalarType max_lenght_variance;
    ScalarType sample_ratio;
    ScalarType CClarkability;
    bool match_valence;
    //ScalarType max_patch_area;
    size_t MinVal;
    size_t MaxVal;
    size_t Concave_Need;

    //bool TraceLoopsBorders;

    std::vector<std::vector<size_t> > Partitions;
    std::vector<int > FacePartitions;
    std::vector<PatchType> PartitionType;
    std::vector<std::vector<size_t> > PartitionCorners;

    std::vector<PatchInfo<ScalarType> > PatchInfos;

private:

    std::vector<TypeVert> VertType;
    std::vector<TypeVert > NodeEmitterTypes;
    std::vector<TypeVert> NodeReceiverTypes;
    ScalarType MaxNarrowWeight;
    bool AllReceivers;
    std::vector<ScalarType> CurrNodeDist;

    //    std::vector<std::pair<ScalarType,size_t> > CandidatesPathLenghts;
    //    std::vector<std::pair<ScalarType,size_t> > CandidatesPathDistances;


    //DATA STRUCTURES FOR THE PATH AND THE NEEDS FOR EACH VERTEX
    std::vector<CandidateTrace> Candidates;
public:
    std::vector<CandidateTrace> ChoosenPaths;
private:
    std::vector<size_t> VerticesNeeds;

    //std::vector<std::vector<ScalarType> > EdgeL;
    std::map<std::pair<size_t,size_t>,ScalarType> EdgeL;
    ScalarType avgEdge;

    void InitAvEdge()
    {
        avgEdge=0;
        size_t Num=0;
        for (size_t i=0;i<Mesh().face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                avgEdge+=(Mesh().face[i].cP0(j)-
                          Mesh().face[i].cP1(j)).Norm();
                Num++;
            }
        avgEdge/=Num;
    }

    void InitEdgeLBorders()
    {
        for (size_t i=0;i<Mesh().face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!vcg::face::IsBorder(Mesh().face[i],j))continue;
                CoordType P0=Mesh().face[i].P0(j);
                CoordType P1=Mesh().face[i].P1(j);
                size_t IndexV0=vcg::tri::Index(Mesh(),Mesh().face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(Mesh(),Mesh().face[i].V1(j));
                std::pair<size_t,size_t> Key(std::min(IndexV0,IndexV1),
                                             std::max(IndexV0,IndexV1));
                assert(EdgeL.count(Key)==0);
                EdgeL[Key]=(P1-P0).Norm();
            }
    }

    ScalarType FieldL(const size_t &IndexNode0,
                      const size_t &IndexNode1)
    {
        CoordType Dir0=VFGraph.NodeDir(IndexNode0);
        CoordType Dir1=VFGraph.NodeDir(IndexNode1);
        CoordType AvgDir=(Dir0+Dir1);
        AvgDir.Normalize();
        size_t indexV0=VertexFieldGraph<MeshType>::NodeVertI(IndexNode0);
        size_t indexV1=VertexFieldGraph<MeshType>::NodeVertI(IndexNode1);
        CoordType Pos0=Mesh().vert[indexV0].RPos;//NodePos(IndexNode0);
        CoordType Pos1=Mesh().vert[indexV1].RPos;
        ScalarType L=(Pos1-Pos0)*AvgDir;
        return L;
        //        assert(IndexNode<NodeJumps.size());
        //        return (NodeJumps[IndexNode]);
    }

    void AddEdgeL(const CandidateTrace &CurrCand)
    {
        if(CurrCand.PathNodes.size()==0)return;
        assert(CurrCand.PathNodes.size()>=2);
        size_t Limit=CurrCand.PathNodes.size()-1;
        size_t Size=CurrCand.PathNodes.size();
        if (CurrCand.IsLoop)Limit++;
        for (size_t j=0;j<Limit;j++)
        {
            size_t N0=CurrCand.PathNodes[j];
            size_t N1=CurrCand.PathNodes[(j+1)%Size];
            size_t IndexV0=VertexFieldGraph<MeshType>::NodeVertI(N0);
            size_t IndexV1=VertexFieldGraph<MeshType>::NodeVertI(N1);
            std::pair<size_t,size_t> Key(std::min(IndexV0,IndexV1),
                                         std::max(IndexV0,IndexV1));
            assert(EdgeL.count(Key)==0);
            ScalarType currL=FieldL(N0,N1);
            assert(EdgeL.count(Key)==0);
            EdgeL[Key]=currL;
        }
    }

    void InitEdgeL()
    {
        EdgeL.clear();
        InitEdgeLBorders();

        for (size_t i=0;i<ChoosenPaths.size();i++)
            AddEdgeL(ChoosenPaths[i]);
    }

    //CLASSIFICATION OF DIFFERENT KINDS OF VERTICES
    void InitVertType(const std::vector<size_t> &ConvexV,
                      const std::vector<size_t> &ConcaveV,
                      const std::vector<size_t> &NarrowV,
                      const std::vector<size_t> &FlatV)
    {

        VertType.clear();
        VertType.resize(Mesh().vert.size(),TVInternal);

        for (size_t i=0;i<ConvexV.size();i++)
        {
            size_t IndexV=ConvexV[i];
            VertType[IndexV]=TVConvex;
        }

        for (size_t i=0;i<ConcaveV.size();i++)
        {
            size_t IndexV=ConcaveV[i];
            VertType[IndexV]=TVConcave;
        }

        for (size_t i=0;i<NarrowV.size();i++)
        {
            size_t IndexV=NarrowV[i];
            VertType[IndexV]=TVNarrow;
        }

        for (size_t i=0;i<FlatV.size();i++)
        {
            size_t IndexV=FlatV[i];
            VertType[IndexV]=TVFlat;
        }
    }

    void InitVertType()
    {
        std::vector<size_t> ConvexV,ConcaveV,FlatV,NarrowV;
        VertexClassifier<MeshType>::FindConvexV(VFGraph,ConvexV);
        VertexClassifier<MeshType>::FindConcaveV(VFGraph,ConcaveV);
        VertexClassifier<MeshType>::FindNarrowV(VFGraph,NarrowV);
        VertexClassifier<MeshType>::FindFlatV(VFGraph,ConvexV,ConcaveV,NarrowV,FlatV);
        InitVertType(ConvexV,ConcaveV,NarrowV,FlatV);
    }

    //THIS INITIALIZE THE EMITTERS FOR EACH VERTEX
    void InitEmitters()
    {
        NodeEmitterTypes.clear();
        NodeReceiverTypes.clear();
        NodeEmitterTypes.resize(VFGraph.NumNodes(),TVNone);
        NodeReceiverTypes.resize(VFGraph.NumNodes(),TVNone);


        std::vector<std::vector<CoordType> > VertFlatDir,VertOrthoDir;
        VertexEmitter<MeshType>::GetOrthoFlatDirections(VFGraph,VertFlatDir,VertOrthoDir);

        for (size_t i=0;i<Mesh().vert.size();i++)
        {

            if (VertType[i]==TVConvex)//convex ones, no emitters
                continue;

            if (VertType[i]==TVInternal)//this will be with sampling
                continue;

            if (VertType[i]==TVFlat)
            {
                size_t Emitter,Receiver;
                VertexEmitter<MeshType>::ComputeFlatEmitterReceivers(VFGraph,VertOrthoDir,VertFlatDir,i,Emitter,Receiver);
                assert(Emitter!=Receiver);
                assert(NodeEmitterTypes[Emitter]==TVNone);
                NodeEmitterTypes[Emitter]=TVFlat;
                assert(NodeReceiverTypes[Receiver]==TVNone);
                NodeReceiverTypes[Receiver]=TVFlat;
            }
            if (VertType[i]==TVConcave)
            {
                std::vector<size_t> Emitter,Receiver;
                VertexEmitter<MeshType>::ComputeConcaveEmitterReceivers(VFGraph,VertOrthoDir,VertFlatDir,i,Emitter,Receiver);
                for (size_t i=0;i<Emitter.size();i++)
                {
                    assert(NodeEmitterTypes[Emitter[i]]==TVNone);
                    NodeEmitterTypes[Emitter[i]]=TVConcave;
                }
                for (size_t i=0;i<Receiver.size();i++)
                {
                    NodeReceiverTypes[Receiver[i]]=TVConcave;
                }
            }

            if (VertType[i]==TVNarrow)
            {
                size_t Emitter,Receiver;
                VertexEmitter<MeshType>::ComputeNarrowEmitterReceivers(VFGraph,VertOrthoDir,VertFlatDir,i,Emitter,Receiver);
                assert(Emitter!=Receiver);
                assert(NodeEmitterTypes[Emitter]==TVNone);
                NodeEmitterTypes[Emitter]=TVNarrow;
                assert(NodeReceiverTypes[Receiver]==TVNone);
                NodeReceiverTypes[Receiver]=TVNarrow;
            }
        }

        //then add internal one for loops and other tracing
        std::vector<size_t> StartingNodes;
        if(sample_ratio<1)
        {
            size_t sampleNum=Mesh().vert.size()*sample_ratio;//floor(sqrt(Mesh().vert.size())+0.5)*10*sample_ratio;
            sampleNum=std::max(sampleNum,(size_t)50);
            //SampleStartingNodes(false,sampleNum,StartingNodes);
            VertexFieldQuery<MeshType>::SamplePoissonNodes(VFGraph,sampleNum,StartingNodes);
        }
        else
        {
            for (size_t i=0;i<Mesh().vert.size();i++)
            {
                std::vector<size_t> IndexN;
                VertexFieldGraph<MeshType>::IndexNodes(i,IndexN);
                if(VFGraph.IsActive(IndexN[0]))
                    StartingNodes.push_back(IndexN[0]);
                if(VFGraph.IsActive(IndexN[1]))
                    StartingNodes.push_back(IndexN[1]);
            }
        }
        for (size_t i=0;i<StartingNodes.size();i++)
        {
            size_t IndexV=VertexFieldGraph<MeshType>::NodeVertI(StartingNodes[i]);
            //std::cout<<"Sampled V "<<IndexV<<std::endl;
            if (VertType[IndexV]!=TVInternal)continue;
            assert(NodeEmitterTypes[StartingNodes[i]]==TVNone);
            NodeEmitterTypes[StartingNodes[i]]=TVInternal;
            //Emitter[IndexV].push_back(StartingNodes[i]);
        }
    }

    //SET AS NON VALID ALL THE TANGENT NODES ON THE SIDE AND THE OTHER THAT CANNOT BE
    //USED FOR TRACING
    void InvalidateTangentNodes()
    {
        VFGraph.SetAllActive();
        //TangentConstraints.clear();
        for (size_t i=0;i<VertType.size();i++)
        {
            if (VertType[i]==TVConvex)
            {
                std::vector<size_t> NodesI;
                VertexFieldGraph<MeshType>::IndexNodes(i,NodesI);
                for (size_t j=0;j<NodesI.size();j++)
                    VFGraph.SetActive(NodesI[j],false);
            }

            if (VertType[i]==TVFlat)//||(VertType[i]==Narrow))
            {
                //get the nodes of a given vertex
                std::vector<size_t> NodesI;
                VertexFieldGraph<MeshType>::IndexNodes(i,NodesI);
                for (size_t j=0;j<NodesI.size();j++)
                {
                    if (NodeEmitterTypes[NodesI[j]]!=TVNone)continue;
                    if (NodeReceiverTypes[NodesI[j]]!=TVNone)continue;
                    VFGraph.SetActive(NodesI[j],false);
                }
            }

            if (VertType[i]==TVConcave)
            {
                //get the nodes of a given vertex
                std::vector<size_t> NodesI;
                VertexFieldGraph<MeshType>::IndexNodes(i,NodesI);
                for (size_t j=0;j<NodesI.size();j++)
                {
                    if (NodeEmitterTypes[NodesI[j]]!=TVNone)continue;
                    if (NodeReceiverTypes[NodesI[j]]!=TVNone)continue;
                    VFGraph.SetActive(NodesI[j],false);
                }
            }
        }
    }


    //EXPANDING CURRENT CANDIDATE PATHS
    void ExpandCandidates()
    {
        std::vector<CandidateTrace> ExpandedCandidates;

        size_t SelfIntN=0;
        for (size_t i=0;i<Candidates.size();i++)
        {
            //            if (Candidates[i].IsLoop)
            //                std::cout<<"IS LOOP"<<std::endl;
            //            else
            //                std::cout<<"NO LOOP"<<std::endl;
            //bool expanded=ExpandPath(Candidates[i].PathNodes,Candidates[i].IsLoop);
            bool expanded=VertexFieldQuery<MeshType>::ExpandPath(VFGraph,Candidates[i].PathNodes,
                                                                 Candidates[i].IsLoop,Drift);
            //std::cout<<"SIZE "<<Candidates[i].PathNodes.size()<<std::endl;
            bool SelfInt=VertexFieldQuery<MeshType>::SelfIntersect(VFGraph,Candidates[i].PathNodes,Candidates[i].IsLoop);
            if (SelfInt){SelfIntN++;continue;}
            if (expanded)
                ExpandedCandidates.push_back(Candidates[i]);
        }
        if (DebugMsg)
            std::cout<<"Self Intersections "<<SelfIntN<<std::endl;

        //all expanded do nothing (already chenged in place)
        if (ExpandedCandidates.size()==Candidates.size())return;
        Candidates=ExpandedCandidates;
    }

    int NumEmitters(size_t IndexV)
    {
        std::vector<size_t> IndexN;
        VertexFieldGraph<MeshType>::IndexNodes(IndexV,IndexN);
        int NumEmit=0;
        for (size_t i=0;i<IndexN.size();i++)
            if (NodeEmitterTypes[IndexN[i]]!=TVNone)NumEmit++;

        return NumEmit;
    }
    //INITILIAZE THE NEED FOR EACH NEEDED KIND OF VERT
    void InitVerticesNeeds()
    {
        VerticesNeeds.resize(Mesh().vert.size(),0);
        for (size_t i=0;i<VertType.size();i++)
            if (VertType[i]==TVNarrow)VerticesNeeds[i]=NARROW_NEED;
        for (size_t i=0;i<VertType.size();i++)
            if (VertType[i]==TVConcave)VerticesNeeds[i]=std::max(NumEmitters(i)-1,(int)Concave_Need);
    }

    //INITIALIZE THE STRUCTURES
    void InitStructures()
    {
        //InitOrthoDirections();

        InitVertType();

        InitEmitters();

        InvalidateTangentNodes();

        InitVerticesNeeds();

    }

    //SIMPLY REMOVE THE NON TRACED PATHS FROM THE CANDIDATES
    void CleanNonTracedCandidates()
    {
        std::vector<CandidateTrace> SwapCandidates;
        for (size_t i=0;i<Candidates.size();i++)
        {
            if ((Candidates[i].Updated) &&
                    (Candidates[i].PathNodes.size()==0))continue;
            SwapCandidates.push_back(Candidates[i]);
        }
        Candidates=SwapCandidates;
    }

    //UPDATE CURRENT CANDIDATES
    bool UpdateCandidates(const std::vector<bool> &CanReceive)
    {
        if (Candidates.size()==0)return false;
        VFGraph.Select(CanReceive);
        //        assert(CanReceive.size()==VFGraph.NumNodes());
        //        VFGraph.ClearSelection();
        //        for (size_t i=0;i<CanReceive.size();i++)
        //        {
        //            if (CanReceive[i])
        //                VFGraph.Select(i);
        //        }

        for (size_t i=0;i<Candidates.size();i++)
        {
            if (Candidates[i].Updated)continue;
            UpdateCandidate(VFGraph,Candidates[i],Drift,MaxNarrowWeight);
        }
        //finally erase the non traced ones
        CleanNonTracedCandidates();
        return (Candidates.size()>0);
    }

    void SortCandidatesByLenghts()
    {
        for (size_t i=0;i<Candidates.size();i++)
        {
            ScalarType currL=VertexFieldQuery<MeshType>::TraceLenght(VFGraph,Candidates[i].PathNodes,
                                                                     Candidates[i].IsLoop);
            Candidates[i].Priority=currL;
        }
        std::sort(Candidates.begin(),Candidates.end());
    }

    void SortCandidatesByDistances()
    {
        for (size_t i=0;i<Candidates.size();i++)
        {
            ScalarType currD=VertexFieldQuery<MeshType>::TraceAVGDistance(VFGraph,Candidates[i].PathNodes);
            Candidates[i].Priority=currD;
            //CandidatesPathDistances.push_back(std::pair<ScalarType,size_t>(currD,i));
        }
        std::sort(Candidates.begin(),Candidates.end());
        //should get the one with higher distance
        std::reverse(Candidates.begin(),Candidates.end());
    }

    void AddChoosen(const CandidateTrace &CurrC)
    {
        ChoosenPaths.push_back(CurrC);
        assert(CurrC.PathNodes.size()>=2);
        UpdateVertNeeds(CurrC.PathNodes);
        AddEdgeL(CurrC);
    }

    void ChooseGreedyByLength(bool UseVertNeeds=true,
                              bool UsePartitionNeeds=false)
    {
        //InitCandidatesPathLenghts();
        SortCandidatesByLenghts();

        //        size_t StartConflPath=0;
        //        if (checkOnlylastConfl)
        //            StartConflPath=ChoosenPaths.size();

        for (size_t i=0;i<Candidates.size();i++)
        {
            //get the current trace from the sorted ones
            //size_t CurrTraceIndex=CandidatesPathLenghts[i].second;
            std::vector<size_t > CurrTrace=Candidates[i].PathNodes;
            //bool IsCurrLoop=Candidates[CurrTraceIndex].IsLoop;

            //get the first vertex to check if it has been already traced or not
            size_t IndexN0=CurrTrace[0];
            size_t IndexN1=CurrTrace.back();

            size_t IndexV0=VertexFieldGraph<MeshType>::NodeVertI(IndexN0);
            size_t IndexV1=VertexFieldGraph<MeshType>::NodeVertI(IndexN1);

            //if it has already a trace then go on
            if ((UseVertNeeds)&&((VerticesNeeds[IndexV0]==0)&&(VerticesNeeds[IndexV1]==0)))continue;


            //bool collide=CollideWithChoosen(CurrTrace,IsCurrLoop,StartConflPath);
            bool collide = CollideWithCandidateSet(VFGraph,Candidates[i],ChoosenPaths);
            if (collide)continue;

            if ((UsePartitionNeeds)&&
                    (!SplitSomeNonOKPartition(VFGraph,Candidates[i],Partitions,FacePartitions,PartitionType)))
                continue;

            //ChoosenPaths.push_back(Candidates[i]);
            //assert(ChoosenPaths.back().PathNodes.size()>=2);
            //UpdateVertNeeds(CurrTrace);
            AddChoosen(Candidates[i]);
        }
    }

    void DeleteCandidates(const std::vector<bool> &To_Delete)
    {
        assert(To_Delete.size()==Candidates.size());
        std::vector<CandidateTrace> CandidatesSwap;
        for (size_t i=0;i<To_Delete.size();i++)
        {
            if (To_Delete[i])continue;
            CandidatesSwap.push_back(Candidates[i]);
        }
        Candidates=CandidatesSwap;
    }

    void UpdateVertNeeds(const std::vector<size_t> &TestTrace)
    {
        for (size_t i=0;i<TestTrace.size();i++)
        {
            size_t IndexN=TestTrace[i];
            size_t IndexV=VertexFieldGraph<MeshType>::NodeVertI(IndexN);

            //            TypeVert VType=VertType[IndexV];
            //            TypeVert NodeEType=NodeEmitterTypes[IndexN];
            //            TypeVert NodeRType=NodeReceiverTypes[IndexN];
            //            if ((VType==NodeEType)||(VType==NodeRType))
            if(VerticesNeeds[IndexV]>0)
                VerticesNeeds[IndexV]--;
        }
    }

    void UpdateVertNeedsFromChoosen()
    {
        for (size_t i=0;i<ChoosenPaths.size();i++)
            UpdateVertNeeds(ChoosenPaths[i].PathNodes);
    }

    bool SolveVertexNeed(const CandidateTrace &TestTrace)
    {
        bool IsCurrLoop=TestTrace.IsLoop;
        assert(!IsCurrLoop);

        size_t IndexN0=TestTrace.PathNodes[0];
        size_t IndexN1=TestTrace.PathNodes.back();

        size_t IndexV0=VertexFieldGraph<MeshType>::NodeVertI(IndexN0);
        size_t IndexV1=VertexFieldGraph<MeshType>::NodeVertI(IndexN1);
        if((VerticesNeeds[IndexV0]==0)&&(VerticesNeeds[IndexV1]==0))
            return false;

        return true;
    }

    bool ChooseNextByDistance(bool UseVertNeeds,
                              bool UsePartitionNeeds)
    {
        std::vector<bool> To_Delete(Candidates.size(),false);
        SortCandidatesByDistances();
        for (size_t i=0;i<Candidates.size();i++)
        {
            //std::vector<size_t > CurrTrace=Candidates[i].PathNodes;

            bool IsCurrLoop=Candidates[i].IsLoop;

            if (UseVertNeeds)
            {
                assert(!IsCurrLoop);
                if (!SolveVertexNeed(Candidates[i]))
                {
                    To_Delete[i]=true;
                    //std::cout<<"Not Needed by Needs "<<std::endl;
                    continue;
                }
            }
            //check if collide
            //bool collide=CollideWithChoosen(CurrTrace,IsCurrLoop);


            bool collide = CollideWithCandidateSet(VFGraph,Candidates[i],ChoosenPaths);
            if (collide)
            {
                To_Delete[i]=true;
                //std::cout<<"Collide Chosen"<<std::endl;
                continue;
            }

            if (UsePartitionNeeds)
            {
                //if (!SplitNonOKPartition(Candidates[i]))
                if (!SplitSomeNonOKPartition(VFGraph,Candidates[i],Partitions,FacePartitions,PartitionType))
                {
                    To_Delete[i]=true;
                    //std::cout<<"Not Needed by Partition "<<std::endl;
                    continue;
                }
            }

            //then the path has been chosen
            //            ChoosenPaths.push_back(Candidates[i]);
            //            assert(ChoosenPaths.back().PathNodes.size()>=2);
            //            //update vertices needs
            //            //if (UseVertNeeds)//CHECK!
            //            UpdateVertNeeds(CurrTrace);

            AddChoosen(Candidates[i]);

            //update distances
            UpdateDistancesWithLastChoosen();

            if (UsePartitionNeeds)
                UpdatePartitionsFromChoosen(true);

            //also update the
            CurrNodeDist=VFGraph.Distances();
            To_Delete[i]=true;
            DeleteCandidates(To_Delete);

            return true;
        }
        DeleteCandidates(To_Delete);
        return false;
    }

    void ChooseGreedyByDistance(bool UseVertNeeds,
                                bool UsePartitionNeeds)
    {

        UpdatePartitionsFromChoosen();
        CurrNodeDist.clear();
        InitNodeDistances();
        while (ChooseNextByDistance(UseVertNeeds,UsePartitionNeeds));
    }



    void GetFlatTangentNodes(std::vector<size_t> &TangentBorderNodes)
    {
        TangentBorderNodes.clear();

        for (size_t i=0;i<VertType.size();i++)
        {
            if (VertType[i]==TVFlat)
            {
                //get the nodes of a given vertex
                std::vector<size_t> NodesI;
                VertexFieldGraph<MeshType>::IndexNodes(i,NodesI);
                for (size_t j=0;j<NodesI.size();j++)
                {
                    if (NodeEmitterTypes[NodesI[j]]!=TVNone)continue;
                    if (NodeReceiverTypes[NodesI[j]]!=TVNone)continue;
                    TangentBorderNodes.push_back(NodesI[j]);
                }
            }
        }
    }

    void GetNodesType(TypeVert TVert,std::vector<size_t> &NodesSet)
    {
        NodesSet.clear();

        for (size_t i=0;i<VertType.size();i++)
        {
            if (VertType[i]==TVert)
            {
                //get the nodes of a given vertex
                std::vector<size_t> NodesI;
                VertexFieldGraph<MeshType>::IndexNodes(i,NodesI);
                NodesSet.insert(NodesSet.end(),NodesI.begin(),NodesI.end());
            }
        }
    }

    void GetEmitterType(const TypeVert EmitType,std::vector<size_t> &NodeEmitType)
    {
        NodeEmitType.clear();
        for (size_t i=0;i<NodeEmitterTypes.size();i++)
            if (NodeEmitterTypes[i]==EmitType)
                NodeEmitType.push_back(i);
    }

    void GetUnsatisfiedEmitterType(const TypeVert EmitType,std::vector<size_t> &NodeEmitType)
    {
        assert((EmitType==TVNarrow)||(EmitType==TVConcave));
        std::vector<size_t> TempEmit;
        GetEmitterType(EmitType,TempEmit);

        NodeEmitType.clear();
        for (size_t i=0;i<TempEmit.size();i++)
        {
            size_t IndexV=VertexFieldGraph<MeshType>::NodeVertI(TempEmit[i]);
            if (VerticesNeeds[IndexV]==0)continue;
            NodeEmitType.push_back(TempEmit[i]);
        }
    }


    void GetReceiverType(const TypeVert ReceiveType,
                         std::vector<size_t> &NodeReceiveType)
    {
        NodeReceiveType.clear();
        for (size_t i=0;i<NodeReceiverTypes.size();i++)
            if (NodeReceiverTypes[i]==ReceiveType)
                NodeReceiveType.push_back(i);
    }

    void GetUnsatisfiedReceiverType(const TypeVert ReceiveType,
                                    std::vector<size_t> &NodeReceiveType)
    {
        assert((ReceiveType==TVNarrow)||(ReceiveType==TVConcave));
        std::vector<size_t> TempReceive;
        GetReceiverType(ReceiveType,TempReceive);

        NodeReceiveType.clear();
        for (size_t i=0;i<TempReceive.size();i++)
        {
            size_t IndexV=VertexFieldGraph<MeshType>::NodeVertI(TempReceive[i]);
            if (VerticesNeeds[IndexV]==0)continue;
            NodeReceiveType.push_back(TempReceive[i]);
        }
    }

public:

    void GetActiveEmittersType(const TypeVert EmitType,std::vector<size_t> &NodeEmitType)
    {
        NodeEmitType.clear();
        //UpdateChosenReceivers();
        for (size_t i=0;i<NodeEmitterTypes.size();i++)
        {
            if (!VFGraph.IsActive(i))continue;
            if (NodeEmitterTypes[i]==EmitType)
                NodeEmitType.push_back(i);
        }
    }

    void GetActiveReceiversType(const TypeVert EmitType,std::vector<size_t> &NodeEmitType)
    {
        NodeEmitType.clear();
        //UpdateChosenReceivers();
        for (size_t i=0;i<NodeReceiverTypes.size();i++)
        {
            if (!VFGraph.IsActive(i))continue;
            if (NodeReceiverTypes[i]==EmitType)
                NodeEmitType.push_back(i);
        }
    }

private:

    void GetTracingConfiguration(const TypeVert FromType,
                                 const TypeVert ToType,
                                 const TraceType TracingType,
                                 std::vector<bool> &CanEmit,
                                 std::vector<bool> &CanReceive,
                                 std::vector<bool> &MustDisable)
    {
        CanEmit.clear();
        CanReceive.clear();
        MustDisable.clear();
        CanEmit.resize(VFGraph.NumNodes(),false);
        CanReceive.resize(VFGraph.NumNodes(),false);
        MustDisable.resize(VFGraph.NumNodes(),false);

        //all convex must be disabled always
        std::vector<size_t> ConvexNodes;
        GetNodesType(TVConvex,ConvexNodes);
        //GetConvexNodes(ConvexNodes);

        //leave it valid only when not all receivers and trace to flat
        if ((!AllReceivers)||(ToType!=TVFlat))
        {
            for (size_t i=0;i<ConvexNodes.size();i++)
                MustDisable[ConvexNodes[i]]=true;
        }

        //by default also disable the traced ones
        std::vector<size_t> ChoosenAndTangentNodes;
        //GetChoosenNodesAndTangent(ChoosenAndTangentNodes);
        GetCandidateNodesNodesAndTangent<MeshType>(ChoosenPaths,ChoosenAndTangentNodes);
        for (size_t i=0;i<ChoosenAndTangentNodes.size();i++)
            MustDisable[ChoosenAndTangentNodes[i]]=true;


        //also disable the border ones
        std::vector<size_t> TangentBorderNodes;
        GetFlatTangentNodes(TangentBorderNodes);
        if ((!AllReceivers)||(ToType!=TVFlat))
        {
            for (size_t i=0;i<TangentBorderNodes.size();i++)
                MustDisable[TangentBorderNodes[i]]=true;
        }

        if ((FromType==TVNarrow)&&(ToType==TVNarrow))
        {
            assert(TracingType!=TraceLoop);

            //get concave nodes
            std::vector<size_t> ConcaveNodes;
            //GetConcaveNodes(ConcaveNodes);
            GetNodesType(TVConcave,ConcaveNodes);

            //then disable them (cannot pass throught them)
            for (size_t i=0;i<ConcaveNodes.size();i++)
                MustDisable[ConcaveNodes[i]]=true;

            //set narrow emitter active and non satisfied
            std::vector<size_t> EmitterNodes;
            GetUnsatisfiedEmitterType(TVNarrow,EmitterNodes);
            for (size_t i=0;i<EmitterNodes.size();i++)
            {
                if (MustDisable[EmitterNodes[i]])continue;
                CanEmit[EmitterNodes[i]]=true;
            }

            //set narrow emitter receiver, active and non satisfied
            std::vector<size_t> ReceiverNodes;
            GetUnsatisfiedReceiverType(TVNarrow,ReceiverNodes);
            for (size_t i=0;i<ReceiverNodes.size();i++)
            {
                if (MustDisable[ReceiverNodes[i]])continue;
                CanReceive[ReceiverNodes[i]]=true;
            }

            return;
        }

        if ((FromType==TVNarrow)&&(ToType==TVConcave))
        {
            assert(TracingType!=TraceLoop);

            //narrow receivers
            std::vector<size_t> NarrowReceivers;
            GetReceiverType(TVNarrow,NarrowReceivers);
            for (size_t i=0;i<NarrowReceivers.size();i++)
                MustDisable[NarrowReceivers[i]]=true;

            //            //concave emitters
            //            std::vector<size_t> ConcaveEmitters;
            //            GetEmitterType(Concave,ConcaveEmitters);
            //            for (size_t i=0;i<ConcaveEmitters.size();i++)
            //                MustDisable[ConcaveEmitters[i]]=true;


            //set narrow emitter active and non satisfied
            std::vector<size_t> EmitterNodes;
            GetUnsatisfiedEmitterType(TVNarrow,EmitterNodes);
            for (size_t i=0;i<EmitterNodes.size();i++)
            {
                if (MustDisable[EmitterNodes[i]])continue;
                CanEmit[EmitterNodes[i]]=true;
            }

            //set narrow emitter receiver, active and non satisfied
            std::vector<size_t> ReceiverNodes;
            GetUnsatisfiedReceiverType(TVConcave,ReceiverNodes);
            for (size_t i=0;i<ReceiverNodes.size();i++)
            {
                if (MustDisable[ReceiverNodes[i]])continue;
                CanReceive[ReceiverNodes[i]]=true;
            }

            //concave emitters disable if not receive
            std::vector<size_t> ConcaveEmitters;
            GetEmitterType(TVConcave,ConcaveEmitters);
            for (size_t i=0;i<ConcaveEmitters.size();i++)
            {
                if (CanReceive[ConcaveEmitters[i]])continue;
                MustDisable[ConcaveEmitters[i]]=true;
            }

            return;
        }

        if ((FromType==TVNarrow)&&(ToType==TVFlat))
        {
            assert(TracingType!=TraceLoop);

            //narrow receivers
            std::vector<size_t> NarrowReceivers;
            GetReceiverType(TVNarrow,NarrowReceivers);

            //concave emitters
            std::vector<size_t> ConcaveNodes;
            //GetConcaveNodes(ConcaveNodes);
            GetNodesType(TVConcave,ConcaveNodes);

            //            if (!AllReceivers)
            //            {
            for (size_t i=0;i<NarrowReceivers.size();i++)
                MustDisable[NarrowReceivers[i]]=true;

            for (size_t i=0;i<ConcaveNodes.size();i++)
                MustDisable[ConcaveNodes[i]]=true;
            //           }

            //set narrow emitter active and non satisfied
            std::vector<size_t> EmitterNodes;
            GetUnsatisfiedEmitterType(TVNarrow,EmitterNodes);
            for (size_t i=0;i<EmitterNodes.size();i++)
            {
                if (MustDisable[EmitterNodes[i]])continue;
                CanEmit[EmitterNodes[i]]=true;
            }

            //set flat emitter receiver, active and non satisfied
            std::vector<size_t> ReceiverNodes;
            GetReceiverType(TVFlat,ReceiverNodes);
            for (size_t i=0;i<ReceiverNodes.size();i++)
            {
                if (MustDisable[ReceiverNodes[i]])continue;
                CanReceive[ReceiverNodes[i]]=true;
            }

            if (AllReceivers)
            {
                for (size_t i=0;i<ConvexNodes.size();i++)
                    CanReceive[ConvexNodes[i]]=true;
                for (size_t i=0;i<TangentBorderNodes.size();i++)
                    CanReceive[TangentBorderNodes[i]]=true;
            }
            return;
        }

        if ((FromType==TVConcave)&&(ToType==TVConcave))
        {
            assert(TracingType!=TraceLoop);

            //narrow receivers
            std::vector<size_t> NarrowReceivers;
            GetReceiverType(TVNarrow,NarrowReceivers);
            for (size_t i=0;i<NarrowReceivers.size();i++)
                MustDisable[NarrowReceivers[i]]=true;

            //narrow emitters
            std::vector<size_t> NarrowEmitters;
            GetEmitterType(TVNarrow,NarrowEmitters);
            for (size_t i=0;i<NarrowEmitters.size();i++)
                MustDisable[NarrowEmitters[i]]=true;


            //set concave emitter active and non satisfied
            std::vector<size_t> EmitterNodes;
            GetUnsatisfiedEmitterType(TVConcave,EmitterNodes);
            for (size_t i=0;i<EmitterNodes.size();i++)
            {
                if (MustDisable[EmitterNodes[i]])continue;
                CanEmit[EmitterNodes[i]]=true;
            }

            //set narrow emitter receiver, active and non satisfied
            std::vector<size_t> ReceiverNodes;
            GetUnsatisfiedReceiverType(TVConcave,ReceiverNodes);
            for (size_t i=0;i<ReceiverNodes.size();i++)
            {
                if (MustDisable[ReceiverNodes[i]])continue;
                CanReceive[ReceiverNodes[i]]=true;
            }
            return;
        }

        if ((FromType==TVConcave)&&(ToType==TVFlat))
        {
            assert(TracingType!=TraceLoop);

            //narrow receivers
            std::vector<size_t> NarrowReceivers;
            GetReceiverType(TVNarrow,NarrowReceivers);
            for (size_t i=0;i<NarrowReceivers.size();i++)
                MustDisable[NarrowReceivers[i]]=true;

            //narrow emitters
            std::vector<size_t> NarrowEmitters;
            GetEmitterType(TVNarrow,NarrowEmitters);
            for (size_t i=0;i<NarrowEmitters.size();i++)
                MustDisable[NarrowEmitters[i]]=true;

            //            //concave receivers
            //            std::vector<size_t> ConcaveReceivers;
            //            GetReceiverType(Concave,ConcaveReceivers);
            //            for (size_t i=0;i<ConcaveReceivers.size();i++)
            //                MustDisable[ConcaveReceivers[i]]=true;


            //set concave emitter active and non satisfied
            std::vector<size_t> EmitterNodes;
            GetUnsatisfiedEmitterType(TVConcave,EmitterNodes);
            for (size_t i=0;i<EmitterNodes.size();i++)
            {
                if (MustDisable[EmitterNodes[i]])continue;
                CanEmit[EmitterNodes[i]]=true;
            }

            //set flat receiver, active and non satisfied
            std::vector<size_t> ReceiverNodes;
            GetReceiverType(TVFlat,ReceiverNodes);
            for (size_t i=0;i<ReceiverNodes.size();i++)
            {
                if (MustDisable[ReceiverNodes[i]])continue;
                CanReceive[ReceiverNodes[i]]=true;
            }

            //concave receivers
            std::vector<size_t> ConcaveReceivers;
            GetReceiverType(TVConcave,ConcaveReceivers);
            for (size_t i=0;i<ConcaveReceivers.size();i++)
            {
                if (CanEmit[ConcaveReceivers[i]])continue;
                MustDisable[ConcaveReceivers[i]]=true;
            }

            if (AllReceivers)
            {
                for (size_t i=0;i<ConvexNodes.size();i++)
                    CanReceive[ConvexNodes[i]]=true;

                for (size_t i=0;i<TangentBorderNodes.size();i++)
                    CanReceive[TangentBorderNodes[i]]=true;
            }
            return;
        }


        if ((FromType==TVFlat)&&(ToType==TVFlat))
        {
            //narrow receivers
            std::vector<size_t> NarrowReceivers;
            GetReceiverType(TVNarrow,NarrowReceivers);
            for (size_t i=0;i<NarrowReceivers.size();i++)
                MustDisable[NarrowReceivers[i]]=true;

            //narrow emitters
            std::vector<size_t> NarrowEmitters;
            GetEmitterType(TVNarrow,NarrowEmitters);
            for (size_t i=0;i<NarrowEmitters.size();i++)
                MustDisable[NarrowEmitters[i]]=true;

            //narrow emitters
            std::vector<size_t> ConcaveNodes;
            //GetConcaveNodes(ConcaveNodes);
            GetNodesType(TVConcave,ConcaveNodes);

            for (size_t i=0;i<ConcaveNodes.size();i++)
                MustDisable[ConcaveNodes[i]]=true;

            //set concave emitter active and non satisfied
            std::vector<size_t> EmitterNodes;
            GetEmitterType(TVFlat,EmitterNodes);
            for (size_t i=0;i<EmitterNodes.size();i++)
            {
                if (MustDisable[EmitterNodes[i]])continue;
                CanEmit[EmitterNodes[i]]=true;
            }

            std::vector<size_t> ReceiverNodes;
            GetReceiverType(TVFlat,ReceiverNodes);
            for (size_t i=0;i<ReceiverNodes.size();i++)
            {
                if (MustDisable[ReceiverNodes[i]])continue;
                CanReceive[ReceiverNodes[i]]=true;
            }
            return;
        }

        if (TracingType==TraceLoop)
        {
            assert(FromType==TVInternal);
            assert(ToType==TVInternal);

            //get concave nodes
            std::vector<size_t> ConcaveNodes;
            //GetConcaveNodes(ConcaveNodes);
            GetNodesType(TVConcave,ConcaveNodes);

            //narrow emitters
            std::vector<size_t> NarrowEmittersNodes;
            GetEmitterType(TVNarrow,NarrowEmittersNodes);

            //narrow receivers
            std::vector<size_t> NarrowReceiversNodes;
            //GetNarrowReceivers(NarrowReceiversNodes);
            GetReceiverType(TVNarrow,NarrowReceiversNodes);

            for (size_t i=0;i<ConcaveNodes.size();i++)
                MustDisable[ConcaveNodes[i]]=true;

            for (size_t i=0;i<NarrowEmittersNodes.size();i++)
                MustDisable[NarrowEmittersNodes[i]]=true;

            for (size_t i=0;i<NarrowReceiversNodes.size();i++)
                MustDisable[NarrowReceiversNodes[i]]=true;

            //set internal emitter
            std::vector<size_t> EmitterNodes;
            GetEmitterType(TVInternal,EmitterNodes);
            for (size_t i=0;i<EmitterNodes.size();i++)
            {
                if (MustDisable[EmitterNodes[i]]==false)
                    CanEmit[EmitterNodes[i]]=true;
            }

            //get also flat ones
            //            if (TraceLoopsBorders)
            //            {
            GetEmitterType(TVFlat,EmitterNodes);
            for (size_t i=0;i<EmitterNodes.size();i++)
            {
                if (MustDisable[EmitterNodes[i]]==false)
                    CanEmit[EmitterNodes[i]]=true;
            }
            //}
            return;
        }
        if (DebugMsg)
        {
            std::cout<<"****************************************"<<std::endl;
            PrintConfiguration(FromType,ToType,TracingType);
            std::cout<<"****************************************"<<std::endl;
        }
        assert(0);
    }

    void InitNodeDistances()
    {
        std::vector<size_t> Sources;

        for (size_t i=0;i<VFGraph.NumNodes();i++)
        {
            if (VFGraph.IsActive(i))continue;
            Sources.push_back(i);
        }
        for (size_t i=0;i<ChoosenPaths.size();i++)
        {
            std::vector<size_t> CurrPath=ChoosenPaths[i].PathNodes;
            Sources.insert(Sources.end(),CurrPath.begin(),CurrPath.end());
            VertexFieldGraph<MeshType>::TangentNodes(CurrPath);
            Sources.insert(Sources.end(),CurrPath.begin(),CurrPath.end());
        }

        //then add the singularities
        if (away_from_singular)
            Sources.insert(Sources.end(),VFGraph.SingNodes.begin(),VFGraph.SingNodes.end());

        //        std::vector<size_t> FlatEmit;
        //        GetFlatEmitters(FlatEmit);
        //        Sources.insert(Sources.end(),FlatEmit.begin(),FlatEmit.end());

        std::sort(Sources.begin(),Sources.end());
        std::vector<size_t>::iterator it;
        it = std::unique (Sources.begin(),Sources.end());
        Sources.resize( std::distance(Sources.begin(),it) );

        if (Sources.size()==0){
            int MinIdx,MaxIdx;

            if (DebugMsg)
                std::cout<<"NO INITIAL FEATURE"<<std::endl;

            GetFarthestCoupleVertx(Mesh(),MinIdx,MaxIdx);

            if (DebugMsg)
                std::cout<<"Initialize with Points "<< MinIdx <<" and "<< MaxIdx <<std::endl;

            Sources.push_back(MinIdx);
            Sources.push_back(MaxIdx);
        }
        std::vector<bool> IsActive;
        VFGraph.IsActiveNodes(IsActive);
        VFGraph.SetAllActive();
        VertexFieldQuery<MeshType>::UpdateDistancesFrom(VFGraph,Sources,Drift);
        VFGraph.SetActiveNodes(IsActive);
        CurrNodeDist=VFGraph.Distances();
    }

    void RemoveNonActive(std::vector<size_t> &Nodes)
    {
        std::vector<size_t> SwapNodes;
        for (size_t i=0;i<Nodes.size();i++)
            if (VFGraph.IsActive(Nodes[i]))
                SwapNodes.push_back(Nodes[i]);
        Nodes=SwapNodes;
    }

    void UpdateDistancesWithLastChoosen()
    {
        std::vector<size_t> Sources;
        std::vector<size_t> CurrPath=ChoosenPaths.back().PathNodes;
        Sources=CurrPath;
        VertexFieldGraph<MeshType>::TangentNodes(CurrPath);
        Sources.insert(Sources.end(),CurrPath.begin(),CurrPath.end());
        RemoveNonActive(Sources);
        VertexFieldQuery<MeshType>::UpdateDistancesFrom(VFGraph,Sources,Drift,&CurrNodeDist);

        //InitDistances();
    }

    void PrintConfiguration(TypeVert FromType,
                            TypeVert ToType,
                            TraceType TrType)
    {
        std::cout<<std::endl<<std::endl<<"** INFO **"<<std::endl;
        if (FromType==TVNarrow)
            std::cout<<"Tracing From Narrow ";
        if (FromType==TVConcave)
            std::cout<<"Tracing From Concave ";
        if (FromType==TVFlat)
            std::cout<<"Tracing From Flat ";
        if (FromType==TVInternal)
            std::cout<<"Tracing From Internal ";
        if (FromType==TVChoosen)
            std::cout<<"Tracing From Choosen ";

        if (ToType==TVNarrow)
            std::cout<<"To Narrow ";
        if (ToType==TVConcave)
            std::cout<<"To Concave ";
        if (ToType==TVFlat)
            std::cout<<"To Flat ";
        if (ToType==TVInternal)
            std::cout<<"To Internal ";
        if (ToType==TVChoosen)
            std::cout<<"To Choosen ";


        if (TrType==TraceDirect)
            std::cout<<"Method Direct "<<std::endl;
        if (TrType==DijkstraReceivers)
            std::cout<<"Method Dijkstra "<<std::endl;
        if (TrType==TraceLoop)
            std::cout<<"Method Loop "<<std::endl;
        std::cout<<"** END **"<<std::endl<<std::endl<<std::endl;
    }

    bool JoinConnection(TypeVert FromType,
                        TypeVert ToType,
                        TraceType TrType,
                        bool UsePartitionNeeds=false)
    {
        if (DebugMsg)
            PrintConfiguration(FromType,ToType,TrType);

        //Candidates.clear();
        assert((FromType!=TVChoosen)&&(ToType!=TVChoosen));
        //        if ((FromType==Choosen)||(ToType==Choosen))
        //            UpdateChosenReceivers();

        std::vector<bool> CanEmit,CanReceive,MustDisable;
        GetTracingConfiguration(FromType,ToType,TrType,CanEmit,CanReceive,MustDisable);

        VFGraph.SetAllActive();
        VFGraph.SetDisabledNodes(MustDisable);
        if ((FromType==TVNarrow)||(ToType==TVNarrow)
                ||(FromType==TVConcave)||(ToType==TVConcave))
            return(TraceFrom(FromType,ToType,TrType,CanEmit,CanReceive,Shortest));
        else
            return(TraceFrom(FromType,ToType,TrType,CanEmit,CanReceive,Fartest));
    }


private:

    void UpdatePartitionType(size_t Index)
    {
        assert(Index<PartitionCorners.size());
        assert(Index<Partitions.size());

        //REMOVE THE REST
        if (PatchInfos[Index].NumEmitters>0)
        {
            PartitionType[Index]=HasEmitter;
            return;
        }

        if (PatchInfos[Index].NumCorners<(int)MinVal)
        {
            PartitionType[Index]=LowCorners;
            return;
        }

        if (PatchInfos[Index].NumCorners>(int)MaxVal)
        {
            PartitionType[Index]=HighCorners;
            return;
        }

        if (PatchInfos[Index].Genus!=1)
        {
            PartitionType[Index]=NonDisk;
            return;
        }

        if ((CClarkability>0)&&(!PatchInfos[Index].CClarkability))
        {
            PartitionType[Index]=MaxCClarkability;
            return;
        }

        if ((match_valence)&&(PatchInfos[Index].ExpectedValence!=PatchInfos[Index].NumCorners))
        {
            PartitionType[Index]=NonMatchValence;
            return;
        }

        PartitionType[Index]=IsOK;
    }

    //TO BE DELETED
    void InitPartitionsType()
    {
        PartitionType.clear();
        PartitionType.resize(Partitions.size(),IsOK);

        InitEdgeL();

        GetPatchInfo(Mesh(),Partitions,PartitionCorners,VerticesNeeds,EdgeL,
                     PatchInfos,avgEdge*CClarkability);

        for (size_t i=0;i<Partitions.size();i++)
            UpdatePartitionType(i);


    }



public:

    void UpdatePartitionsFromChoosen(bool UpdateType=true)
    {
        bool IsOk=SelectMeshPatchBorders(Mesh(),ChoosenPaths);//SelectBorders();

        if (!IsOk)
        {
            MeshType traceMesh;
            std::vector<bool> Selected(ChoosenPaths.size(),false);
            //Selected[i]=true;
            MeshTraces(VFGraph,ChoosenPaths,Selected,traceMesh);
            //vcg::tri::io::ExporterPLY<MeshType>::Save(VFGraph.Mesh(),"double_direction_domain.ply");
            vcg::tri::io::ExporterPLY<MeshType>::Save(traceMesh,"double_direction_error.ply",vcg::tri::io::Mask::IOM_FACECOLOR);
        }

        assert(IsOk);
        std::vector<size_t> StartF;
        for (size_t i=0;i<Mesh().face.size();i++)
            StartF.push_back(i);

        RetrievePatchesFromSelEdges(Mesh(),StartF,Partitions);


        DerivePerFacePartition(Mesh(),Partitions,FacePartitions);

        FindPartitionsCorners(VFGraph,VertType,ChoosenPaths,Partitions,PartitionCorners);

        if (UpdateType)
            InitPartitionsType();
    }

public:


    void CopyParametersFrom(const PatchTracer &PTracer)
    {
        split_on_removal=PTracer.split_on_removal;
        //avoid_increase_valence=PTracer.avoid_increase_valence;
        //avoid_collapse_irregular=PTracer.avoid_collapse_irregular;
        //max_lenght_distortion=PTracer.max_lenght_distortion;
        //max_lenght_variance=PTracer.max_lenght_variance;
        sample_ratio=PTracer.sample_ratio;
        //max_patch_area=PTracer.max_patch_area;
        MinVal=PTracer.MinVal;
        MaxVal=PTracer.MaxVal;
        DebugMsg=PTracer.DebugMsg;
        away_from_singular=PTracer.away_from_singular;
        CClarkability=PTracer.CClarkability;
        match_valence=PTracer.match_valence;
    }


private:

    bool HasIncompleteEmitter()
    {
        for (size_t i=0;i<PartitionType.size();i++)
            if (PartitionType[i]==HasEmitter)return true;
        return false;
    }

    bool HasTerminated()
    {
        for (size_t i=0;i<PartitionType.size();i++)
            if (PartitionType[i]!=IsOK)return false;
        return true;
    }

    bool PathHasConcaveNarrowVert(size_t IndexPath)
    {
        for (size_t j=0;j<ChoosenPaths[IndexPath].PathNodes.size();j++)
        {
            size_t IndexN=ChoosenPaths[IndexPath].PathNodes[j];
            size_t IndexV=VertexFieldGraph<MeshType>::NodeVertI(IndexN);
            if ((VertType[IndexV]==TVNarrow)||(VertType[IndexV]==TVConcave))
                return true;
        }
        return false;
    }

    bool PathHasTJunction(size_t IndexPath)
    {
        std::vector<size_t> DirV[2];
        DirV[0].resize(Mesh().vert.size(),0);
        DirV[1].resize(Mesh().vert.size(),0);
        for (size_t i=0;i<ChoosenPaths.size();i++)
        {
            if (ChoosenPaths[i].PathNodes.size()==0)continue;
            if (ChoosenPaths[i].IsLoop)continue;
            size_t NodeI0=ChoosenPaths[i].PathNodes[0];
            size_t NodeI1=ChoosenPaths[i].PathNodes.back();
            size_t Dir0=VertexFieldGraph<MeshType>::NodeDirI(NodeI0)%2;
            size_t Dir1=VertexFieldGraph<MeshType>::NodeDirI(NodeI1)%2;
            size_t VertI0=VertexFieldGraph<MeshType>::NodeVertI(NodeI0);
            size_t VertI1=VertexFieldGraph<MeshType>::NodeVertI(NodeI1);
            DirV[Dir0][VertI0]+=1;
            DirV[Dir1][VertI1]+=1;
        }

        for (size_t i=0;i<ChoosenPaths[IndexPath].PathNodes.size();i++)
        {
            size_t NodeI=ChoosenPaths[IndexPath].PathNodes[i];
            size_t VertI=VertexFieldGraph<MeshType>::NodeVertI(NodeI);
            size_t Dir=VertexFieldGraph<MeshType>::NodeDirI(NodeI)%2;
            size_t CrossDir=(Dir+1)%2;
            //if cross with a sharp feature then not check
            if (Mesh().vert[VertI].IsB())continue;
            if (DirV[CrossDir][VertI]==1)return true;
            if ((DirV[Dir][VertI]==2)&&(DirV[CrossDir][VertI]==0))return true;
        }
        return false;
    }

    bool HasPathDeadEnd()
    {
        std::vector<size_t> PathV;
        PathV.resize(Mesh().vert.size(),0);
        for (size_t i=0;i<ChoosenPaths.size();i++)
        {
            if (ChoosenPaths[i].PathNodes.size()==0)continue;
            size_t limit=ChoosenPaths[i].PathNodes.size()-1;
            if (ChoosenPaths[i].IsLoop)limit++;
            size_t node_num=ChoosenPaths[i].PathNodes.size();
            for (size_t j=0;j<limit;j++)
            {
                size_t NodeI0=ChoosenPaths[i].PathNodes[j];
                size_t NodeI1=ChoosenPaths[i].PathNodes[(j+1)%node_num];
                size_t VertI0=VertexFieldGraph<MeshType>::NodeVertI(NodeI0);
                size_t VertI1=VertexFieldGraph<MeshType>::NodeVertI(NodeI1);
                PathV[VertI0]++;
                PathV[VertI1]++;
            }
        }
        for (size_t i=0;i<PathV.size();i++)
        {
            if (Mesh().vert[i].IsB())continue;
            if (PathV[i]==1)return true;
        }
        return false;
    }

    //TO BE DELETED
    void GetCurrentConfigurationAround(const std::vector<vcg::face::Pos<FaceType> > &FacesPath,
                                       bool twoSides,
                                       size_t &NonProper,
                                       size_t &Valence3,
                                       size_t &Valence5,
                                       size_t &Valence6)
    {
        NonProper=0;
        Valence3=0;
        Valence5=0;
        Valence6=0;
        //get the indexes of faces
        std::vector<size_t> IdxFaces;
        for (size_t i=0;i<FacesPath.size();i++)
        {
            IdxFaces.push_back(vcg::tri::Index(Mesh(),FacesPath[i].F()));
            if (twoSides)
                IdxFaces.push_back(vcg::tri::Index(Mesh(),FacesPath[i].FFlip()));
        }
        //std::cout<<"2-Selected Path Pos "<<IdxFaces.size()<<std::endl;
        //then retrieve partitions
        //RetrievePartitionsFaces(IdxFaces);
        RetrievePatchesFromSelEdges(Mesh(),IdxFaces,Partitions);
        //std::cout<<"There are "<<PartitionType.size()<<" partitions"<<std::endl;
        //find corners
        //FindPartitionsCorners();
        FindPartitionsCorners<MeshType>(VFGraph,VertType,ChoosenPaths,Partitions,PartitionCorners);
        //find type
        InitPartitionsType();
        //std::cout<<"There are "<<PartitionType.size()<<" partitions"<<std::endl;
        for (size_t i=0;i<PartitionType.size();i++)
        {
            switch(PartitionType[i]) {
            case LowCorners:
                NonProper++;
                break;
            case HighCorners:
                NonProper++;
                break;
            case NonDisk:
                NonProper++;
                break;
            case HasEmitter:
                NonProper++;
                break;
            default : break;
            }
        }
        for (size_t i=0;i<PartitionCorners.size();i++)
        {
            if (PartitionCorners[i].size()==3)Valence3++;
            if (PartitionCorners[i].size()==5)Valence5++;
            if (PartitionCorners[i].size()==6)Valence6++;
        }
    }

    void UpdatePatchAround(const std::vector<vcg::face::Pos<FaceType> > &FacesPath)
    {
        //get the indexes of faces
        std::vector<size_t> IdxFaces;
        for (size_t i=0;i<FacesPath.size();i++)
        {
            IdxFaces.push_back(vcg::tri::Index(Mesh(),FacesPath[i].F()));
            IdxFaces.push_back(vcg::tri::Index(Mesh(),FacesPath[i].FFlip()));
        }
        //then retrieve partitions
        RetrievePatchesFromSelEdges(Mesh(),IdxFaces,Partitions);
        //find corners
        FindPartitionsCorners<MeshType>(VFGraph,VertType,ChoosenPaths,Partitions,PartitionCorners);
        //find type
        InitPartitionsType();
    }

    bool RemoveIfPossible(size_t IndexPath)
    {
        if (ChoosenPaths[IndexPath].PathNodes.size()==0)return false;
        if (ChoosenPaths[IndexPath].Unremovable){std::cout<<"Unremoovable"<<std::endl;return false;}
        //if it includes a concave or narrow then cannot remove
        if (PathHasConcaveNarrowVert(IndexPath))return false;
        //check if have t junction in the middle
        if (PathHasTJunction(IndexPath))return false;
        //CHECK ENDPOINTS!
        assert(IndexPath<ChoosenPaths.size());

        //get the old configuration
        std::vector<vcg::face::Pos<FaceType> > FacesPath;
        VFGraph.GetNodesPos(ChoosenPaths[IndexPath].PathNodes,ChoosenPaths[IndexPath].IsLoop,FacesPath);

        UpdatePatchAround(FacesPath);
        std::vector<PatchInfo<ScalarType> > PatchInfos0=PatchInfos;
        std::vector<std::vector<size_t> > FacePatches0=Partitions;

        //test removal
        CandidateTrace OldTr=ChoosenPaths[IndexPath];
        ChoosenPaths[IndexPath].PathNodes.clear();

        //check Tjunctions
        if (HasPathDeadEnd())
        {
            //restore
            ChoosenPaths[IndexPath]=OldTr;
            Mesh().SelectPos(FacesPath,true);
            return false;
        }
        //deselect
        Mesh().SelectPos(FacesPath,false);
        UpdatePatchAround(FacesPath);
        std::vector<PatchInfo<ScalarType> > PatchInfos1=PatchInfos;
        std::vector<std::vector<size_t> > FacePatches1=Partitions;

        for (size_t i=0;i<PatchInfos1.size();i++)
        {
            if ((PatchInfos1[i].NumCorners<MIN_ADMITTIBLE)
                    ||(PatchInfos1[i].NumCorners>MAX_ADMITTIBLE))
            {
                //restore
                ChoosenPaths[IndexPath]=OldTr;
                Mesh().SelectPos(FacesPath,true);
                return false;
            }
        }

        bool CanRemove=true;
        //std::vector<typename MeshType::ScalarType> QThresold;
        CanRemove=BetterConfiguaration(Mesh(),FacePatches0,FacePatches1,
                                       PatchInfos0,PatchInfos1,MinVal,
                                       MaxVal,CClarkability,avgEdge,
                                       match_valence);

        if (!CanRemove)
        {
            //restore
            ChoosenPaths[IndexPath]=OldTr;
            Mesh().SelectPos(FacesPath,true);
            return false;
        }
        return true;
    }

    void RemoveEmptyPaths()
    {
        std::vector<CandidateTrace> ChoosenPathsSwap;
        //remove the ones with no nodes
        for (size_t i=0;i<ChoosenPaths.size();i++)
        {
            if (ChoosenPaths[i].PathNodes.size()==0)continue;
            ChoosenPathsSwap.push_back(ChoosenPaths[i]);
        }
        ChoosenPaths=ChoosenPathsSwap;
    }

    bool RemoveIteration()
    {
        bool HasRemoved=false;
        if (DebugMsg)
            std::cout<<"REMOVING ITERATION"<<std::endl;

        for (int i=ChoosenPaths.size()-1;i>=0;i--)
        {
            //std::cout<<"Removing "<<i<<" if Possible"<<std::endl;
            HasRemoved|=RemoveIfPossible(i);
        }

        MergeContiguousPaths(ChoosenPaths);

        if (DebugMsg)
            std::cout<<"DONE!"<<std::endl;

        //        if (HasRemoved)
        //        {
        //        MergeContiguousPaths();
        RemoveEmptyPaths();
        return HasRemoved;
    }

    void SplitIntoSubPathsBySel(const CandidateTrace &ToSplit,
                                std::vector<CandidateTrace> &Portions)
    {
        Portions.clear();
        int StartI=0;
        if (ToSplit.IsLoop)
        {
            StartI=-1;
            for (size_t j=0;j<ToSplit.PathNodes.size();j++)
            {
                size_t IndexV=VertexFieldGraph<MeshType>::NodeVertI(ToSplit.PathNodes[j]);
                if (Mesh().vert[IndexV].IsS())
                {
                    StartI=j;
                    break;
                }
            }

            //no crossing at all
            if (StartI==-1)
            {
                Portions.push_back(ToSplit);
                return;
            }
        }
        assert(StartI>=0);
        assert(StartI<ToSplit.PathNodes.size());
        size_t numNodes=ToSplit.PathNodes.size();
        std::vector<std::vector<size_t> > SubPaths;
        SubPaths.resize(1);
        int CurrI=StartI;
        do
        {
            size_t CurrN=ToSplit.PathNodes[CurrI];
            size_t IndexV=VertexFieldGraph<MeshType>::NodeVertI(CurrN);
            size_t NextI=(CurrI+1)%numNodes;

            //then add a new subpath
            if ((Mesh().vert[IndexV].IsS())&&(CurrI!=StartI))
            {
                SubPaths.back().push_back(CurrN);
                SubPaths.resize(SubPaths.size()+1);
            }
            SubPaths.back().push_back(CurrN);
            CurrI=NextI;
        }while (CurrI!=StartI);

        //add the first one in case is loop
        if (ToSplit.IsLoop)
            SubPaths.back().push_back(ToSplit.PathNodes[StartI]);

        //in case not loop and last one was selected then remove the last
        if ((!ToSplit.IsLoop)&&(SubPaths.back().size()==1))
            SubPaths.pop_back();

        //        std::cout<<"Max Size "<<SubPaths.size()<<std::endl;
        //        if (ToSplit.IsLoop)
        //            std::cout<<"Loop "<<std::endl;
        for (size_t i=0;i<SubPaths.size();i++)
        {
            //            std::cout<<"Indeex "<<i<<std::endl;
            //            std::cout<<"Size "<<SubPaths[i].size()<<std::endl;
            assert(SubPaths[i].size()>=2);
            size_t Node0=SubPaths[i][0];
            size_t Node1=SubPaths[i].back();
            size_t IndexV0=VertexFieldGraph<MeshType>::NodeVertI(Node0);
            size_t IndexV1=VertexFieldGraph<MeshType>::NodeVertI(Node1);
            TypeVert T0=VertType[IndexV0];
            TypeVert T1=VertType[IndexV1];
            Portions.push_back(CandidateTrace(T0,T1,ToSplit.TracingMethod,Node0));
            Portions.back().PathNodes=SubPaths[i];
            Portions.back().IsLoop=false;
            Portions.back().Updated=true;
        }
        //case of a loop with a single intersection, then keep the original
        if ((Portions.size()==1)&&(ToSplit.IsLoop))
        {
            Portions.clear();
            Portions.push_back(ToSplit);
        }
        if (ToSplit.Unremovable)
            for (size_t i=0;i<Portions.size();i++)
                Portions[i].Unremovable=true;
    }

    void SelectCrossIntersections()
    {
        vcg::tri::UpdateSelection<MeshType>::VertexClear(Mesh());
        vcg::tri::UpdateQuality<MeshType>::VertexConstant(Mesh(),0);
        //count occourence
        for (size_t i=0;i<ChoosenPaths.size();i++)
            for (size_t j=0;j<ChoosenPaths[i].PathNodes.size();j++)
            {
                size_t IndexV=VertexFieldGraph<MeshType>::NodeVertI(ChoosenPaths[i].PathNodes[j]);
                Mesh().vert[IndexV].Q()+=1;
            }
        for (size_t i=0;i<Mesh().vert.size();i++)
        {
            if (Mesh().vert[i].Q()>=2)Mesh().vert[i].SetS();
            if (Mesh().vert[i].IsB() && (Mesh().vert[i].Q()>0))Mesh().vert[i].SetS();
        }
    }


public:

    void SplitIntoSubPaths()
    {
        SelectCrossIntersections();
        std::vector<CandidateTrace> NewTraces;
        for (size_t i=0;i<ChoosenPaths.size();i++)
        {
            std::vector<CandidateTrace> Portions;
            SplitIntoSubPathsBySel(ChoosenPaths[i],Portions);
            NewTraces.insert(NewTraces.end(),Portions.begin(),Portions.end());
        }
        ChoosenPaths.clear();
        ChoosenPaths=NewTraces;

        for (size_t i=0;i<ChoosenPaths.size();i++)
        {
            //std::cout<<ChoosenPaths[i].PathNodes.size()<<std::endl;
            assert(ChoosenPaths[i].PathNodes.size()>=2);
        }

        UpdatePartitionsFromChoosen();
        vcg::tri::UpdateSelection<MeshType>::Clear(Mesh());
        ColorByPartitions();
    }

public:

    MeshType &Mesh()
    {
        return VFGraph.Mesh();
    }

    //TO BE DELETED
    struct PatchInfoType
    {
        size_t LowC;
        size_t HighC;
        size_t NonDiskLike;
        size_t HasEmit;
        size_t NumPatchs;
        size_t SizePatches[8];
    };

    //TO BE DELETED
    void GetInfo(PatchInfoType &PInfo)
    {
        PInfo.LowC=0;
        PInfo.HighC=0;
        PInfo.NonDiskLike=0;
        PInfo.HasEmit=0;
        PInfo.NumPatchs=PartitionType.size();
        PInfo.SizePatches[0]=0;
        PInfo.SizePatches[1]=0;
        PInfo.SizePatches[2]=0;
        PInfo.SizePatches[3]=0;
        PInfo.SizePatches[4]=0;
        PInfo.SizePatches[5]=0;
        PInfo.SizePatches[6]=0;
        PInfo.SizePatches[7]=0;

        for (size_t i=0;i<PartitionType.size();i++)
        {
            switch(PartitionType[i]) {
            case LowCorners  :
                PInfo.LowC++;
                break;
            case HighCorners  :
                PInfo.HighC++;
                break;
            case NonDisk  :
                PInfo.NonDiskLike++;
                break;
            case HasEmitter  :
                PInfo.HasEmit++;
                break;
            default:
                break;
                //            case MoreSing  :
                //                HasMoreSing++;
                //                break;
                //            case IsOK  :
                //                break;
            }
        }
        for (size_t i=0;i<PartitionCorners.size();i++)
        {
            size_t numC=PartitionCorners[i].size();
            numC=std::min(numC,(size_t)7);
            PInfo.SizePatches[numC]++;
        }
    }

    void WriteInfo()
    {

        PatchInfoType PInfo;
        GetInfo(PInfo);
        std::cout<<"***FINAL STATS***"<<std::endl;
        std::cout<<"* Num Patches "<<PInfo.NumPatchs<<std::endl;
        std::cout<<"* Low Valence Patches "<<PInfo.LowC<<std::endl;
        std::cout<<"* High Valence Patches "<<PInfo.HighC<<std::endl;
        std::cout<<"* Non Disk Like Patches "<<PInfo.NonDiskLike<<std::endl;
        std::cout<<"* With Emitters Patches "<<PInfo.HasEmit<<std::endl;
        for (size_t i=0;i<8;i++)
            std::cout<<"* Patch with  "<<i<<" corners are: "<<PInfo.SizePatches[i]<<std::endl;
        //std::cout<<"* With More Singularities "<<HasMoreSing<<std::endl;
    }

    //    ScalarType PatchDistortion(size_t IndexP)
    //    {

    //    }
    void SmoothPatches(size_t Steps=3,typename MeshType::ScalarType Damp=0.5)
    {
        SmoothMeshPatches(Mesh(),FacePartitions,Steps,Damp);
    }

    void SetAllRemovable()
    {
        for (size_t i=0;i<ChoosenPaths.size();i++)
            ChoosenPaths[i].Unremovable=false;
    }

    void RemovePaths()//bool DoSmooth=true)
    {
        //max_patch_area=MeshArea(Mesh())*0.02;
        std::vector<std::vector<vcg::face::Pos<FaceType> > > PathPos;

        //select pos
        vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(Mesh());
        GetPathPos(VFGraph,ChoosenPaths,PathPos);
        Mesh().SelectPos(PathPos,true);

//        if (DoSmooth)
//            SmoothPatches(10);

        if (DebugMsg)
            std::cout<<"Removing..."<<std::endl;
        while (RemoveIteration()){}

        UpdatePartitionsFromChoosen(true);
        ColorByPartitions();

        if (DebugMsg)
            WriteInfo();

        //MergeContiguousPaths();
        //ColorByPatchQuality();
    }

    void GetPatchMesh(const size_t &IndexPatch,
                      MeshType &PatchMesh)
    {
        GetMeshFromPatch(Mesh(),IndexPatch,Partitions,PatchMesh);
    }

    void GetVisualCornersPos(std::vector<CoordType> &PatchCornerPos)
    {
        PatchCornerPos.clear();
        //        for (size_t i=0;i<PartitionCorners.size();i++)
        //            for (size_t j=0;j<PartitionCorners[i].size();j++)
        //            {
        //                size_t IndexV=PartitionCorners[i][j];
        //                PatchCornerPos.push_back(Mesh().vert[IndexV].P());
        //            }
        for (size_t i=0;i<Partitions.size();i++)
        {
            std::set<size_t> CurrPartV(PartitionCorners[i].begin(),
                                       PartitionCorners[i].end());
            for (size_t j=0;j<Partitions[i].size();j++)
            {
                size_t IndexF=Partitions[i][j];
                CoordType bary=(Mesh().face[IndexF].P(0)+
                                Mesh().face[IndexF].P(1)+
                                Mesh().face[IndexF].P(2))/3;
                for (size_t e=0;e<3;e++)
                {
                    size_t IndexV=vcg::tri::Index(Mesh(),Mesh().face[IndexF].V(e));
                    CoordType VertPos=Mesh().vert[IndexV].P();
                    if (CurrPartV.count(IndexV)==0)continue;
                    PatchCornerPos.push_back(bary*0.5+VertPos*0.5);
                }
            }
        }
    }


    void Init(ScalarType _Drift,bool _DebugMsg=true)
    {
        DebugMsg=true;//_DebugMsg;
        Drift=_Drift;


        InitStructures();
        InitEdgeL();
        ChoosenPaths.clear();
        //ChoosenIsLoop.clear();
        MaxNarrowWeight=sqrt(TotArea(Mesh()))*MAX_NARROW_CONST*Drift;
        UpdatePartitionsFromChoosen();
        ColorByPartitions();
        InitAvEdge();

    }

    size_t CopyPathsFrom(PatchTracer<MeshType> &Ptr,
                         std::vector<size_t> &VertMap)
    {
        std::vector<std::vector<size_t> > CurrV;
        std::vector<std::vector<size_t> > CurrDir;
        std::vector<bool> IsLoop;
        Ptr.GetCurrVertDir(CurrV,CurrDir,IsLoop);

        //get the vertices type
        assert(VertMap.size()==Mesh().vert.size());

        //initialize the inv map
        std::map<size_t,size_t> InvVertMap;
        for (size_t i=0;i<VertMap.size();i++)
            InvVertMap[VertMap[i]]=i;

        std::vector<std::vector<std::pair<size_t,size_t> > > PathNodes(CurrV.size());
        std::vector<std::vector<std::pair<size_t,size_t> > > PathDirs(CurrV.size());

        //std::cout<<"SIZE CHOSEN 0 "<<ChoosenPaths.size()<<std::endl;
        for (size_t i=0;i<CurrV.size();i++)
        {
            size_t Limit=CurrV[i].size()-1;
            if (IsLoop[i])Limit++;
            for (size_t j=0;j<Limit;j++)
            {
                size_t IndexV0=CurrV[i][j];
                size_t IndexV1=CurrV[i][(j+1)%CurrV[i].size()];
                size_t IndexD0=CurrDir[i][j];
                size_t IndexD1=CurrDir[i][(j+1)%CurrDir[i].size()];

                if (InvVertMap.count(IndexV0)==0)continue;//the vertex is not in the submesh
                if (InvVertMap.count(IndexV1)==0)continue;//the vertex is not in the submesh

                //get indexes in the submesh
                IndexV0=InvVertMap[IndexV0];
                IndexV1=InvVertMap[IndexV1];

                assert(IndexV0<Mesh().vert.size());
                assert(IndexV1<Mesh().vert.size());

                int IndexF,IndexE;
                Mesh().WichFaceEdge(IndexV0,IndexV1,IndexF,IndexE);
                if (IndexF==-1)continue;//such edge is not in the submesh

                assert(IndexE>=0);
                assert(IndexE<3);

                //if is in the border then does not add it
                if (vcg::face::IsBorder(Mesh().face[IndexF],IndexE))continue;

                std::pair<size_t,size_t> EdgeV(IndexV0,IndexV1);
                std::pair<size_t,size_t> EdgeD(IndexD0,IndexD1);

                //else add to the chosen
                PathNodes[i].push_back(EdgeV);
                PathDirs[i].push_back(EdgeD);
            }
        }

        size_t insertedPath=0;
        for (size_t i=0;i<PathNodes.size();i++)
        {
            if (PathNodes[i].size()==0)continue;

            size_t IndexV=PathNodes[i][0].first;
            assert(IndexV<Mesh().vert.size());
            size_t IndexD=PathDirs[i][0].first;
            size_t IndexN0=VFGraph.IndexNode(IndexV,IndexD);
            ChoosenPaths.push_back(CandidateTrace(TVNone,TVNone,TraceDirect,IndexN0));
            insertedPath++;
            ChoosenPaths.back().PathNodes.push_back(IndexN0);
            //            size_t test=VFGraph.NodeVertI(ChoosenPaths.back().PathNodes.back());
            //            assert(test<Mesh().vert.size());
            for (size_t j=1;j<PathNodes[i].size();j++)
            {
                size_t IndexV=PathNodes[i][j].first;
                assert(IndexV<Mesh().vert.size());
                size_t IndexD=PathDirs[i][j].first;
                size_t IndexVPred=PathNodes[i][j-1].second;
                size_t IndexDPred=PathDirs[i][j-1].second;
                //in this case start a net path
                if (IndexV!=IndexVPred)
                {
                    //add the last
                    size_t IndexVEnd=IndexVPred;
                    size_t IndexDEnd=IndexDPred;
                    assert(IndexVEnd<Mesh().vert.size());
                    size_t IndexNE=VFGraph.IndexNode(IndexVEnd,IndexDEnd);
                    ChoosenPaths.back().PathNodes.push_back(IndexNE);
                    //start a new path
                    size_t IndexN0=VFGraph.IndexNode(IndexV,IndexD);
                    ChoosenPaths.push_back(CandidateTrace(TVNone,TVNone,TraceDirect,IndexN0));
                    insertedPath++;
                }
                //adf the new element at the end
                size_t IndexN=VFGraph.IndexNode(IndexV,IndexD);
                ChoosenPaths.back().PathNodes.push_back(IndexN);
                //                test=VFGraph.NodeVertI(ChoosenPaths.back().PathNodes.back());
                //                assert(test<Mesh().vert.size());
            }
            //add the last one
            size_t IndexVEnd=PathNodes[i].back().second;
            //            std::cout<<"Vert 0 "<<IndexVEnd<<std::endl;

            assert(IndexVEnd<Mesh().vert.size());
            size_t IndexDEnd=PathDirs[i].back().second;
            size_t IndexNE=VFGraph.IndexNode(IndexVEnd,IndexDEnd);

            //            size_t test1=VFGraph.NodeVertI(IndexNE);
            //            std::cout<<"Vert 1.5 "<<test1<<std::endl;

            ChoosenPaths.back().PathNodes.push_back(IndexNE);

        }

        //        MeshType traceMesh;
        //        std::vector<bool> Selected(ChoosenPaths.size(),false);
        //        //Selected[i]=true;
        //        MeshTraces(VFGraph,ChoosenPaths,Selected,traceMesh);
        //        //vcg::tri::io::ExporterPLY<MeshType>::Save(VFGraph.Mesh(),"double_direction_domain.ply");
        //        vcg::tri::io::ExporterPLY<MeshType>::Save(traceMesh,"double_direction_error.ply",vcg::tri::io::Mask::IOM_FACECOLOR);
        //        std::cout<<"BBB"<<std::endl;

        return insertedPath;
    }

    void CopyFrom(PatchTracer<MeshType> &Ptr,
                  std::vector<size_t> &VertMap,
                  size_t IndexPatch)
    {
        Drift=Ptr.Drift;
        DebugMsg=Ptr.DebugMsg;
        MaxNarrowWeight=Ptr.MaxNarrowWeight;
        avgEdge=Ptr.avgEdge;


        //get the vertices type
        assert(VertMap.size()==Mesh().vert.size());

        //set the vert type
        VertType.resize(Mesh().vert.size(),TVNone);
        VerticesNeeds.resize(Mesh().vert.size(),0);

        //copy the original vert type
        for (size_t i=0;i<VertMap.size();i++)
        {
            size_t IndexV=VertMap[i];
            assert(IndexV<Ptr.Mesh().vert.size());
            assert(IndexV<Ptr.VertType.size());
            VertType[i]=Ptr.VertType[IndexV];
            VerticesNeeds[i]=Ptr.VerticesNeeds[IndexV];
        }


        //update coherently for the new patch
        for (size_t i=0;i<VertType.size();i++)
        {
            if ((VertType[i]!=TVNarrow)&&
                    (VertType[i]!=TVConcave)&&
                    Mesh().vert[i].IsB())
                VertType[i]=TVFlat;

            //            if (VertType[i]==Internal)
            //                VertType[i]=None;

            if ((VertType[i]==TVNarrow)&&(VerticesNeeds[i]==0))
                VertType[i]=TVFlat;

            if ((VertType[i]==TVConcave)&&(VerticesNeeds[i]==0))
                VertType[i]=TVFlat;
        }

        //then also set the corners as convex
        std::set<size_t> CornerSet(Ptr.PartitionCorners[IndexPatch].begin(),
                                   Ptr.PartitionCorners[IndexPatch].end());
        for (size_t i=0;i<VertType.size();i++)
        {
            size_t IndexV=VertMap[i];
            if (CornerSet.count(IndexV)>0)
                VertType[i]=TVConvex;
        }


        //get configuration on borders
        std::vector<size_t> FlatEmitters,FlatReceivers,ChosenEmitters,
                ChosenReceivers,FlatTangent,ChosenTangent;

        Ptr.GetPatchBorderNodes(IndexPatch,FlatEmitters,FlatReceivers,ChosenEmitters,
                                ChosenReceivers,FlatTangent,ChosenTangent);

        //initialize the emitters
        NodeEmitterTypes.clear();
        NodeReceiverTypes.clear();
        NodeEmitterTypes.resize(VFGraph.NumNodes(),TVNone);
        NodeReceiverTypes.resize(VFGraph.NumNodes(),TVNone);

        FlatEmitters.insert(FlatEmitters.end(),ChosenEmitters.begin(),ChosenEmitters.end());
        FlatReceivers.insert(FlatReceivers.end(),ChosenReceivers.begin(),ChosenReceivers.end());

        std::set<std::pair<size_t,size_t> > FlattenEmitSet;
        std::set<std::pair<size_t,size_t> > FlattenReceiveSet;

        for (size_t i=0;i<FlatEmitters.size();i++)
        {
            size_t IndexV,IndexDir;
            VertexFieldGraph<MeshType>::VertDir(FlatEmitters[i],IndexV,IndexDir);
            FlattenEmitSet.insert(std::pair<size_t,size_t>(IndexV,IndexDir));
        }
        for (size_t i=0;i<FlatReceivers.size();i++)
        {
            size_t IndexV,IndexDir;
            VertexFieldGraph<MeshType>::VertDir(FlatReceivers[i],IndexV,IndexDir);
            FlattenReceiveSet.insert(std::pair<size_t,size_t>(IndexV,IndexDir));
        }

        //then adjust concave or convex
        for (size_t i=0;i<VertType.size();i++)
        {
            if (VertType[i]==TVFlat)
            {
                size_t OrigV=VertMap[i];
                for (size_t Dir=0;Dir<4;Dir++)
                {
                    std::pair<size_t,size_t> Key(OrigV,Dir);
                    size_t IndexN=VertexFieldGraph<MeshType>::IndexNode(i,Dir);
                    if (FlattenReceiveSet.count(Key)>0)
                    {
                        NodeReceiverTypes[IndexN]=TVFlat;
                        //assert(FlattenEmitSet.count(Key)==0);
                    }else
                        if (FlattenEmitSet.count(Key)>0)
                        {
                            NodeEmitterTypes[IndexN]=TVFlat;
                            //assert(FlattenReceiveSet.count(Key)==0);
                        }
                }
            }
            else
                //in this case copy from original
                if ((VertType[i]==TVNarrow)||(VertType[i]==TVConcave))
                {
                    std::vector<size_t> Nodes0;
                    VertexFieldGraph<MeshType>::IndexNodes(i,Nodes0);

                    size_t IndexV=VertMap[i];
                    std::vector<size_t> Nodes1;
                    VertexFieldGraph<MeshType>::IndexNodes(IndexV,Nodes1);

                    assert(Nodes0.size()==4);
                    assert(Nodes1.size()==4);
                    for (size_t i=0;i<Nodes0.size();i++)
                    {
                        NodeEmitterTypes[Nodes0[i]]=Ptr.NodeEmitterTypes[Nodes1[i]];
                        NodeReceiverTypes[Nodes0[i]]=Ptr.NodeReceiverTypes[Nodes1[i]];
                    }
                }
        }
    }

    void GetPatchNodes(const size_t &IndexPatch,
                       std::vector<size_t> &PatchNodes)
    {
        for (size_t i=0;i<Partitions[IndexPatch].size();i++)
        {
            size_t FaceI=Partitions[IndexPatch][i];
            for (size_t j=0;j<3;j++)
            {
                size_t VertI=vcg::tri::Index(Mesh(),Mesh().face[FaceI].V(j));
                std::vector<size_t> NodeI;
                VertexFieldGraph<MeshType>::IndexNodes(VertI,NodeI);
                PatchNodes.insert(PatchNodes.end(),NodeI.begin(),NodeI.end());
            }
        }
        std::sort(PatchNodes.begin(),PatchNodes.end());
        auto last=std::unique(PatchNodes.begin(),PatchNodes.end());
        PatchNodes.erase(last, PatchNodes.end());
    }

    enum ChooseMode{Fartest,Shortest};

    bool TraceFrom(TypeVert FromType,
                   TypeVert ToType,
                   TraceType TrType,
                   const std::vector<bool> &CanEmit,
                   const std::vector<bool> &CanReceive,
                   const ChooseMode ChMode,
                   bool UsePartitionNeeds=false)//,
    //bool checkOnlylastConfl=false)
    {
        Candidates.clear();

        if (DebugMsg)
            std::cout<<"Adding candidates (Trace From)"<<std::endl;

        for (size_t i=0;i<VFGraph.NumNodes();i++)
        {
            //not the same kind
            if (!CanEmit[i])continue;

            //debug
            //            if (FromType==TVNarrow)
            //                std::cout<<"BOIADE "<<i<<std::endl;

            //should be active
            assert(VFGraph.IsActive(i));

            Candidates.push_back(CandidateTrace(FromType,ToType,TrType,i));
        }
        if (DebugMsg)
        {
            std::cout<<"There are "<<Candidates.size()<<" Candidates "<<std::endl;
            std::cout<<"Updating candidates"<<std::endl;
        }
        UpdateCandidates(CanReceive);

        if (DebugMsg)
            std::cout<<"Before Expansion there are "<<Candidates.size()<<" candidates"<<std::endl;

        ExpandCandidates();

        if (Candidates.size()==0)return false;

        int size0=ChoosenPaths.size();

        if (DebugMsg)
            std::cout<<"After Expansion there are "<<Candidates.size()<<" candidates"<<std::endl;

        bool UseNodeVal=((FromType==TVNarrow)||(FromType==TVConcave));
        if (ChMode==Fartest)
            ChooseGreedyByDistance(UseNodeVal,UsePartitionNeeds);
        else
        {
            //ChoosenPaths=std::vector<CandidateTrace>(Candidates.begin(),Candidates.end());
            //ChooseGreedyByLengthVertNeeds(UseNodeVal);
            ChooseGreedyByLength(UseNodeVal,UsePartitionNeeds);
        }
        //ChooseGreedyByLengthVertNeeds(false,checkOnlylastConfl);

        int size1=ChoosenPaths.size();

        if (DebugMsg)
            std::cout<<"Choosen "<<size1-size0<<std::endl;

        return  ((size1-size0)>0);
    }

    void GetPatchBorderNodes(const size_t &IndexPatch,
                             std::vector<size_t> &FlatEmitters,
                             std::vector<size_t> &FlatReceivers,
                             std::vector<size_t> &ChosenEmitters,
                             std::vector<size_t> &ChosenReceivers,
                             std::vector<size_t> &FlatTangent,
                             std::vector<size_t> &ChosenTangent)
    {
        FlatEmitters.clear();
        FlatReceivers.clear();
        ChosenEmitters.clear();
        ChosenReceivers.clear();
        FlatTangent.clear();
        ChosenTangent.clear();

        std::vector<size_t> PatchNodes;
        GetPatchNodes(IndexPatch,PatchNodes);
        std::set<size_t> PatchNodesSet(PatchNodes.begin(),PatchNodes.end());


        //        GetFlatEmitters(FlatEmitters);
        //        GetFlatReceivers(FlatReceivers);
        GetActiveEmittersType(TVFlat,FlatEmitters);
        GetActiveReceiversType(TVFlat,FlatReceivers);

        //just filters the one in the patch
        FlatEmitters=FilterFromSet(FlatEmitters,PatchNodesSet);
        FlatReceivers=FilterFromSet(FlatReceivers,PatchNodesSet);

        GetFlatTangentNodes(FlatTangent);
        FlatTangent=FilterFromSet(FlatTangent,PatchNodesSet);

        //then get the one from choosen
        std::vector<size_t> ChoosenNodes;
        //GetChoosenNodes(ChoosenNodes);
        GetCandidateNodes(ChoosenPaths,ChoosenNodes);
        ChoosenNodes=FilterFromSet(ChoosenNodes,PatchNodesSet);

        //GetChoosenNodesAndTangent(ChosenTangent);
        GetCandidateNodesNodesAndTangent<MeshType>(ChoosenPaths,ChosenTangent);
        ChosenTangent=FilterFromSet(ChosenTangent,PatchNodesSet);

        //then for each one get the ortho nodes
        for (size_t i=0;i<ChoosenNodes.size();i++)
        {
            size_t OrthoN0,OrthoN1;
            VertexFieldGraph<MeshType>::OrthoNode(ChoosenNodes[i],OrthoN0,OrthoN1);
            assert(OrthoN0!=OrthoN1);
            std::vector<size_t> NodeNeigh0,NodeNeigh1;

            VFGraph.GetNodeNeigh(OrthoN0,NodeNeigh0);
            NodeNeigh0=FilterFromSet(NodeNeigh0,PatchNodesSet);

            VFGraph.GetNodeNeigh(OrthoN1,NodeNeigh1);
            NodeNeigh1=FilterFromSet(NodeNeigh1,PatchNodesSet);

            if (NodeNeigh0.size()>NodeNeigh1.size())
            {
                ChosenEmitters.push_back(OrthoN0);
                ChosenReceivers.push_back(OrthoN1);
            }
            else
            {
                ChosenEmitters.push_back(OrthoN1);
                ChosenReceivers.push_back(OrthoN0);
            }
        }
    }

    void ColorByPartitions()
    {
        ColorMeshByPartitions(Mesh(),Partitions);
    }

    void ColorByExpValence()
    {
        ColorMeshByExpValence(Mesh(),PartitionCorners,Partitions);
    }

    void ColorByValence()
    {
        ColorMeshByValence(Mesh(),PartitionCorners,Partitions,MinVal,MaxVal);
    }

    void ColorByTopology()
    {
        ColorMeshByTopology(Mesh(),Partitions,PartitionType);
    }

    void ColorByLenghtVariance()
    {
        ColorByVarianceLenght(Mesh(),Partitions,PartitionCorners,EdgeL);
    }

    void ColorByLenghtDistortion()
    {
        ColorByDistortionLenght(Mesh(),Partitions,PartitionCorners,EdgeL);
    }

    void ColorByArapDistortion()
    {
        ColorByUVDistortionFaces(Mesh(),Partitions,PartitionCorners,Arap,false,false);
    }

    void ColorByCClarkability()
    {
        InitEdgeL();
        ColorByCatmullClarkability(Mesh(),Partitions,PartitionCorners,
                                   EdgeL,CClarkability,avgEdge);
    }

    size_t UnsatisfiedNum()
    {
        size_t UnsatisfiedNum=0;
        for (size_t i=0;i<VerticesNeeds.size();i++)
            UnsatisfiedNum+=VerticesNeeds[i];
        return UnsatisfiedNum;
    }

    void JoinNarrowStep(bool UpdatePartition=true)
    {
        bool Joined=true;
        size_t NumPath0=ChoosenPaths.size();
        do
        {
            Joined=false;

            //NARROW TO NARROW
            //std::cout<<"0a"<<std::endl;
            Joined|=JoinConnection(TVNarrow,TVNarrow,DijkstraReceivers);

            //std::cout<<"1a"<<std::endl;
            //NARROW TO CONCAVE
            Joined|=JoinConnection(TVNarrow,TVConcave,DijkstraReceivers);

            //std::cout<<"2a"<<std::endl;
            //NARROW TO FLAT
            Joined|=JoinConnection(TVNarrow,TVFlat,TraceDirect);
            Joined|=JoinConnection(TVNarrow,TVFlat,DijkstraReceivers);

            //std::cout<<"3a"<<std::endl;
            //            //NARROW TO TRACED
            //            Joined|=JoinConnection(Narrow,Choosen,TraceDirect);
            //            Joined|=JoinConnection(Narrow,Choosen,DijkstraReceivers);

            if (DebugMsg)
                std::cout<<"Still "<<UnsatisfiedNum()<<" Non Connected"<<std::endl;
        }
        while (Joined);
        size_t NumPath1=ChoosenPaths.size();
        if (NumPath1==NumPath0)return;

        if(UpdatePartition)
        {
            UpdatePartitionsFromChoosen();
            ColorByPartitions();
        }
    }

    void JoinNarrow(bool UpdatePartition=true)
    {
//        MaxNarrowWeight/=100;
//        JoinNarrowStep(UpdatePartition);
//        MaxNarrowWeight*=100;
        JoinNarrowStep(UpdatePartition);
    }

    void JoinConcaveStep(bool UpdatePartition=true)
    {
        bool Joined=true;
        size_t NumPath0=ChoosenPaths.size();
        do
        {
            Joined=false;
            //
            //CONCAVE TO CONCAVE
            Joined|=JoinConnection(TVConcave,TVConcave,DijkstraReceivers);

            //CONCAVE TO FLAT
            Joined|=JoinConnection(TVConcave,TVFlat,TraceDirect);
            Joined|=JoinConnection(TVConcave,TVFlat,DijkstraReceivers);

            //            //CONCAVE TO TRACED
            //            Joined|=JoinConnection(Concave,Choosen,TraceDirect);
            //            Joined|=JoinConnection(Concave,Choosen,DijkstraReceivers);
            if (DebugMsg)
                std::cout<<"Still "<<UnsatisfiedNum()<<" Non Connected"<<std::endl;
        }
        while (Joined);
        //JoinConnection(Concave,Flat,DijkstraReceivers);

        size_t NumPath1=ChoosenPaths.size();
        if (NumPath1==NumPath0)return;

        if(UpdatePartition)
        {
            UpdatePartitionsFromChoosen();
            ColorByPartitions();
        }
    }

    void JoinConcave(bool UpdatePartition=true)
    {
//        MaxNarrowWeight/=100;
//        JoinConcaveStep(UpdatePartition);
//        MaxNarrowWeight*=100;
        JoinConcaveStep(UpdatePartition);
    }

    void JoinBoundaries(bool UpdatePartition=true,
                        bool UsePartitionNeeds=false)
    {
        //        std::vector<bool> IsActive;
        //        VFGraph.IsActiveNodes(IsActive);

        //        bool Joined=true;
        //        size_t NumPath0=ChoosenPaths.size();
        //        do
        //        {
        //            Joined=false;
        //            Joined|=JoinConnection(Flat,Flat,DijkstraReceivers);
        //            //Joined|=JoinConnection(Flat,Choosen,TraceDirect);
        //        }while (Joined);


        //        VFGraph.SetActiveNodes(IsActive);

        //        size_t NumPath1=ChoosenPaths.size();
        //        if (NumPath1==NumPath0)return;

        if (DebugMsg)
            std::cout<<"**TRACING BORDERS ***"<<std::endl;

        bool Joined=true;
        do{
            Joined=false;
            Joined|=JoinConnection(TVFlat,TVFlat,TraceDirect,UsePartitionNeeds);
            Joined|=JoinConnection(TVFlat,TVFlat,DijkstraReceivers,UsePartitionNeeds);
            //            InitCandidates(Flat,Flat,TraceDirect);
            //            InitCandidates(Flat,Flat,DijkstraReceivers);

            //            ChooseGreedyByDistance(false,UsePartitionNeeds);
        }while (Joined);

        if(UpdatePartition)
        {
            UpdatePartitionsFromChoosen();
            ColorByPartitions();
        }
    }


    void TraceLoops(bool UpdatePartition=true, bool UsePartitionNeeds=false)
    {

        if (DebugMsg)
            std::cout<<"**TRACING LOOPS ***"<<std::endl;

        InitCandidates(TVInternal,TVInternal,TraceLoop);
        ChooseGreedyByDistance(false,UsePartitionNeeds);
        if(UpdatePartition)
        {
            UpdatePartitionsFromChoosen();
            ColorByPartitions();
        }
    }

    void BatchRemoval(bool do_smooth=true)
    {
        if (do_smooth)
        {
            UpdatePartitionsFromChoosen(false);
            SmoothPatches(2);
        }
        RemovePaths();//false);
        if (!split_on_removal){
            //FixValences();
            if (DebugMsg)
                WriteInfo();
            return;
        }else
        {
            SplitIntoSubPaths();
            RemovePaths();//false);
            //FixValences();
        }
        if (DebugMsg)
            WriteInfo();
    }


    void FixCorners(size_t IndexPatch)
    {

        //save the original vertex index
        for (size_t i=0;i<Mesh().vert.size();i++)
            Mesh().vert[i].Q()=i;

        //extract the submesh
        vcg::tri::UpdateSelection<MeshType>::FaceClear(Mesh());
        vcg::tri::UpdateSelection<MeshType>::VertexClear(Mesh());

        for (size_t i=0;i<Partitions[IndexPatch].size();i++)
            Mesh().face[Partitions[IndexPatch][i]].SetS();

        vcg::tri::UpdateSelection<MeshType>::VertexFromFaceLoose(Mesh());
        MeshType TestMesh;
        vcg::tri::Append<MeshType,MeshType>::Mesh(TestMesh,Mesh(),true);
        TestMesh.UpdateAttributes();

        std::set<size_t> Corners(PartitionCorners[IndexPatch].begin(),
                                 PartitionCorners[IndexPatch].end());

        //compute sum angle
        std::vector<ScalarType > angleSumH(TestMesh.vert.size(),0);

        for(size_t i=0;i<TestMesh.face.size();i++)
        {
            for(int j=0;j<TestMesh.face[i].VN();++j)
            {
                size_t IndexV=vcg::tri::Index(TestMesh,TestMesh.face[i].V(j));
                angleSumH[IndexV] += face::WedgeAngleRad(TestMesh.face[i],j);
            }
        }

        std::vector<std::pair<ScalarType,size_t> > AngleVert;

        for(size_t i=0;i<TestMesh.vert.size();i++)
        {
            if(!TestMesh.vert[i].IsB())continue;
            AngleVert.push_back(std::pair<ScalarType,size_t>(angleSumH[i],TestMesh.vert[i].Q()));
        }
        std::sort(AngleVert.begin(),AngleVert.end());
        assert(AngleVert.size()>=3);

        if (Corners.size()<MinVal)
        {
            for (size_t i=0;i<AngleVert.size();i++)
            {
                if (Corners.count(AngleVert[i].second)>0)continue;
                Corners.insert(AngleVert[i].second);
                if (Corners.size()==MinVal)break;
            }
            assert(Corners.size()==MinVal);
        }
        else
        {
            assert(Corners.size()>MaxVal);
            std::reverse(AngleVert.begin(),AngleVert.end());
            for (size_t i=0;i<AngleVert.size();i++)
            {
                int IndexV=AngleVert[i].second;
                if (Corners.count(IndexV)==0)continue;
                if ((VertType[IndexV]==TVNarrow)
                        ||(VertType[IndexV]==TVConcave)
                        ||(VertType[IndexV]==TVConvex))continue;
                Corners.erase(IndexV);
                if (Corners.size()==MaxVal)break;
            }
            if (Corners.size()>MaxVal)
            {
                for (size_t i=0;i<AngleVert.size();i++)
                {
                    int IndexV=AngleVert[i].second;
                    if (Corners.count(IndexV)==0)continue;
                    Corners.erase(IndexV);
                    if (Corners.size()==MaxVal)break;
                }
            }
            assert(Corners.size()==MaxVal);
        }

        PartitionCorners[IndexPatch]=std::vector<size_t>(Corners.begin(),Corners.end());
    }

    void FixValences()
    {
        size_t NeedFix=0;
        for (size_t i=0;i<PartitionCorners.size();i++)
            if ((PartitionCorners[i].size()<MIN_ADMITTIBLE)||
                    (PartitionCorners[i].size()>MAX_ADMITTIBLE))
            {
                NeedFix++;
                FixCorners(i);
            }
        //if (DebugMsg)
        std::cout<<"FINAL Fixed "<<NeedFix<<std::endl;
    }

    void BatchAddLoops(bool ForceReceivers,
                       bool AddOnlyNeeded,
                       bool InterleaveRemoval,
                       bool FinalRemoval,
                       bool SmoothOnRemoval=true)
    {
        if (ForceReceivers)
        {
            AllReceivers=true;
            MaxNarrowWeight/=100;
            JoinNarrow(false);
            JoinConcave(false);
            MaxNarrowWeight*=100;
            AllReceivers=false;
        }
        JoinNarrow(false);
        JoinConcave(false);

        //then update the partitions
        if (AddOnlyNeeded)
            UpdatePartitionsFromChoosen();

        Candidates.clear();

        //        ScalarType val0=max_lenght_distortion;
        //        ScalarType val1=max_lenght_variance;

        //        max_lenght_distortion=-1;
        //        max_lenght_variance=-1;

        if (DebugMsg)
            std::cout<<"**TRACING LOOPS ***"<<std::endl;

        TraceLoops(false,AddOnlyNeeded);

        if (InterleaveRemoval)
            BatchRemoval(SmoothOnRemoval);

        if (DebugMsg)
            std::cout<<"Num Chosen "<<ChoosenPaths.size()<<std::endl;

        //        max_lenght_distortion=val0;
        //        max_lenght_variance=val1;

        if (AddOnlyNeeded)
            UpdatePartitionsFromChoosen();

        if (DebugMsg)
            std::cout<<"**TRACING BORDERS ***"<<std::endl;

        JoinBoundaries(false,AddOnlyNeeded);

        if (FinalRemoval)
            BatchRemoval(SmoothOnRemoval);
    }

    void InitCandidates(TypeVert FromType,
                        TypeVert ToType,
                        TraceType TrType)
    {

        std::vector<bool> CanEmit,CanReceive,MustDisable;
        GetTracingConfiguration(FromType,ToType,TrType,CanEmit,CanReceive,MustDisable);

        VFGraph.SetAllActive();
        VFGraph.SetDisabledNodes(MustDisable);

        if (DebugMsg)
            std::cout<<"Adding candidates (Init)"<<std::endl;

        for (size_t i=0;i<VFGraph.NumNodes();i++)
        {
            //not the same kind
            if (!CanEmit[i])continue;

            //should be active
            assert(VFGraph.IsActive(i));

            Candidates.push_back(CandidateTrace(FromType,ToType,TrType,i));
        }
        if (DebugMsg)
            std::cout<<"There are "<<Candidates.size()<<" Candidates "<<std::endl;
        UpdateCandidates(CanReceive);

        if (DebugMsg)
            std::cout<<"Before Expansion there are "<<Candidates.size()<<" candidates"<<std::endl;

        ExpandCandidates();

        if (DebugMsg)
            std::cout<<"After Expansion there are "<<Candidates.size()<<" candidates"<<std::endl;

    }

    void  GetCurrCandidates(std::vector<std::vector<size_t> > &CurrCandidates)
    {
        CurrCandidates.clear();
        for (size_t i=0;i<Candidates.size();i++)
            CurrCandidates.push_back(Candidates[i].PathNodes);
    }

    void  GetCurrCandidatesIsLoop(std::vector<bool> &CurrCandidatesIsLoop)
    {
        CurrCandidatesIsLoop.clear();
        for (size_t i=0;i<Candidates.size();i++)
            CurrCandidatesIsLoop.push_back(Candidates[i].IsLoop);
    }

    void  GetCurrChosen(std::vector<std::vector<size_t> > &CurrChosen)
    {
        CurrChosen.clear();
        for (size_t i=0;i<ChoosenPaths.size();i++)
            CurrChosen.push_back(ChoosenPaths[i].PathNodes);
    }

    void  GetCurrVertDir(std::vector<std::vector<size_t> > &CurrV,
                         std::vector<std::vector<size_t> > &CurrDir,
                         std::vector<bool> &IsLoop)
    {
        CurrV.clear();
        CurrDir.clear();
        IsLoop.clear();
        CurrV.resize(ChoosenPaths.size());
        CurrDir.resize(ChoosenPaths.size());
        for (size_t i=0;i<ChoosenPaths.size();i++)
        {
            VFGraph.NodeVertI(ChoosenPaths[i].PathNodes,CurrV[i]);
            VFGraph.NodeDirI(ChoosenPaths[i].PathNodes,CurrDir[i]);
            IsLoop.push_back(ChoosenPaths[i].IsLoop);
        }
    }

    void GetUnsatisfied(std::vector<size_t> &Remaining)
    {
        Remaining.clear();
        for (size_t i=0;i<VerticesNeeds.size();i++)
        {
            if (VerticesNeeds[i]==0)continue;
            //assert((VertType[i]==Narrow)||(VertType[i]==Concave));
            std::vector<size_t> NodesI;
            VertexFieldGraph<MeshType>::IndexNodes(i,NodesI);
            for (size_t j=0;j<NodesI.size();j++)
            {
                if (!VFGraph.IsActive(NodesI[j]))continue;
                if((NodeEmitterTypes[NodesI[j]]==TVNarrow)||
                        (NodeEmitterTypes[NodesI[j]]==TVConcave))
                {
                    Remaining.push_back(NodesI[j]);
                }
            }
        }
    }

    void GetUnsolvedPartitions(std::vector<std::vector<size_t> > &UnsolvedPartition,
                               std::vector<PatchType> &UnsolvedType)
    {
        UnsolvedPartition.clear();
        UnsolvedType.clear();

        UpdatePartitionsFromChoosen(true);
        for (size_t i=0;i<Partitions.size();i++)
        {
            if (PartitionType[i]==IsOK)continue;
            UnsolvedPartition.push_back(Partitions[i]);
            UnsolvedType.push_back(PartitionType[i]);
        }
    }

    void GetCurrChosenIsLoop(std::vector<bool> &ChosenIsLoop)
    {
        ChosenIsLoop.clear();
        GetIsLoop(ChoosenPaths,ChosenIsLoop);
    }

    void GetUnsolvedPartitionsIndex(std::vector<size_t > &UnsolvedPartitionIndex,
                                    std::vector<PatchType> &PatchTypes)
    {
        UnsolvedPartitionIndex.clear();
        PatchTypes.clear();

        for (size_t i=0;i<Partitions.size();i++)
        {
            if (PartitionType[i]==IsOK)continue;
            UnsolvedPartitionIndex.push_back(i);
            PatchTypes.push_back(PartitionType[i]);
        }
    }

    void SetChoosenFromVertDir(const std::vector<std::vector<size_t> > &VertIdx,
                               const std::vector<std::vector<size_t> > &VertDir,
                               const std::vector<bool> &IsLoop)
    {
        ChoosenPaths.clear();
        assert(VertIdx.size()==VertDir.size());
        assert(VertIdx.size()==IsLoop.size());

        for (size_t i=0;i<VertIdx.size();i++)
        {
            CandidateTrace CTrace;
            for (size_t j=0;j<VertIdx[i].size();j++)
            {
                size_t IndexV=VertIdx[i][j];
                size_t IndexDir=VertDir[i][j];
                size_t IndexN=VertexFieldGraph<MeshType>::IndexNode(IndexV,IndexDir);
                CTrace.PathNodes.push_back(IndexN);
            }

            CTrace.FromType=VertType[VertIdx[i][0]];
            CTrace.ToType=VertType[VertIdx[i].back()];
            CTrace.TracingMethod=TraceDirect;
            CTrace.InitNode=CTrace.PathNodes[0];
            CTrace.IsLoop=IsLoop[i];
            CTrace.Updated=true;
            ChoosenPaths.push_back(CTrace);
            assert(ChoosenPaths.back().PathNodes.size()>=2);
        }
        UpdateVertNeedsFromChoosen();
        UpdatePartitionsFromChoosen();
        ColorByPartitions();
        InitEdgeL();
    }

    void getCornerSharp(std::vector<size_t> &CornerSharp)
    {
        CornerSharp.clear();
        for (size_t i=0;i<VertType.size();i++)
        {
            if ((VertType[i]==TVNarrow)||
                    (VertType[i]==TVConcave)||
                    (VertType[i]==TVConvex))
                CornerSharp.push_back(i);
        }
    }

    void TestGetNodes(const TypeVert FromType,
                      const TypeVert ToType,
                      const TraceType TracingType,
                      std::vector<size_t> &Emitter,
                      std::vector<size_t> &Receiver,
                      std::vector<size_t> &Disabled)
    {
        Emitter.clear();
        Receiver.clear();
        Disabled.clear();
        std::vector<bool> CanEmit,CanReceive,MustDisable;
        GetTracingConfiguration(FromType,ToType,TracingType,CanEmit,CanReceive,MustDisable);

        for (size_t i=0;i<CanEmit.size();i++)
            if (CanEmit[i])Emitter.push_back(i);
        for (size_t i=0;i<CanReceive.size();i++)
            if (CanReceive[i])Receiver.push_back(i);
        for (size_t i=0;i<MustDisable.size();i++)
            if (MustDisable[i])Disabled.push_back(i);
    }

    struct PathInfo
    {
        int PatchNum;
        int ValNum[10];
    };

    void ComputePatchesUV()
    {
        MeshType splittedUV;
        ParametrizePatches(Mesh(),splittedUV,Partitions,PartitionCorners,Arap,false,true,true,true);
    }

    void SubdivideIrrPatches()
    {
        std::vector<bool> MustSplit(Partitions.size(),false);
        for (size_t i=0;i<PartitionCorners.size();i++)
            if (PartitionCorners[i].size()!=4)
                MustSplit[i]=true;

        std::vector<std::vector<ScalarType> > SplitSide(PartitionCorners.size());
        for (size_t i=0;i<SplitSide.size();i++)
            SplitSide[i].resize(PartitionCorners[i].size(),0.5);


        PatchSplitter<MeshType> PSplit(Mesh());
        std::vector<std::vector<size_t> > NewFacePaches,NewCorners;
        PSplit.Subdivide(Partitions,PartitionCorners,MustSplit,SplitSide,NewFacePaches,NewCorners);

        Partitions=NewFacePaches;
        PartitionCorners=NewCorners;

        DerivePerFacePartition(Mesh(),Partitions,FacePartitions);

        //SelectMeshPatchBorders(Mesh(),FacePartitions);
        Mesh().SelectSharpFeatures();

        ColorByPartitions();
        ChoosenPaths.clear();

        assert(PartitionCorners.size()==Partitions.size());

        //        MeshType testMesh;
        //        for (size_t i=0;i<PartitionCorners.size();i++)
        //            for (size_t j=0;j<PartitionCorners[i].size();j++)
        //                vcg::tri::Allocator<MeshType>::AddVertex(testMesh,Mesh().vert[PartitionCorners[i][j]].P());

        //        vcg::tri::io::ExporterPLY<MeshType>::Save(testMesh,"test_corners.ply",vcg::tri::io::Mask::IOM_FACECOLOR);

    }

    PatchTracer(VertexFieldGraph<MeshType> &_VFGraph):VFGraph(_VFGraph)
    {
        split_on_removal=true;
        //avoid_increase_valence=true;
        //avoid_collapse_irregular=false;
        away_from_singular=true;
        match_valence=true;
        //        max_lenght_distortion=-1;//1.2;
        //        max_lenght_variance=-1;//2;
        CClarkability=1;
        sample_ratio=0.1;
        MinVal=3;
        MaxVal=5;
        AllReceivers=false;
        Concave_Need=1;
        //max_patch_area=MeshArea(Mesh())*0.5;
        //TraceLoopsBorders=true;
    }
};

#endif
