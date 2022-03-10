#ifndef PATCH_TRACER
#define PATCH_TRACER

size_t Time_FirstTrace=0;

size_t Time_InitSubPatches0=0;
size_t Time_InitSubPatches1=0;
size_t Time_InitSubPatches2=0;

//size_t Time_InitSubPatches2_0=0;
//size_t Time_InitSubPatches2_1=0;
//size_t Time_InitSubPatches2_2=0;
//size_t Time_InitSubPatches2_3=0;
//size_t Time_InitSubPatches2_4=0;
//size_t Time_InitSubPatches2_5=0;
//size_t Time_InitSubPatches2_6=0;
//size_t Time_InitSubPatches2_7=0;

size_t Time_InitSubPatches3=0;
size_t Time_TraceSubPatches=0;
size_t Time_UpdatePartitionsTotal=0;
size_t Time_UpdatePartitionsLazy=0;
size_t Time_Removal=0;
size_t Time_InitMetaMesh=0;
size_t Time_Collapse_Step0=0;
size_t Time_Collapse_Step1=0;
size_t Time_Collapse_Step2=0;
size_t Time_Collapse_Step3=0;
size_t Time_Collapse_Step4=0;

//#include "GL_vert_field_graph.h"
#include "vertex_emitter.h"
#include "patch_manager.h"
#include "metamesh.h"


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
#include <vcg/complex/algorithms/update/selection.h>
#include "vertex_emitter.h"
#include "vertex_classifier.h"
#include "candidate_path.h"
#include "patch_manager.h"
#include "patch_subdivider.h"
#include "edge_direction_table.h"
//#include <vcg/complex/algorithms/crease_cut.h>


//namespace std {

////template <typename ScalarType>
////struct hash<typename vcg::Point3<ScalarType> > {
////    std::size_t operator()(typename vcg::Point3<ScalarType> &p) const {
////        return std::hash<ScalarType>()(p.X()) ^ std::hash<ScalarType>()(p.Y()) ^ std::hash<ScalarType>()(p.Z());
////    }
////};

////   struct hash<vcg::Point3f> {
////        std::size_t operator()(const vcg::Point3f &p) const {
////            return std::hash<float>()(p.X()) ^ std::hash<float>()(p.Y()) ^ std::hash<float>()(p.Z());
////        }
////    };

////template<typename X, typename Y>
////struct hash<std::pair<X, Y> > {
////    std::size_t operator()(const std::pair<X, Y> &pair) const {
////        const size_t _HASH_P0 = 73856093u;
////        const size_t _HASH_P1 = 19349663u;
////        //const size_t _HASH_P2 = 83492791u;

////        //return size_t(p.V(0))*_HASH_P0 ^  size_t(p.V(1))*_HASH_P1 ^  size_t(p.V(2))*_HASH_P2;
////        return std::hash<X>()(pair.first)*_HASH_P0 ^ std::hash<Y>()(pair.second)*_HASH_P1;
////    }
////};

////template<typename X>
////struct hash<X> {
////    std::size_t operator()(const EdgeVert &EV) const {
////        const size_t _HASH_P0 = 73856093u;
////        const size_t _HASH_P1 = 19349663u;
////        const size_t _HASH_P2 = 83492791u;

////        //return size_t(p.V(0))*_HASH_P0 ^  size_t(p.V(1))*_HASH_P1 ^  size_t(p.V(2))*_HASH_P2;
////        return std::hash<X>()(EV.EV0)*_HASH_P0 ^ std::hash<Y>()(EV.EV1)*_HASH_P1^ std::hash<Y>()(EV.EV2)*_HASH_P2;
////    }
////};

//}

//struct EdgeVertKeyHasher
//{
//  std::size_t operator()(const EdgeVert& k) const
//  {
//    using std::size_t;
//    using std::hash;
//    using std::string;

//    return ((hash<string>()(k.EV0)
//             ^ (hash<string>()(k.EV1) << 1)) >> 1)
//             ^ (hash<int>()(k.CurrV) << 1);
//  }
//};

//enum TraceType{TraceDirect,DijkstraReceivers,TraceLoop};
enum PatchType{LowCorners,HighCorners,NonDisk,
               HasEmitter,MaxCClarkability,NonMatchValence,
               NoQualityMatch,IsOK};//MoreSing,IsOK};

enum PriorityMode{PrioModeLoop,PrioModBorder,PrioModBlend};

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

template <class FaceType>
bool IsSelectedMeshPos(const vcg::face::Pos<FaceType> &Pos)
{
    FaceType *f=Pos.F();
    int e=Pos.E();
    return(f->IsFaceEdgeS(e));
}

template <class FaceType>
void SelectMeshPos(vcg::face::Pos<FaceType> &Pos)
{
    FaceType *f=Pos.F();
    int e=Pos.E();
    f->SetFaceEdgeS(e);
    f->SetF(e);

    Pos.FlipF();

    f=Pos.F();
    e=Pos.E();
    f->SetFaceEdgeS(e);
    f->SetF(e);
}



template <class MeshType>
void RetrievePosSeqFromSelEdges(MeshType &mesh,std::vector<std::vector<vcg::face::Pos<typename MeshType::FaceType> > > &PosSeq)
{
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef vcg::face::Pos<FaceType> PosType;
    PosSeq.clear();

    //first select the ones with multiple selection or border
    vcg::tri::UpdateQuality<MeshType>::VertexConstant(mesh,0);
    for (size_t i=0;i<mesh.face.size();i++)
        for (size_t j=0;j<3;j++)
        {
            if (vcg::face::IsBorder(mesh.face[i],j))continue;
            if (!mesh.face[i].IsFaceEdgeS(j))continue;

            VertexType *v0=mesh.face[i].V0(j);
            VertexType *v1=mesh.face[i].V1(j);
            if (v0>v1)continue;

            v0->Q()+=1;
            v1->Q()+=1;
        }

    vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);
    for (size_t i=0;i<mesh.vert.size();i++)
    {
        if (mesh.vert[i].Q()==0)continue;
        if (mesh.vert[i].Q()==2)continue;
        mesh.vert[i].SetS();
    }

    //then retrieve the sequences
    std::set<std::pair<VertexType*,VertexType*> > ExploredE;
    for (size_t i=0;i<mesh.face.size();i++)
        for (size_t j=0;j<3;j++)
        {
            if (vcg::face::IsBorder(mesh.face[i],j))continue;
            if (!mesh.face[i].IsFaceEdgeS(j))continue;
            if (!mesh.face[i].V(j)->IsS())continue;

            VertexType *v0=mesh.face[i].V0(j);
            VertexType *v1=mesh.face[i].V1(j);
            std::pair<VertexType*,VertexType*> key(std::min(v0,v1),std::max(v0,v1));
            if (ExploredE.count(key)>0)continue;
            ExploredE.insert(key);

            PosType CurrPos(&mesh.face[i],j);
            CurrPos.FlipV();

            PosSeq.resize(PosSeq.size()+1);
            bool has_complete=false;
            do{
                has_complete=CurrPos.V()->IsS();
                PosSeq.back().push_back(CurrPos);

                VertexType *v0=CurrPos.V();
                VertexType *v1=CurrPos.VFlip();
                std::pair<VertexType*,VertexType*> key(std::min(v0,v1),std::max(v0,v1));
                ExploredE.insert(key);

                if (!has_complete)
                    CurrPos.NextEdgeS();

            }while (!has_complete);
        }
    std::cout<<"Retrieved "<<PosSeq.size()<<" paths"<<std::endl;
}

template <class MeshType>
bool SelectMeshPatchBorders(const VertexFieldGraph<MeshType> &VFGraph,
                            const std::vector<CandidateTrace> &TraceSet)
{
    typedef typename MeshType::FaceType FaceType;
    typedef typename vcg::face::Pos<FaceType> PosType;
    //std::set<std::pair<size_t,size_t> > BorderPatches;
    //int t0=clock();
    vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(VFGraph.Mesh());
    vcg::tri::UpdateFlags<MeshType>::FaceClearF(VFGraph.Mesh());

    //int t1=clock();
    //first add borders
    for (size_t i=0;i<VFGraph.Mesh().face.size();i++)
        for (size_t j=0;j<3;j++)
        {
            if (!VFGraph.Mesh().face[i].IsB(j))continue;
            VFGraph.Mesh().face[i].SetFaceEdgeS(j);
            VFGraph.Mesh().face[i].SetF(j);
        }
    //int t2=clock();
    //std::unordered_set<std::pair<size_t,size_t> > BorderEdges;
    for (size_t i=0;i<TraceSet.size();i++)
    {
        if (TraceSet[i].PathNodes.size()==0)continue;
        std::vector<PosType> NodePos;
        VFGraph.GetNodesPos(TraceSet[i].PathNodes,TraceSet[i].IsLoop,NodePos);
        for (size_t j=0;j<NodePos.size();j++)
        {
            if (IsSelectedMeshPos(NodePos[j]))return false;//conflict collinear segs
            SelectMeshPos(NodePos[j]);
        }
    }
    //int t3=clock();
    //        std::cout<<"Update "<<std::endl;
    //        std::cout<<"T0: "<<t1-t0<<std::endl;
    //        std::cout<<"T1: "<<t2-t1<<std::endl;
    //        std::cout<<"T2: "<<t3-t2<<std::endl;


    return true;
}


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

//int timeCopyMap=0;

template <class MeshType>
void GetEdgeDirVertMap(const VertexFieldGraph<MeshType> &VFGraph,
                       const std::vector<CandidateTrace> &TraceSet,
                       //std::unordered_map<EdgeVert,size_t,EdgeVertKeyHasher> &EdgeDirVert)
                       std::map<EdgeVert,size_t> &EdgeDirVert)
{
    //EdgeDirVert.clear();
    //int t0=clock();
    EdgeDirVert=VFGraph.EdgeBorderDir;
    //int t1=clock();
    //timeCopyMap+=t1-t0;
    //for each edge set the direction per vertex
    //std::cout<<"Num added Traces"<<TraceSet.size()<<std::endl;
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

            //            if (EdgeDirVert.count(EdgeKey0)>0)
            //            {
            //                std::cout<<"WARNING DOUBLE EDGE"<<std::endl;
            //                MeshType traceMesh;
            //                std::vector<bool> Selected(TraceSet.size(),false);
            //                Selected[i]=true;
            //                MeshTraces(VFGraph,TraceSet,Selected,traceMesh);
            //                vcg::tri::io::ExporterPLY<MeshType>::Save(VFGraph.Mesh(),"double_direction_domain.ply",vcg::tri::io::Mask::IOM_FACECOLOR);
            //                vcg::tri::io::ExporterPLY<MeshType>::Save(traceMesh,"double_direction_error.ply",vcg::tri::io::Mask::IOM_FACECOLOR);
            //                assert(0);
            //            }

            //            if (EdgeDirVert.count(EdgeKey1)>0)
            //            {
            //                std::cout<<"WARNING DOUBLE EDGE"<<std::endl;
            //                MeshType traceMesh;
            //                std::vector<bool> Selected(TraceSet.size(),false);
            //                Selected[i]=true;
            //                MeshTraces(VFGraph,TraceSet,Selected,traceMesh);
            //                vcg::tri::io::ExporterPLY<MeshType>::Save(VFGraph.Mesh(),"double_direction_domain.ply",vcg::tri::io::Mask::IOM_FACECOLOR);
            //                vcg::tri::io::ExporterPLY<MeshType>::Save(traceMesh,"double_direction_error.ply",vcg::tri::io::Mask::IOM_FACECOLOR);
            //                assert(0);
            //            }
            assert(EdgeDirVert.count(EdgeKey0)==0);
            assert(EdgeDirVert.count(EdgeKey1)==0);

            EdgeDirVert[EdgeKey0]=DirV0;
            EdgeDirVert[EdgeKey1]=((DirV1+2)%4);//put the inverse cause look internally the interval

        }
        //int t2=clock();
        //                std::cout<<"Edge Map"<<std::endl;
        //                std::cout<<"T0: "<<t1-t0<<std::endl;
        //                std::cout<<"T1: "<<t2-t1<<std::endl;
    }

}

template <class MeshType>
void FindPerVertDirs(const VertexFieldGraph<MeshType> &VFGraph,
                     const std::vector<size_t> &Partition,
                     std::map<EdgeVert,size_t> &EdgeDirVert,
                     //std::unordered_map<EdgeVert,size_t,EdgeVertKeyHasher> &EdgeDirVert,
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
            if (EdgeDirVert.count(EdgeKey0)>0)//||(VFGraph.EdgeBorderDir.count(EdgeKey0)>0))
            {
                size_t EdgeDir=EdgeDirVert[EdgeKey0];
                DirVert[IndexV0].push_back(EdgeDir);
            }
            ////            else
            ////            {
            //                //check also on borders
            //                if (VFGraph.EdgeBorderDir.count(EdgeKey0)>0)
            //                {
            //                    size_t EdgeDir=EdgeDirVert[EdgeKey0];
            //                    DirVert[IndexV0].push_back(EdgeDir);
            //                }
            ////            }

            EdgeVert EdgeKey1(MinV,MaxV,IndexV1);
            if (EdgeDirVert.count(EdgeKey1)>0)//||(VFGraph.EdgeBorderDir.count(EdgeKey1)>0))
            {
                size_t EdgeDir=EdgeDirVert[EdgeKey1];
                DirVert[IndexV1].push_back(EdgeDir);
            }
            ////            else
            ////            {
            //                //check also on borders
            //                if (VFGraph.EdgeBorderDir.count(EdgeKey1)>0)
            //                {
            //                    size_t EdgeDir=EdgeDirVert[EdgeKey1];
            //                    DirVert[IndexV1].push_back(EdgeDir);
            //                }
            //            }
        }
    }
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
        //        std::cout<<"A:"<<IndexF0<<std::endl;
        //        std::cout<<"B:"<<FacePartitions.size()<<std::endl;
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
        bool OnCorner=false;
        int ExpVal=PatchManager<MeshType>::ExpectedValence(mesh,Partitions[i],PartitionCorners[i],OnCorner);
        vcg::Color4b CurrCol;
        if (ExpVal==-1)
            CurrCol=vcg::Color4b::Gray;
        else
        {
            //            if (PartitionCorners[i].size()==ExpVal)
            //                CurrCol=vcg::Color4b::Green;
            if (((int)PartitionCorners[i].size()!=ExpVal)&&(!OnCorner))
                CurrCol=vcg::Color4b::Red;
            else
                CurrCol=vcg::Color4b::Green;
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
        vcg::Color4b CurrCol=vcg::Color4b::LightGreen;

        //if (PartitionCorners[i].size()==MinVal)
        if (PartitionCorners[i].size()==4)
            CurrCol=vcg::Color4b::Gray;

        if (PartitionCorners[i].size()==3)
            CurrCol=vcg::Color4b::Yellow;

        if (PartitionCorners[i].size()==5)
            CurrCol=vcg::Color4b::Blue;

        //if (PartitionCorners[i].size()==MaxVal)
        if (PartitionCorners[i].size()==6)
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
        case NonMatchValence:CurrCol=vcg::Color4b::DarkBlue;break;
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
//,bool DebugMsg=false)
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
    if (SelfInt)
        return false;


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
class EmptyMeshQuality
{
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::ScalarType  ScalarType;
    typedef typename MeshType::FacePointer FacePointer;

public:
    EmptyMeshQuality(){}

    ScalarType operator()(MeshType &m) const
    {(void)m;return 0;}
};

template <class TriMeshType,class PatchQualityFunctor=EmptyMeshQuality<TriMeshType> >
class PatchTracer
{
public:
    typedef PatchTracer<TriMeshType,PatchQualityFunctor> MyTracerType;
    typedef TriMeshType MeshType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename vcg::face::Pos<FaceType> PosType;

private:

    VertexFieldGraph<MeshType> &VFGraph;

public:
    //vertices types and ortho directions, used to select sources
    ScalarType Drift;

    bool split_on_removal;
    bool DebugMsg;
    //bool FirstBorder;
    PriorityMode PrioMode;
    //bool avoid_increase_valence;
    //bool avoid_collapse_irregular;
    bool away_from_singular;
    //    ScalarType max_lenght_distortion;
    //    ScalarType max_lenght_variance;
    ScalarType sample_ratio;
    ScalarType CClarkability;
    bool match_valence;
    bool check_quality_functor;
    //    ScalarType maxQThr;
    //    ScalarType minQThr;
    //ScalarType max_patch_area;
    size_t MinVal;
    size_t MaxVal;
    size_t Concave_Need;
    size_t subInt;
    bool AllowDarts;
    bool CheckUVIntersection;
    bool ContinuosCheckUVInt;
    bool CheckTJunction;
    bool AllowSelfGluedPatch;
    bool CheckQuadrangulationLimits;
    bool AllowRemoveConcave;
    //bool TraceLoopsBorders;

    std::vector<std::vector<size_t> > Partitions;
    std::vector<int> FacePartitions;
    std::vector<PatchType> PartitionType;
    std::vector<std::vector<size_t> > PartitionCorners;

    std::vector<PatchInfo<ScalarType> > PatchInfos;


    MetaMesh<MeshType> MMesh;

private:

    std::vector<TypeVert> VertType;
    std::vector<TypeVert > NodeEmitterTypes;
    std::vector<TypeVert> NodeReceiverTypes;

    bool AllReceivers;
    std::vector<ScalarType> CurrNodeDist;

    //    std::vector<std::pair<ScalarType,size_t> > CandidatesPathLenghts;
    //    std::vector<std::pair<ScalarType,size_t> > CandidatesPathDistances;

    std::vector<bool> Traceable;

    //DATA STRUCTURES FOR THE PATH AND THE NEEDS FOR EACH VERTEX
    std::vector<CandidateTrace> Candidates;

    std::vector<CandidateTrace> DiscardedCandidates;

    std::vector<ScalarType> PrioVect;

public:

    ScalarType MaxNarrowWeight;
    std::vector<CandidateTrace> ChoosenPaths;

    //    //static size_t Time_InitSubP;
    //    static size_t Time_InitSubPatches;
    //    static size_t Time_TraceSubPatches;
    //    static size_t Time_RetrieveSubPatches;
    ////    static size_t Time_AddingEmitters;
    ////    static size_t Time_AddingLoops;
    ////    static size_t Time_AddingBorder;


    //                std::cout<<"t0:"<<lupd_time_menage_sel<<std::endl;
    //                std::cout<<"t1:"<<lupd_time_getting_diff<<std::endl;
    //                std::cout<<"t2:"<<lupd_time_changed_part<<std::endl;
    //                std::cout<<"t3:"<<lupd_time_store_part<<std::endl;
    //                std::cout<<"t4:"<<lupd_time_update_part_type<<std::endl;
    //                std::cout<<"t5:"<<lupd_time_update_part<<std::endl;
private:
    std::vector<size_t> VerticesNeeds;

    //std::vector<std::vector<ScalarType> > EdgeL;
    std::map<std::pair<size_t,size_t>,ScalarType> EdgeL;
    //std::unordered_map<std::pair<size_t,size_t>,ScalarType> EdgeL;
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

public:

    void InitEdgeL()
    {
        //std::cout<<"** INIT EDGE LENGHT **"<<std::endl;

        EdgeL.clear();
        InitEdgeLBorders();

        for (size_t i=0;i<ChoosenPaths.size();i++)
            AddEdgeL(ChoosenPaths[i]);

        //std::cout<<"** DONE **"<<std::endl;
    }

private:
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

public:

    void SampleLoopEmitters(bool filter_border,size_t fixed_num=0)
    {
        //then add internal one for loops and other tracing
        std::vector<size_t> StartingNodes;
        size_t sampleNum=MIN_SAMPLES;
        size_t minNumHard=MIN_SAMPLES_HARD;

        assert(sample_ratio<=1);
        assert(sample_ratio>0);

        if (fixed_num==0)
        {
            sampleNum=Mesh().vert.size()*sample_ratio;//floor(sqrt(Mesh().vert.size())+0.5)*10*sample_ratio;
            sampleNum=std::max(sampleNum,(size_t)MIN_SAMPLES);
            sampleNum=std::min(sampleNum,(size_t)MAX_SAMPLES);
        }
        else
        {
            sampleNum=fixed_num;
            minNumHard=sampleNum;
        }

        if (filter_border)
        {
            VertexFieldQuery<MeshType>::SamplePoissonNodesBorderFiltering(VFGraph,sampleNum,StartingNodes,minNumHard,DebugMsg);
        }
        else
        {
            VertexFieldQuery<MeshType>::SamplePoissonNodes(VFGraph,sampleNum,StartingNodes,minNumHard,DebugMsg);
        }
        //then check each connected component has samples
        if (DebugMsg)
        {
            std::cout<<" TARGET "<<sampleNum *2 <<" NODES "<<std::endl;
            std::cout<<" SAMPLED INITIAL "<<StartingNodes.size()<<" NODES "<<std::endl;
        }

        for (size_t i=0;i<StartingNodes.size();i++)
        {
            size_t IndexV=VertexFieldGraph<MeshType>::NodeVertI(StartingNodes[i]);
            if (VertType[IndexV]!=TVInternal)continue;
            NodeEmitterTypes[StartingNodes[i]]=TVInternal;
        }
    }

private:

    //THIS INITIALIZE THE EMITTERS FOR EACH VERTEX
    void InitEmitters()
    {
        NodeEmitterTypes.clear();
        NodeReceiverTypes.clear();
        NodeEmitterTypes.resize(VFGraph.NumNodes(),TVNone);
        NodeReceiverTypes.resize(VFGraph.NumNodes(),TVNone);


        std::vector<std::vector<CoordType> > VertFlatDir,VertOrthoDir;
        VertexEmitter<MeshType>::GetOrthoFlatDirections(VFGraph.Mesh(),VertFlatDir,VertOrthoDir);

        for (size_t i=0;i<Mesh().vert.size();i++)
        {

            if (VertType[i]==TVConvex)//convex ones, no emitters
                continue;

            if (VertType[i]==TVInternal)//this will be with sampling
                continue;

            if (VertType[i]==TVFlat)
            {
                size_t Emitter,Receiver;
                VertexEmitter<MeshType>::ComputeFlatEmitterReceivers(VFGraph,VertOrthoDir,i,Emitter,Receiver);
                //VertexEmitter<MeshType>::ComputeFlatEmitterReceivers(VFGraph,VertOrthoDir,VertFlatDir,i,Emitter,Receiver);
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
                //VertexEmitter<MeshType>::ComputeNarrowEmitterReceivers(VFGraph,VertOrthoDir,VertFlatDir,i,Emitter,Receiver);
                VertexEmitter<MeshType>::ComputeNarrowEmitterReceivers(VFGraph,VertFlatDir,i,Emitter,Receiver);
                assert(Emitter!=Receiver);
                assert(NodeEmitterTypes[Emitter]==TVNone);
                NodeEmitterTypes[Emitter]=TVNarrow;
                assert(NodeReceiverTypes[Receiver]==TVNone);
                NodeReceiverTypes[Receiver]=TVNarrow;
            }
        }

        SampleLoopEmitters(false);

    }

    //remove the connection with singularities that are not concave
    //or narrow and we should not trace from them
    void InvalidateNonConcaveSing()
    {
        std::set<size_t> InvalidatedNodes;
        for (size_t i=0;i<VertType.size();i++)
        {
            if (VertType[i]!=TVFlat)continue;
            if (!VFGraph.IsSingVert[i])continue;
            //get all nodes associated
            std::vector<size_t> Nodes;
            VertexFieldGraph<MeshType>::IndexNodes(i,Nodes);
            InvalidatedNodes.insert(Nodes.begin(),Nodes.end());
        }

        VFGraph.RemoveConnections(InvalidatedNodes);
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
            if (SelfInt){
                SelfIntN++;
                DiscardedCandidates.push_back(Candidates[i]);
                continue;
            }
            if (expanded)
                ExpandedCandidates.push_back(Candidates[i]);
            else
                DiscardedCandidates.push_back(Candidates[i]);
        }
        if (DebugMsg)
            std::cout<<"Self Intersections "<<SelfIntN<<std::endl;

        //all expanded do nothing (already chenged in place)
        if (ExpandedCandidates.size()==Candidates.size())return;
        Candidates=ExpandedCandidates;
    }

    //    size_t NumEmitters(MeshType &mesh,
    //                       std::vector<size_t> &PatchFaces,
    //                       std::vector<size_t> &VerticesNeeds)
    //    {
    //        size_t ret=0;
    //        vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);
    //        for (size_t i=0;i<PatchFaces.size();i++)
    //            for (size_t j=0;j<3;j++)
    //            {
    //                size_t IndexV=vcg::tri::Index(mesh,mesh.face[PatchFaces[i]].V(j));
    //                if (mesh.vert[IndexV].IsV())continue;
    //                mesh.vert[IndexV].SetV();
    //                ret+=VerticesNeeds[IndexV];
    //            }
    //        return ret;
    //    }

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
            if (VertType[i]==TVConcave)
                VerticesNeeds[i]=std::max(NumEmitters(i)-1,(int)Concave_Need);
    }

    //INITIALIZE THE STRUCTURES
    void InitInternalStructures()
    {
        //InitOrthoDirections();
        if (DebugMsg)
        {
            std::cout<<"** STARTING INITIALIZATION **"<<std::endl;
            std::cout<<"* Initializing Vertex Type"<<std::endl;
        }

        InitVertType();

        if (DebugMsg)
            std::cout<<"* Setting Emitters type"<<std::endl;


        InitEmitters();

        if (DebugMsg)
            std::cout<<"* Invalidating Tangent Nodes"<<std::endl;

        InvalidateTangentNodes();

        if (DebugMsg)
            std::cout<<"* Invalidating Singularities"<<std::endl;

        InvalidateNonConcaveSing();

        if (DebugMsg)
            std::cout<<"* Setting Vertices needs"<<std::endl;
        InitVerticesNeeds();

        if (DebugMsg)
            std::cout<<"** INITIALIZATION DONE **"<<std::endl;

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
        size_t selectedV=VFGraph.Select(CanReceive);
        if (DebugMsg)
            std::cout<<"Num Receivers "<<selectedV<<std::endl;

        for (size_t i=0;i<Candidates.size();i++)
        {
            if (Candidates[i].Updated)continue;

            //            if (DebugMsg)
            //                std::cout<<"Updating Candidate "<<std::endl;

            UpdateCandidate(VFGraph,Candidates[i],Drift,MaxNarrowWeight);//,DebugMsg);
        }
        //finally erase the non traced ones
        CleanNonTracedCandidates();
        if (DebugMsg)
            std::cout<<"Done Updating "<<selectedV<<std::endl;
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

    void AddChoosen(CandidateTrace &CurrC)
    {
        //        MeshType outMesh;
        //        CandidateTrace TestC=CurrC;
        //        TestC.PathNodes.resize(4);
        //        TestC.IsLoop=false;
        //        PatchManager<MeshType>::MeshTrace(VFGraph,TestC,VFGraph.Mesh().bbox.Diag()*0.001,outMesh);
        //        vcg::tri::io::ExporterPLY<MeshType>::Save(outMesh,"test_trace.ply");
        //        vcg::tri::io::ExporterPLY<MeshType>::Save(VFGraph.Mesh(),"test_domain.ply");


        //        VertexFieldQuery<MeshType>::CutEarPath(VFGraph,CurrC.PathNodes,CurrC.IsLoop);

        ChoosenPaths.push_back(CurrC);
        assert(CurrC.PathNodes.size()>=2);

        UpdateVertNeeds(CurrC.PathNodes);

        //add on the table
        AddEdgeNodes<MeshType>(CurrC.PathNodes,CurrC.IsLoop,EDirTable);

        AddEdgeL(CurrC);
    }

    void ChooseGreedyByLength(bool UseVertNeeds=true,
                              bool UsePartitionNeeds=false)
    {
        //assert(!UsePartitionNeeds);//need to be implemented
        //InitCandidatesPathLenghts();
        SortCandidatesByLenghts();

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
            if ((UseVertNeeds)&&((VerticesNeeds[IndexV0]==0)&&(VerticesNeeds[IndexV1]==0)))
            {
                DiscardedCandidates.push_back(Candidates[i]);
                continue;
            }


            //bool collide=CollideWithChoosen(CurrTrace,IsCurrLoop,StartConflPath);
            bool collide = CollideWithCandidateSet(VFGraph,Candidates[i],ChoosenPaths);
            if (collide)
            {
                DiscardedCandidates.push_back(Candidates[i]);
                continue;
            }

            if ((UsePartitionNeeds)&&
                    (!SplitSomeNonOKPartition(VFGraph,Candidates[i],Partitions,FacePartitions,PartitionType)))
            {
                DiscardedCandidates.push_back(Candidates[i]);
                continue;
            }

            //ChoosenPaths.push_back(Candidates[i]);
            //assert(ChoosenPaths.back().PathNodes.size()>=2);
            //UpdateVertNeeds(CurrTrace);
            AddChoosen(Candidates[i]);

            if (UsePartitionNeeds)
            {
                //std::vector<size_t> LastAdded(1,ChoosenPaths.size()-1);
                //LazyUpdatePartitions(LastAdded);
                LazyUpdatePartitions();
            }
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
                              bool UsePartitionNeeds,
                              size_t &time_sort,
                              size_t &time_collide,
                              size_t &time_solve_check,
                              size_t &time_update_dist,
                              size_t &time_add)
    {
        std::vector<bool> To_Delete(Candidates.size(),false);
        int t0=clock();
        SortCandidatesByDistances();
        int t1=clock();
        time_sort+=t1-t0;

        size_t lupd_time_menage_sel=0;
        size_t lupd_time_getting_diff=0;
        size_t lupd_time_changed_part=0;
        size_t lupd_time_store_part=0;
        size_t lupd_time_update_part_type=0;
        size_t lupd_time_update_part=0;

        for (size_t i=0;i<Candidates.size();i++)
        {
            bool IsCurrLoop=Candidates[i].IsLoop;

            if (UseVertNeeds)
            {
                assert(!IsCurrLoop);
                if (!SolveVertexNeed(Candidates[i]))
                {
                    To_Delete[i]=true;
                    continue;
                }
            }
            //check if collide
            //bool collide=CollideWithChoosen(CurrTrace,IsCurrLoop);

            t0=clock();
            bool collide = CollideWithCandidateSet(VFGraph,Candidates[i],ChoosenPaths);
            t1=clock();
            time_collide+=t1-t0;

            if (collide)
            {
                To_Delete[i]=true;
                continue;
            }

            if (UsePartitionNeeds)
            {
                if (!SplitSomeNonOKPartition(VFGraph,Candidates[i],Partitions,FacePartitions,PartitionType))
                {
                    To_Delete[i]=true;
                    //std::cout<<"Not Needed by Partition "<<std::endl;
                    continue;
                }
            }

            t0=clock();
            AddChoosen(Candidates[i]);
            t1=clock();
            time_add+=t1-t0;

            //update distances
            t0=clock();
            UpdateDistancesWithLastChoosen();

            t1=clock();
            time_update_dist+=t1-t0;

            if (UsePartitionNeeds)
            {
                t0=clock();
                //set the last one to update
                //std::vector<size_t> LastAdded(1,ChoosenPaths.size()-1);
                //LazyUpdatePartitions(LastAdded);

                LazyUpdatePartitions(lupd_time_menage_sel,
                                     lupd_time_getting_diff,
                                     lupd_time_changed_part,
                                     lupd_time_store_part,
                                     lupd_time_update_part_type,
                                     lupd_time_update_part);

                //                std::cout<<"t0:"<<lupd_time_menage_sel<<std::endl;
                //                std::cout<<"t1:"<<lupd_time_getting_diff<<std::endl;
                //                std::cout<<"t2:"<<lupd_time_changed_part<<std::endl;
                //                std::cout<<"t3:"<<lupd_time_store_part<<std::endl;
                //                std::cout<<"t4:"<<lupd_time_update_part_type<<std::endl;
                //                std::cout<<"t5:"<<lupd_time_update_part<<std::endl;

                //WriteInfo();
                //UpdateDistancesWithLastChoosen();
                //                UpdatePartitionsFromChoosen(true);
                t1=clock();
                time_solve_check+=t1-t0;
            }

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
        if (DebugMsg)
            std::cout<<"Updating Partitions Greedy Distance"<<std::endl;

        //THIS SHOULD SPEED UP
        LazyUpdatePartitions();
        //UpdatePartitionsFromChoosen(true);

        CurrNodeDist.clear();
        if (DebugMsg)
            std::cout<<"Init All Nodes Distance test"<<std::endl;
        InitNodeDistances();
        if (DebugMsg)
            std::cout<<"Done"<<std::endl;
        size_t time_sort=0;
        size_t time_collide=0;
        size_t time_solve_check=0;
        size_t time_update_dist=0;
        size_t time_add=0;
        //std::cout<<"adding candidates by bigger distance"<<std::endl;
        while (ChooseNextByDistance(UseVertNeeds,UsePartitionNeeds,
                                    time_sort,time_collide,time_solve_check,
                                    time_update_dist,time_add))
        {}
        //        std::cout<<"done"<<std::endl;
        //        //            step=(step+1)%100;
        //        //            if (step!=0)continue;
        //        std::cout<<"Time Sort "<<time_sort/(ScalarType)CLOCKS_PER_SEC<<std::endl;
        //        std::cout<<"Time Collide "<<time_collide/(ScalarType)CLOCKS_PER_SEC<<std::endl;
        //        std::cout<<"Time Solve Chek "<<time_solve_check/(ScalarType)CLOCKS_PER_SEC<<std::endl;
        //        std::cout<<"Time Update Dist "<<time_update_dist/(ScalarType)CLOCKS_PER_SEC<<std::endl;
        //        std::cout<<"Time Update Add "<<time_add/(ScalarType)CLOCKS_PER_SEC<<std::endl;

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

    TypeVert GetTypeOfNode(size_t IndexN)
    {
        size_t IndexV=VertexFieldGraph<MeshType>::NodeVertI(IndexN);
        return (VertType[IndexV]);
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

    void GetVertexType(TypeVert TVert,std::vector<size_t> &VertexSet)
    {
        VertexSet.clear();

        for (size_t i=0;i<VertType.size();i++)
            if (VertType[i]==TVert)
                VertexSet.push_back(i);
    }

    void GetEmitterType(const TypeVert EmitType,std::vector<size_t> &NodeEmitType)
    {
        NodeEmitType.clear();
        for (size_t i=0;i<NodeEmitterTypes.size();i++)
            if (NodeEmitterTypes[i]==EmitType)
                NodeEmitType.push_back(i);
    }

    void GetAllEmitter(std::vector<size_t> &NodeEmitType)
    {
        for (size_t i=0;i<NodeEmitterTypes.size();i++)
            if (NodeEmitterTypes[i]!=TVNone)
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

    void AddSingleEdgePath(const size_t &IndexN0,
                           const size_t &IndexN1)
    {
        CandidateTrace CurrC;
        CurrC.FromType=GetTypeOfNode(IndexN0);
        CurrC.ToType=GetTypeOfNode(IndexN1);
        CurrC.IsLoop=false;
        CurrC.PathNodes.push_back(IndexN0);
        CurrC.PathNodes.push_back(IndexN1);
        CurrC.Priority=0;
        CurrC.Unremovable=false;
        CurrC.Updated=true;
        ChoosenPaths.push_back(CurrC);
    }

    int NumEmitterType(const TypeVert EmitType)
    {
        size_t NumE=0;
        for (size_t i=0;i<NodeEmitterTypes.size();i++)
            if (NodeEmitterTypes[i]==EmitType)
                NumE++;
        return NumE;
    }

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


    bool MergeIfPossible(const std::vector<size_t> &NodesSeq0,
                         const std::vector<size_t> &NodesSeq1,
                         std::vector<size_t> &NewNodeSeq,
                         bool &NewIsLoop)
    {
        size_t Node0_Begin=NodesSeq0[0];
        size_t Node0_End=NodesSeq0.back();
        size_t Node1_Begin=NodesSeq1[0];
        size_t Node1_End=NodesSeq1.back();
        std::vector<size_t> NewTrace0,NewTrace1;
        NewIsLoop=false;
        if (Node0_End==Node1_Begin)
        {
            NewTrace0=NodesSeq0;
            //remove the duplicated
            NewTrace0.pop_back();
            NewTrace1=NodesSeq1;
            //then append trace0 with trace1
            NewNodeSeq=NewTrace0;
            NewNodeSeq.insert(NewNodeSeq.end(),NewTrace1.begin(),NewTrace1.end());

            if (NewNodeSeq[0]==NewNodeSeq.back())
            {
                NewNodeSeq.pop_back();
                NewIsLoop=true;
            }

            return true;
        }
        if (Node1_End==Node0_Begin)
        {
            NewTrace0=NodesSeq1;
            //remove the duplicated
            NewTrace0.pop_back();
            NewTrace1=NodesSeq0;

            //then append trace0 with trace1
            NewNodeSeq=NewTrace0;
            NewNodeSeq.insert(NewNodeSeq.end(),NewTrace1.begin(),NewTrace1.end());

            if (NewNodeSeq[0]==NewNodeSeq.back())
            {
                NewNodeSeq.pop_back();
                NewIsLoop=true;
            }

            return true;
        }
        if (Node0_End==VertexFieldGraph<MeshType>::TangentNode(Node1_End))
        {
            NewTrace0=NodesSeq0;
            //remove the duplicated
            NewTrace0.pop_back();
            NewTrace1=NodesSeq1;
            //invert them
            VertexFieldGraph<MeshType>::TangentNodes(NewTrace1);
            std::reverse(NewTrace1.begin(),NewTrace1.end());
            //then append trace0 with trace1
            NewNodeSeq=NewTrace0;
            NewNodeSeq.insert(NewNodeSeq.end(),NewTrace1.begin(),NewTrace1.end());

            if (NewNodeSeq[0]==NewNodeSeq.back())
            {
                NewNodeSeq.pop_back();
                NewIsLoop=true;
            }

            return true;
        }
        if (VertexFieldGraph<MeshType>::TangentNode(Node0_Begin)==Node1_Begin)
        {
            NewTrace0=NodesSeq0;
            //invert them
            VertexFieldGraph<MeshType>::TangentNodes(NewTrace0);
            std::reverse(NewTrace0.begin(),NewTrace0.end());
            //remove the duplicated
            NewTrace0.pop_back();
            NewTrace1=NodesSeq1;

            //then append trace0 with trace1
            NewNodeSeq=NewTrace0;
            NewNodeSeq.insert(NewNodeSeq.end(),NewTrace1.begin(),NewTrace1.end());

            if (NewNodeSeq[0]==NewNodeSeq.back())
            {
                NewNodeSeq.pop_back();
                NewIsLoop=true;
            }

            return true;
        }
        return false;
    }

    size_t TestedMerged;
    size_t DoneMerged;

    bool MergeContiguousPathStep(std::vector<CandidateTrace> &TraceSet)
    {
        if (TraceSet.size()<2)return false;
        bool HasMerged=false;
        for (size_t i=0;i<TraceSet.size()-1;i++)
            for (size_t j=(i+1);j<TraceSet.size();j++)
            {
                if (i==j)continue;
                if (TraceSet[i].IsLoop)continue;
                if (TraceSet[j].IsLoop)continue;
                if(TraceSet[i].PathNodes.size()==0)continue;
                if(TraceSet[j].PathNodes.size()==0)continue;

                std::vector<size_t> NewNodeSeq;
                bool NewIsLoop=false;
                bool has_merged=MergeIfPossible(TraceSet[i].PathNodes,TraceSet[j].PathNodes,NewNodeSeq,NewIsLoop);
                TestedMerged++;
                if (has_merged)
                {
                    TraceSet[j].PathNodes.clear();
                    TraceSet[i].PathNodes=NewNodeSeq;
                    TraceSet[i].IsLoop=NewIsLoop;
                    TraceSet[i].FromType=GetTypeOfNode(TraceSet[i].PathNodes[0]);
                    TraceSet[i].ToType=GetTypeOfNode(TraceSet[i].PathNodes.back());
                    DoneMerged++;
                    HasMerged=true;
                    //return true;
                }
            }
        //return false;
        return HasMerged;
    }



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
            //            GetEmitterType(TVFlat,EmitterNodes);
            //            for (size_t i=0;i<EmitterNodes.size();i++)
            //            {
            //                if (MustDisable[EmitterNodes[i]]==false)
            //                    CanEmit[EmitterNodes[i]]=true;
            //            }
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
            return(TraceFrom(FromType,ToType,TrType,CanEmit,CanReceive,Shortest,UsePartitionNeeds));
        else
            return(TraceFrom(FromType,ToType,TrType,CanEmit,CanReceive,Fartest,UsePartitionNeeds));
    }

public:

    void InitEdgeDirTable()
    {
        //        //std::cout<<"*** INITIALIZING TABLE ***"<<std::endl;
        //        std::vector<size_t> ConvexIdx,ConcaveIdx;

        //        for (size_t i=0;i<VertType.size();i++)
        //        {
        //            if (VertType[i]==TVConvex)
        //                ConvexIdx.push_back(i);
        //            if (VertType[i]==TVConcave)
        //                ConcaveIdx.push_back(i);
        //        }
        EDirTable.Init(VertType);
        //add the borders
        AddBorder<MeshType>(VFGraph,EDirTable);

        //add the other traced nodes
        std::vector<std::vector<size_t> > Nodes;
        std::vector<bool> IsLoop;
        for (size_t i=0;i<ChoosenPaths.size();i++)
        {
            Nodes.push_back(ChoosenPaths[i].PathNodes);
            IsLoop.push_back(ChoosenPaths[i].IsLoop);
        }
        AddEdgeNodes<MeshType>(Nodes,IsLoop,EDirTable);

    }

    EdgeDirectionTable EDirTable;

private:




    void UpdatePartitionType(size_t Index)
    {
        assert(Index<PartitionCorners.size());
        assert(Index<Partitions.size());

        //REMOVE THE REST
        if (PatchInfos[Index].Genus!=1)
        {
            PartitionType[Index]=NonDisk;
            return;
        }

        if (PatchInfos[Index].NumEmitters>0)//&&(!AllowRemoveConcave))
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

        if (check_quality_functor)
        {
            PatchQualityFunctor PFunct;
            MeshType m;
            bool UseInternalCuts=(AllowDarts||AllowSelfGluedPatch);
            GetPatchMesh(Index,m,UseInternalCuts);

            //PFunct.ContinuousCheckSelfInt()=ContinuosCheckUVInt;
            PatchInfos[Index].QualityFunctorValue=PFunct(m);
            //PFunct.ContinuousCheckSelfInt()=true;
            if (PatchInfos[Index].QualityFunctorValue>0)
            {
                PartitionType[Index]=NoQualityMatch;
                return;
            }
        }
        //        if ((check_quality_functor)&&(PatchInfos[Index].QualityFunctorValue>0))
        //        {
        //            PartitionType[Index]=NoQualityMatch;
        //            return;
        //        }
        //        //optional QualityFunctor
        //        if (check_quality_functor)
        //        {
        //            PatchQualityFunctor PFunct;
        //            MeshType m;
        //            GetPatchMesh(Index,m);
        //            bool IsOKQ=PFunct(m);
        //            if (!IsOKQ)
        //            {
        //                PartitionType[Index]=NoQualityMatch;
        //                return;
        //            }
        //        }
        //            std::cout<<"isOK"<<std::endl;
        //        else
        //            std::cout<<"NotOK"<<std::endl;

        PartitionType[Index]=IsOK;
    }

public:

    void InitPartitionsType()
    {
        PartitionType.clear();
        PartitionType.resize(Partitions.size(),IsOK);
        //int t0=clock();
        //InitEdgeL();
        //int t1=clock();
        //if the singularity have been used then ccability of valence 4 is not computed

        //optional QualityFunctor


        PatchManager<MeshType>::GetPatchInfo(Mesh(),Partitions,PartitionCorners,VerticesNeeds,EdgeL,
                                             PatchInfos,avgEdge*CClarkability,match_valence,
                                             AllowDarts,AllowSelfGluedPatch,
                                             CheckQuadrangulationLimits);
        //update quality functor
        //        if (check_quality_functor)
        //        {
        //            PatchQualityFunctor PFunct;
        //            for (size_t i=0;i<Partitions.size();i++)
        //            {
        //                MeshType m;
        //                bool UseInternalCuts=(AllowDarts||AllowSelfGluedPatch);
        //                GetPatchMesh(i,m,UseInternalCuts);
        //                PatchInfos[i].QualityFunctorValue=PFunct(m);
        //            }
        //        }
        //int t2=clock();

        for (size_t i=0;i<Partitions.size();i++)
            UpdatePartitionType(i);

        //int t3=clock();
        //                std::cout<<"Time Init Lenghts "<<t1-t0<<std::endl;
        //                std::cout<<"Time Get PInfo "<<t2-t1<<std::endl;
        //                std::cout<<"Time Update Type "<<t3-t2<<std::endl;
    }


    std::vector<std::vector<bool> > EdgeSel0;
    std::vector<std::vector<bool> > EdgeSel1;

public:

    void AddSurroundingPos(VertexType *v,std::vector<vcg::face::Pos<FaceType> > &Pos)
    {
        std::vector<FaceType*> faceVec;
        std::vector<int> index;
        vcg::face::VFStarVF(v,faceVec,index);
        for (size_t i=0;i<faceVec.size();i++)
            Pos.push_back(vcg::face::Pos<FaceType>(faceVec[i],index[i]));
    }

    void LazyUpdatePartitions(size_t &time_menage_sel,
                              size_t &time_getting_diff,
                              size_t &time_changed_part,
                              size_t &time_store_part,
                              size_t &time_update_part_type,
                              size_t &time_update_part)
    {
        int t0=clock();
        //get the old border nodes
        PatchManager<MeshType>::SaveEdgeSel(Mesh(),EdgeSel0);
        //select with the new ones
        //SelectMeshPatchBorders(Mesh(),ChoosenPaths);
        SelectMeshPatchBorders(VFGraph,ChoosenPaths);
        //get the new ones
        PatchManager<MeshType>::SaveEdgeSel(Mesh(),EdgeSel1);
        int t1=clock();

        //then get the difference in terms of faces
        std::vector<vcg::face::Pos<FaceType> > changedF;
        assert(EdgeSel0.size()==EdgeSel1.size());
        for (size_t i=0;i<EdgeSel0.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                if (EdgeSel0[i][j]==EdgeSel1[i][j])continue;
                FaceType *f=&Mesh().face[i];
                vcg::face::Pos<FaceType> currPos(f,j);
                changedF.push_back(currPos);

                //if the vertex is narrow or concave then add the ones around
                size_t IndexV0=vcg::tri::Index(Mesh(),currPos.V());
                size_t IndexV1=vcg::tri::Index(Mesh(),currPos.VFlip());
                if ((VertType[IndexV0]==TVNarrow)||
                        (VertType[IndexV0]==TVConcave))
                {
                    AddSurroundingPos(currPos.V(),changedF);
                    //                    VertexType *v=currPos.V();
                    //                    vcg::face::VFStarVF(v,faceVec,index);
                }
                if ((VertType[IndexV1]==TVNarrow)||
                        (VertType[IndexV1]==TVConcave))
                {
                    AddSurroundingPos(currPos.VFlip(),changedF);
                    //                    VertexType *v=currPos.V();
                    //                    vcg::face::VFStarVF(v,faceVec,index);
                }
            }
        }
        int t2=clock();
        //std::cout<<"Changed "<<changedF.size()<<" Faces"<<std::endl;
        //get the partitions need to be updated
        std::set<int> ChangedPartitions;
        for (size_t i=0;i<changedF.size();i++)
        {
            size_t IndexF=vcg::tri::Index(Mesh(),changedF[i].F());
            int indexPatch=FacePartitions[IndexF];
            assert(indexPatch<(int)Partitions.size());
            assert(indexPatch>=0);
            ChangedPartitions.insert(indexPatch);
        }
        int t3=clock();

        //then save the old one not changed
        std::vector<std::vector<size_t> > PartitionsSwap;
        std::vector<PatchType> PartitionTypeSwap;
        std::vector<std::vector<size_t> > PartitionCornersSwap;
        std::vector<PatchInfo<ScalarType> > PatchInfosSwap;
        assert(Partitions.size()==PartitionType.size());
        assert(Partitions.size()==PartitionCorners.size());
        assert(Partitions.size()==PatchInfos.size());

        //        std::cout<<"Changed Partitions"<<std::endl;
        //        for (std::set<int>::iterator It=ChangedPartitions.begin();
        //             It!=ChangedPartitions.end();It++)
        //        {
        //            std::cout<<(*It)<<",";
        //        }
        //        std::cout<<std::endl;

        for (size_t i=0;i<Partitions.size();i++)
        {
            if (ChangedPartitions.count(i)>0)continue;//need to be updated
            PartitionsSwap.push_back(Partitions[i]);
            PartitionTypeSwap.push_back(PartitionType[i]);
            PartitionCornersSwap.push_back(PartitionCorners[i]);
            PatchInfosSwap.push_back(PatchInfos[i]);
        }
        int t4=clock();

        //update around the changed ones
        UpdatePatchAround(changedF);

        int t5=clock();
        //        std::cout<<"First "<<Partitions.size()<<" Has Changed"<<std::endl;
        //        std::cout<<"Then next "<<PartitionsSwap.size()<<" Has Put Back"<<std::endl;
        //append partitions to old ones
        Partitions.insert(Partitions.end(),PartitionsSwap.begin(),PartitionsSwap.end());
        PartitionType.insert(PartitionType.end(),PartitionTypeSwap.begin(),PartitionTypeSwap.end());
        PartitionCorners.insert(PartitionCorners.end(),PartitionCornersSwap.begin(),PartitionCornersSwap.end());
        PatchInfos.insert(PatchInfos.end(),PatchInfosSwap.begin(),PatchInfosSwap.end());

        PatchManager<MeshType>::DerivePerFacePartition(Mesh(),Partitions,FacePartitions);
        int t6=clock();
        time_menage_sel+=t1-t0;
        time_getting_diff+=t2-t1;
        time_changed_part+=t3-t2;
        time_store_part+=t4-t3;
        time_update_part_type+=t5-t4;
        time_update_part+=t6-t5;

        Time_UpdatePartitionsLazy+=t6-t0;
    }


    void LazyUpdatePartitions()
    {
        size_t lupd_time_menage_sel=0;
        size_t lupd_time_getting_diff=0;
        size_t lupd_time_changed_part=0;
        size_t lupd_time_store_part=0;
        size_t lupd_time_update_part_type=0;
        size_t lupd_time_update_part=0;
        LazyUpdatePartitions(lupd_time_menage_sel,
                             lupd_time_getting_diff,
                             lupd_time_changed_part,
                             lupd_time_store_part,
                             lupd_time_update_part_type,
                             lupd_time_update_part);
    }

    void UpdatePartitionsFromChoosen(bool UpdateType=true)
    {

        //std::cout<<"Updating Partitions"<<std::endl;
        //        if ((UpdateType)&&
        //            (EdgeSel0.size()==Mesh().face.size())&&
        //            (this->FirstBorder))
        //        {
        //            LazyUpdatePartitions();
        //            return;
        //        }

        int t0=clock();
        //bool IsOk=SelectMeshPatchBorders(Mesh(),ChoosenPaths);//SelectBorders();
        bool IsOk= SelectMeshPatchBorders(VFGraph,ChoosenPaths);
        if (!IsOk)
        {
            MeshType traceMesh;
            std::vector<bool> Selected(ChoosenPaths.size(),false);
            //Selected[i]=true;
            PatchManager<MeshType>::MeshTraces(VFGraph,ChoosenPaths,Selected,traceMesh);
            //vcg::tri::io::ExporterPLY<MeshType>::Save(VFGraph.Mesh(),"double_direction_domain.ply");
            vcg::tri::io::ExporterPLY<MeshType>::Save(traceMesh,"double_direction_error.ply",vcg::tri::io::Mask::IOM_FACECOLOR);
        }
        assert(IsOk);
        std::vector<size_t> StartF;
        for (size_t i=0;i<Mesh().face.size();i++)
            StartF.push_back(i);

        RetrievePatchesFromSelEdges(Mesh(),StartF,Partitions);

        PatchManager<MeshType>::DerivePerFacePartition(Mesh(),Partitions,FacePartitions);

        //        int t1=clock();

        //        int test_t0=clock();
        //        FindPartitionsCorners(VFGraph,VertType,ChoosenPaths,Partitions,PartitionCorners);
        //        int test_t1=clock();
        FindCorners<MeshType>(EDirTable,VFGraph.Mesh(),Partitions,PartitionCorners);

        if ((AllowDarts)||(AllowSelfGluedPatch))
            UpdateCornersWithInternalCuts();

        if (UpdateType)
        {
            InitPartitionsType();
        }
        //WriteInfo();
        int t3=clock();
        Time_UpdatePartitionsTotal+=t3-t0;

        //std::cout<<"Ended Updating Partitions"<<std::endl;

        //        std::cout<<"** Timing Update Partitions **"<<std::endl;
        //        std::cout<<"Time Derive Patch "<<(t1-t0)/(ScalarType)CLOCKS_PER_SEC<<std::endl;
        //        std::cout<<"Time Find Corners "<<(t2-t1)/(ScalarType)CLOCKS_PER_SEC<<std::endl;
        //        std::cout<<"Time Update Partitions "<<(t3-t2)/(ScalarType)CLOCKS_PER_SEC<<std::endl<<std::endl;
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

        Traceable=std::vector<bool>(VFGraph.NumNodes(),true);
        //        InitEdgeDirTable();

        EdgeSel0.clear();
        EdgeSel0.resize(Mesh().face.size(),std::vector<bool>(3,false));
        EdgeSel1.clear();
        EdgeSel1.resize(Mesh().face.size(),std::vector<bool>(3,false));
        //InitTraceableBorders();
    }


    bool HasIncompleteEmitter()
    {
        for (size_t i=0;i<PartitionType.size();i++)
            if (PartitionType[i]==HasEmitter)return true;
        return false;
    }

private:


    bool HasTerminated()
    {
        for (size_t i=0;i<PartitionType.size();i++)
            if (PartitionType[i]!=IsOK)return false;
        return true;
    }

    bool PathHasConcaveNarrowVert(size_t IndexPath)
    {
        //if (AllowRemoveConcave)return false;
        for (size_t j=0;j<ChoosenPaths[IndexPath].PathNodes.size();j++)
        {
            size_t IndexN=ChoosenPaths[IndexPath].PathNodes[j];
            size_t IndexV=VertexFieldGraph<MeshType>::NodeVertI(IndexN);
            if (!AllowRemoveConcave)
            {
                if ((VertType[IndexV]==TVNarrow)||(VertType[IndexV]==TVConcave))
                    return true;
            }
            else
            {
                if ((VertType[IndexV]==TVNarrow)&&(!AllowDarts))
                    return true;
            }
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
            //next check not needed
            if (AllowDarts)continue;
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


    void UpdatePatchAround(const std::vector<vcg::face::Pos<FaceType> > &FacesPath)
    {
        //int t0=clock();
        //get the indexes of faces
        std::vector<size_t> IdxFaces;
        for (size_t i=0;i<FacesPath.size();i++)
        {
            IdxFaces.push_back(vcg::tri::Index(Mesh(),FacesPath[i].F()));
            IdxFaces.push_back(vcg::tri::Index(Mesh(),FacesPath[i].FFlip()));
        }
        //then retrieve partitions
        RetrievePatchesFromSelEdges(Mesh(),IdxFaces,Partitions);
        //int t1=clock();
        //std::cout<<"Updating "<<Partitions.size()<<" part"<<std::endl;

        PatchManager<MeshType>::DerivePerFacePartition(Mesh(),Partitions,FacePartitions);
        //find corners
        //int t2=clock();
        //        FindPartitionsCorners<MeshType>(VFGraph,VertType,ChoosenPaths,Partitions,PartitionCorners);

        FindCorners<MeshType>(EDirTable,VFGraph.Mesh(),Partitions,PartitionCorners);
        if ((AllowDarts)||(AllowSelfGluedPatch))
            UpdateCornersWithInternalCuts();
        //        int test_t2=clock();
        //        std::cout<<"Time 0 "<<test_t1-test_t0<<std::endl;
        //        std::cout<<"Time 1 "<<test_t2-test_t1<<std::endl;
        //        assert(PartitionCorners==PartitionCorners2);

        //find type
        //int t3=clock();
        InitPartitionsType();
        //int t4=clock();

        //        std::cout<<"t0:"<<t1-t0<<std::endl;
        //        std::cout<<"t1:"<<t2-t1<<std::endl;
        //        std::cout<<"t2:"<<t3-t2<<std::endl;
        //        std::cout<<"t3:"<<t4-t3<<std::endl;
    }

    size_t RMZeroSize;
    size_t RMUnremoveable;
    size_t RMHasConcaveNarrow;
    size_t RMHasTJunction;
    size_t RMHasDeadEnd;
    size_t RMOutOfCornerNum;
    size_t RMNotProfitable;

    bool RemoveIfPossible(size_t IndexPath)//,bool Do_Really_Remove=true)
    {
        //std::cout<<"0"<<std::endl;
        //int t0=clock();
        if (ChoosenPaths[IndexPath].PathNodes.size()==0)
        {
            RMZeroSize++;
            return false;
        }

        if (ChoosenPaths[IndexPath].Unremovable)
        {
            //std::cout<<"Unremoveable"<<std::endl;
            RMUnremoveable++;
            return false;
        }
        //if it includes a concave or narrow then cannot remove
        //        if (!AllowDarts)
        //        {
        if (PathHasConcaveNarrowVert(IndexPath))
        {
            RMHasConcaveNarrow++;
            return false;
        }

        //check if have t junction in the middle
        if ((CheckTJunction)&&(PathHasTJunction(IndexPath)))
        //if (PathHasTJunction(IndexPath))
        {
            RMHasTJunction++;
            return false;
        }
        //}

        //std::cout<<"1"<<std::endl;

        //CHECK ENDPOINTS!
        assert(IndexPath<ChoosenPaths.size());

        //get the old configuration
        std::vector<vcg::face::Pos<FaceType> > FacesPath;
        VFGraph.GetNodesPos(ChoosenPaths[IndexPath].PathNodes,ChoosenPaths[IndexPath].IsLoop,FacesPath);

        //int t1=clock();
        //first time no need to update the quality
        bool swap_quality_funct=check_quality_functor;
        check_quality_functor=false;
        UpdatePatchAround(FacesPath);
        check_quality_functor=swap_quality_funct;

        //int t2=clock();
        //UpdatePartitionsFromChoosen(true);
        std::vector<PatchInfo<ScalarType> > PatchInfos0=PatchInfos;
        //std::vector<std::vector<size_t> > FacePatches0=Partitions;

        //std::cout<<"2"<<std::endl;

        //test removal
        CandidateTrace OldTr=ChoosenPaths[IndexPath];
        ChoosenPaths[IndexPath].PathNodes.clear();
        //remove from the table
        RemoveEdgeNodes<MeshType>(OldTr.PathNodes,OldTr.IsLoop,EDirTable);

        //std::cout<<"3"<<std::endl;
        //check Tjunctions
        if ((!AllowDarts) && HasPathDeadEnd())
        {
            //restore
            RMHasDeadEnd++;
            ChoosenPaths[IndexPath]=OldTr;
            Mesh().SelectPos(FacesPath,true);
            AddEdgeNodes<MeshType>(OldTr.PathNodes,OldTr.IsLoop,EDirTable);
            return false;
        }

        //deselect
        Mesh().SelectPos(FacesPath,false);

        //int t3=clock();
        //UpdatePartitionsFromChoosen(true);
        UpdatePatchAround(FacesPath);
        //int t4=clock();
        std::vector<PatchInfo<ScalarType> > PatchInfos1=PatchInfos;
        //std::vector<std::vector<size_t> > FacePatches1=Partitions;

        if (CheckQuadrangulationLimits)
        {
            for (size_t i=0;i<PatchInfos1.size();i++)
            {
                //((!AllowDarts)&&
                if ((PatchInfos1[i].NumCorners<(int)MIN_ADMITTIBLE)
                        ||(PatchInfos1[i].NumCorners>(int)MAX_ADMITTIBLE))
                {
                    RMOutOfCornerNum++;
                    //restore
                    ChoosenPaths[IndexPath]=OldTr;
                    Mesh().SelectPos(FacesPath,true);
                    AddEdgeNodes<MeshType>(OldTr.PathNodes,OldTr.IsLoop,EDirTable);
                    return false;
                }
            }
        }

        bool CanRemove=true;
        //int t5=clock();
        //std::vector<typename MeshType::ScalarType> QThresold;
        CanRemove=PatchManager<MeshType>::BetterConfiguration(PatchInfos0,PatchInfos1,MinVal,
                                                              MaxVal,CClarkability,avgEdge,
                                                              match_valence,!AllowRemoveConcave,
                                                              DebugMsg);

        //int t6=clock();
        //CanRemove&=Do_Really_Remove;
        if (!CanRemove)
        {
            //restore
            RMNotProfitable++;
            ChoosenPaths[IndexPath]=OldTr;
            Mesh().SelectPos(FacesPath,true);
            AddEdgeNodes<MeshType>(OldTr.PathNodes,OldTr.IsLoop,EDirTable);
            return false;
        }
        //        if (!Do_Really_Remove)
        //        {
        //            //restore
        //            RMNotProfitable++;
        //            ChoosenPaths[IndexPath]=OldTr;
        //            Mesh().SelectPos(FacesPath,true);
        //            AddEdgeNodes<MeshType>(OldTr.PathNodes,OldTr.IsLoop,EDirTable);
        //        }
        //        //END DEBUG CODE, NOT USEFUL

        //        std::cout<<"T0: "<<t1-t0<<std::endl;
        //        std::cout<<"T1: "<<t2-t1<<std::endl;
        //        std::cout<<"T2: "<<t3-t2<<std::endl;
        //        std::cout<<"T3: "<<t4-t3<<std::endl;
        //        std::cout<<"T4: "<<t5-t4<<std::endl;
        //        std::cout<<"T5: "<<t6-t5<<std::endl;
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

    bool RemoveIteration(bool OnlyOne=false)
    {
        RMZeroSize=0;
        RMUnremoveable=0;
        RMHasConcaveNarrow=0;
        RMHasTJunction=0;
        RMHasDeadEnd=0;
        RMOutOfCornerNum=0;
        RMNotProfitable=0;

        bool HasRemoved=false;
        if (DebugMsg)
            std::cout<<"REMOVING ITERATION"<<std::endl;

        //InitEdgeL();
        for (int i=ChoosenPaths.size()-1;i>=0;i--)
        {
            //std::cout<<"test removing "<<i<<std::endl;
            HasRemoved|=RemoveIfPossible(i);
            if (OnlyOne && HasRemoved)break;
        }

        if (HasRemoved)
        {
            RemoveEmptyPaths();
            if (!AllowDarts)
                MergeContiguousPaths(ChoosenPaths);
        }

        if (DebugMsg)
            std::cout<<"DONE!"<<std::endl;

        return HasRemoved;
    }

    //    bool RemoveGraduallyDart()
    //    {
    //        CandidateTrace OldTr=ChoosenPaths[IndexPath];
    //        ChoosenPaths[IndexPath].PathNodes.clear();

    ////        //remove from the table
    ////        RemoveEdgeNodes<MeshType>(OldTr.PathNodes,OldTr.IsLoop,EDirTable);


    //        std::vector<vcg::face::Pos<FaceType> > FacesPath;
    //        VFGraph.GetNodesPos(ChoosenPaths[IndexPath].PathNodes,ChoosenPaths[IndexPath].IsLoop,FacesPath);

    //        //get the mesh
    //        PatchQualityFunctor PFunct;
    //        MeshType m;
    //        GetPatchMesh(Index,m,true);
    //        do
    //        {

    //            //compute the parametrization
    //            ScalarType QVal=PFunct(m);
    //            //evaluate the distortion

    //        }
    ////        PatchInfos[Index].QualityFunctorValue=PFunct(m);
    ////        if (PatchInfos[Index].QualityFunctorValue>0)
    ////        {
    ////            PartitionType[Index]=NoQualityMatch;
    ////            return;
    ////        }
    //    }


    //    void GetNextEdgeFromDistortion(std::vector<vcg::face::Pos<FaceType> > &PathSeq,
    //                                   std::vector<vcg::face::Pos<FaceType> > &Removed,
    //                                   size_t IntSize=5)
    //    {
    //        //initial test, need to be substituted with distortion
    //        if (PathSeq.size()<=IntSize)
    //        {
    //            Removed=IntSize;
    //            PathSeq.clear();
    //            return;
    //        }

    //        for (size_t i=0;i<IntSize;i++)
    //        {
    //            Removed.push_back(PathSeq.back());
    //            PathSeq.pop_back();
    //        }
    //        assert(PathSeq.size()>0);
    //        assert(Removed.size()==IntSize);
    //    }

    //    void SetFaceEdgeS(vcg::face::Pos<FaceType> &Pos,bool value)
    //    {
    //        assert(!Pos.IsBorder());
    //        assert(!Pos.IsEdgeS());

    //        FaceType *f=Pos.F();
    //        int E=Pos.E();
    //        assert(E>=0);
    //        assert(E<3);
    //        if (value)
    //            f->SetFaceEdgeS(E);

    //        //go on the other side
    //        vcg::face::Pos<FaceType> Opp=Pos;
    //        Opp.FlipF();
    //        assert(!Opp.IsEdgeS());
    //        f=Opp.F();
    //        E=Opp.E();
    //        assert(E>=0);
    //        assert(E<3);

    //        if (value)
    //            f->SetFaceEdgeS(E);
    //    }

    //    void SelectFaceEdgeS(std::vector<vcg::face::Pos<FaceType> > &PosSeq,bool value)
    //    {
    //        for (size_t i=0;i<PosSeq.size();i++)
    //          SetFaceEdgeS(PosSeq[i],value);
    //    }



    //    void RemoveGradually(MeshType &m,std::vector<vcg::face::Pos<FaceType> > &PathSeq)
    //    {
    //       //initial test, the sequence should not be selected
    //       for (size_t i=0;i<PathSeq.size();i++)
    //       {
    //           assert(!PathSeq[i].IsBorder());
    //           assert(!PathSeq[i].IsEdgeS());
    //           //assert(!PathSeq[i].IsEdgeS());
    //       }

    //       std::vector<vcg::face::Pos<FaceType> > CurrSeq=PathSeq;
    //       bool has_completed=false;
    //       do{
    //           //cut the mesh along selection
    //           vcg::tri::CutMeshAlongSelectedFaceEdges(m);
    //           //parametrize
    //           ScalarType QVal=PFunct(m);
    //           if (QVal==0)
    //               has_completed=true;
    //           else
    //           {
    //               std::vector<vcg::face::Pos<FaceType> > To_Select;
    //               GetNextEdgeFromDistortion(CurrSeq,To_Select);
    //               SelectFaceEdgeS(To_Select);
    //           }
    //           //in this case had completely re-added all the cut
    //           if (CurrSeq.size()==0)
    //               has_completed=true;
    //       }while (!has_completed);

    //    }

    //    void RemoveGradually(size_t IndexPath)
    //    {
    //        //get the pos of the current Path
    //        std::vector<vcg::face::Pos<FaceType> > FacesPath;
    //        VFGraph.GetNodesPos(ChoosenPaths[IndexPath].PathNodes,ChoosenPaths[IndexPath].IsLoop,FacesPath);
    //        SelectFaceEdgeS(FacesPath,false);
    //        //to be continued

    //        //then get the submesh

    //    }


    //    void RemoveDartIteration()
    //    {
    //        bool HasRemoved=false;
    //        if (DebugMsg)
    //            std::cout<<"REMOVING DART ITERATION"<<std::endl;

    //        //InitEdgeL();

    ////        //set to no test distortion
    ////        this->check_quality_functor=false;
    ////        for (int i=ChoosenPaths.size()-1;i>=0;i--)
    ////        {
    ////            //if (ChoosenPaths[i].IsLoop)continue;
    ////             //bool CurrentRemoved=RemoveIfPossible(i,false);
    ////            //HasRemoved|=CurrentRemoved;

    ////            //if cannot be removed is because of topologycal condition so
    ////            //do not worth to continue and split
    //////            if (!CurrentRemoved)
    //////                ChoosenPaths[i].Unremovable=true;
    //////            else
    //////            {
    ////                std::cout<<"Possible to remove"<<std::endl;
    ////                ChoosenPaths[i].Unremovable=false;
    ////            //}

    ////            //RemoveGradually(i);
    ////        }

    //        SetAllRemovable();

    //        SplitIntoIntervals();

    //        //set distortion again
    //        //this->check_quality_functor=true;

    //        //then try the removal
    //        BatchRemovalOnMesh(false);

    //        if (HasRemoved)
    //        {
    //            RemoveEmptyPaths();
    //        }

    //        if (DebugMsg)
    //            std::cout<<"DONE!"<<std::endl;

    //        //return HasRemoved;
    //    }

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
        //std::cout<<"Size "<<ToSplit.PathNodes.size()<<std::endl;
        assert(StartI<(int)ToSplit.PathNodes.size());
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

    //    void SplitIntoIntervals(CandidateTrace &CTrace,
    //                            std::vector<CandidateTrace> &SplitPaths,
    //                            size_t &sizeEdges=5)
    //    {

    //    }

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
        //THIS SHOULD SPEED UP
        UpdatePartitionsFromChoosen(false);
        //UpdatePartitionsFromChoosen();

        vcg::tri::UpdateSelection<MeshType>::Clear(Mesh());
        ColorByPartitions();
    }

    void SplitIntoIntervals(std::vector<CandidateTrace> &To_Split)
    {
        //first split ito sub paths along the intersection
        //SplitIntoSubPaths();

        //then select each path every sizeEdges step
        vcg::tri::UpdateSelection<MeshType>::VertexClear(Mesh());
        for (size_t i=0;i<To_Split.size();i++)
        {
            size_t sizeEdges=floor(0.5+(((ScalarType)To_Split[i].PathNodes.size())/(ScalarType)subInt));
            //if (sizeEdges<=1)continue;
            sizeEdges=std::max(sizeEdges,(size_t)1);
            if (To_Split[i].Unremovable)continue;
            if (To_Split[i].PathNodes.size()<=sizeEdges)continue;
            for (size_t j=0;j<To_Split[i].PathNodes.size();j+=sizeEdges)
            {
                size_t IndexV=VertexFieldGraph<MeshType>::NodeVertI(To_Split[i].PathNodes[j]);
                Mesh().vert[IndexV].SetS();
            }
            //select the last one
            size_t IndexV=VertexFieldGraph<MeshType>::NodeVertI(To_Split[i].PathNodes.back());
            Mesh().vert[IndexV].SetS();
        }

        std::vector<CandidateTrace> NewTraces;
        for (size_t i=0;i<To_Split.size();i++)
        {
            std::vector<CandidateTrace> Portions;
            SplitIntoSubPathsBySel(To_Split[i],Portions);
            NewTraces.insert(NewTraces.end(),Portions.begin(),Portions.end());
        }

        To_Split.clear();
        To_Split=NewTraces;

        for (size_t i=0;i<To_Split.size();i++)
            assert(To_Split[i].PathNodes.size()>=2);


        //UpdatePartitionsFromChoosen(false);
        vcg::tri::UpdateSelection<MeshType>::VertexClear(Mesh());
        //ColorByPartitions();
    }

    void MergeContiguousPaths()
    {
        MergeContiguousPaths(ChoosenPaths);
    }

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
        size_t NotMatchQ;
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
        PInfo.NotMatchQ=0;
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
            case NoQualityMatch :
                PInfo.NotMatchQ++;
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
        std::cout<<"* Non Match Quality Patches "<<PInfo.NotMatchQ<<std::endl;
        std::cout<<"* With Emitters Patches "<<PInfo.HasEmit<<std::endl;
        for (size_t i=0;i<8;i++)
            std::cout<<"* Patch with  "<<i<<" corners are: "<<PInfo.SizePatches[i]<<std::endl;
        //std::cout<<"* With More Singularities "<<HasMoreSing<<std::endl;
    }

    //    ScalarType PatchDistortion(size_t IndexP)
    //    {

    //    }
    void SmoothPatches(size_t Steps=3,typename MeshType::ScalarType Damp=0.5,bool CheckSurfaceFolds=true)
    {
        if (CheckSurfaceFolds)
            PatchManager<MeshType>::SmoothMeshPatchesFromEdgeSel(Mesh(),Steps,Damp,0.2);
        //PatchManager<MeshType>::SmoothMeshPatches(Mesh(),FacePartitions,Steps,Damp);
        else
            PatchManager<MeshType>::SmoothMeshPatchesFromEdgeSel(Mesh(),Steps,Damp,-1);

        //PatchManager<MeshType>::SmoothMeshPatches(Mesh(),FacePartitions,Steps,Damp,-1);
    }

    void SetAllUnRemoveable()
    {
        for (size_t i=0;i<ChoosenPaths.size();i++)
            ChoosenPaths[i].Unremovable=true;
    }

    void SetAllRemovable()
    {
        for (size_t i=0;i<ChoosenPaths.size();i++)
            ChoosenPaths[i].Unremovable=false;
    }

    //    ScalarType GetPriority(size_t IndexPath)
    //    {
    //        assert(PrioVect.size()==Mesh().vert.size());
    //        ScalarType Avg=0;
    //        for (size_t i=0;i<ChoosenPaths[IndexPath].PathNodes.size();i++)
    //        {
    //            size_t IndexV=VertexFieldGraph<MeshType>::NodeVertI(ChoosenPaths[IndexPath].PathNodes[i]);
    //            assert(IndexV<PrioVect.size());
    //            Avg+=PrioVect[IndexV];
    //        }
    //        Avg/=ChoosenPaths[IndexPath].PathNodes.size();
    //        return Avg;
    //    }

    ScalarType GetPriority(CandidateTrace &CurrCand)
    {
        assert(PrioVect.size()==Mesh().vert.size());
        ScalarType Avg=0;
        for (size_t i=0;i<CurrCand.PathNodes.size();i++)
        {
            size_t IndexV=VertexFieldGraph<MeshType>::NodeVertI(CurrCand.PathNodes[i]);
            assert(IndexV<PrioVect.size());
            Avg+=PrioVect[IndexV];
        }
        Avg/=CurrCand.PathNodes.size();
        return Avg;
    }

    void SortPathByPriority(std::vector<CandidateTrace> &To_Sort)
    {
        std::cout<<"Sorting By Priority"<<std::endl;
        assert(PrioVect.size()==Mesh().vert.size());
        std::vector<std::pair<ScalarType,int> > PathPrio;
        for (size_t i=0;i<To_Sort.size();i++)
        {
            ScalarType currV=GetPriority(To_Sort[i]);
            PathPrio.push_back(std::pair<ScalarType,int>(currV,i));
        }

        std::sort(PathPrio.begin(),PathPrio.end());
        std::reverse(PathPrio.begin(),PathPrio.end());
        std::vector<CandidateTrace> SwapPaths;
        for (size_t i=0;i<PathPrio.size();i++)
        {
            int IndexP=PathPrio[i].second;
            assert(IndexP<(int)To_Sort.size());
            SwapPaths.push_back(To_Sort[IndexP]);
        }
        To_Sort=SwapPaths;
    }

    void RemovePaths()//bool DoSmooth=true)
    {
        if (PrioVect.size()>0)
            SortPathByPriority(ChoosenPaths);

        std::vector<std::vector<vcg::face::Pos<FaceType> > > PathPos;

        //select pos
        vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(Mesh());
        GetPathPos(VFGraph,ChoosenPaths,PathPos);
        Mesh().SelectPos(PathPos,true);

        if (DebugMsg)
            std::cout<<"Removing..."<<std::endl;

        //int s=0;
        while (RemoveIteration()){}

        //CutEarPath();

        //THIS SHOULD SPEED UP
        //UpdatePartitionsFromChoosen(true);
        UpdatePartitionsFromChoosen(false);

        //        ColorByPartitions();

        if (DebugMsg)
            WriteInfo();
    }

    void GetPatchMesh(const size_t &IndexPatch,
                      MeshType &PatchMesh,
                      bool InternalCuts)
    {

        PatchManager<MeshType>::GetMeshFromPatch(Mesh(),IndexPatch,Partitions,PatchMesh,InternalCuts);
    }

    void UpdateCornersWithInternalCuts()
    {
        if ((!AllowDarts)&&(!AllowSelfGluedPatch))return;

        for (size_t i=0;i<Partitions.size();i++)
        {
            std::vector<size_t> CornerValence;
            std::vector<size_t> NewCorners;
            PatchManager<MeshType>::GetCornerValuesFromInternalFeatures(Mesh(),Partitions[i],
                                                                        PartitionCorners[i],
                                                                        CornerValence);
            assert(CornerValence.size()==PartitionCorners[i].size());
            for (size_t j=0;j<CornerValence.size();j++)
            {
                size_t currValence=CornerValence[j];
                assert(currValence<=2);
                assert(currValence>=0);
                size_t currCorner=PartitionCorners[i][j];
                if (currValence==0)continue;
                //if the corner is splitted then count it twice
                NewCorners.resize(NewCorners.size()+currValence,currCorner);
                //push_back(PartitionCorners[i][j]);
                // if (CornerValence==1)
                //     NewCorners.push_back(PartitionCorners[i][j]);
            }
            PartitionCorners[i]=NewCorners;
        }
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

    void GetTraceableFlatNodes(std::vector<size_t> &TraceableNodes)
    {
        TraceableNodes.clear();
        for (size_t i=0;i<Traceable.size();i++)
        {
            TypeVert T=GetTypeOfNode(i);
            if (T!=TVFlat)continue;
            if (!Traceable[i])continue;
            TraceableNodes.push_back(i);
        }
    }

    void InitTraceableBorders()//size_t minsuBSteps=MIN_BORDER_SAMPLE,
    //size_t maxsuBSteps=MAX_BORDER_SAMPLE)
    {
        if (DebugMsg)
            std::cout<<"** INIT TRACEABLE BORDERS **"<<std::endl;

        std::vector<size_t> ConcaveV;
        GetVertexType(TVConcave,ConcaveV);

        std::vector<size_t> ConvexV;
        GetVertexType(TVConvex,ConvexV);

        std::vector<size_t> NarrowV;
        GetVertexType(TVNarrow,NarrowV);

        std::vector<size_t> Corners=ConcaveV;
        Corners.insert(Corners.end(),ConvexV.begin(),ConvexV.end());
        Corners.insert(Corners.end(),NarrowV.begin(),NarrowV.end());

        std::vector<std::vector<size_t> > BorderSequences;
        PatchManager<MeshType>::GetBorderSequences(Mesh(),Corners,BorderSequences);

        //set by defailt all Flat as non traceable
        std::vector<size_t> BorderFlat;
        GetNodesType(TVFlat,BorderFlat);
        for (size_t i=0;i<BorderFlat.size();i++)
        {
            size_t IndexN=BorderFlat[i];
            Traceable[IndexN]=false;
        }

        //for each non OK patch set the useful vertices
        std::vector<bool> UsefulV(Mesh().vert.size(),false);

        if (Partitions.size()==0)
            UsefulV=std::vector<bool>(Mesh().vert.size(),true);

        for (size_t i=0;i<Partitions.size();i++)
        {
            if (PartitionType[i]==IsOK)continue;

            for (size_t j=0;j<Partitions[i].size();j++)
            {
                FaceType *f=&Mesh().face[Partitions[i][j]];
                for (size_t k=0;k<3;k++)
                {
                    size_t IndexV=vcg::tri::Index(Mesh(),f->V(k));
                    UsefulV[IndexV]=true;
                }
            }
        }

        ScalarType totLen=0;
        for (size_t i=0;i<BorderSequences.size();i++)
        {
            ScalarType angle,len;
            PatchManager<MeshType>::GetSequencesLenghtAngle(Mesh(),BorderSequences[i],angle,len);
            totLen+=len;
        }
        ScalarType lenInt=totLen/1000;
        //std::cout<<"There are "<<BorderSequences.size()<<" border sequences"<<std::endl;
        for (size_t i=0;i<BorderSequences.size();i++)
        {
            ScalarType angle,len;
            PatchManager<MeshType>::GetSequencesLenghtAngle(Mesh(),BorderSequences[i],angle,len);
            //            std::cout<<"Size "<<BorderSequences[i].size()<<std::endl;
            ScalarType angleInt=((1.75*M_PI)/(MAX_BORDER_SAMPLE+1));
            //std::cout<<"Angle "<<angle<<std::endl;
            size_t DivisionAngle=floor(0.5+angle/angleInt);
            //ScalarType lenInt=Mesh().bbox.Diag()*0.1;//(Mesh().bbox.Diag()*0.25)/MAX_BORDER_SAMPLE;
            size_t DivisionLen=floor(0.5+len/lenInt);
            size_t Division=std::max(DivisionAngle,DivisionLen);
            //Division=std::min(Division,(size_t)MAX_BORDER_SAMPLE+1);
            Division=std::max(Division,(size_t)MIN_BORDER_SAMPLE+1);
            size_t Step=floor(0.5+((ScalarType)BorderSequences[i].size()/(ScalarType)Division));
            Step=std::max(Step,(size_t)1);
            for (size_t j=0;j<BorderSequences[i].size();j+=Step)
            {
                //std::cout<<BorderSequences[i][j]<<" ";
                size_t currVIndex=BorderSequences[i][j];
                if (VertType[currVIndex]!=TVFlat)continue;
                if (!UsefulV[currVIndex])continue;
                assert(currVIndex<Mesh().vert.size());
                std::vector<size_t> Nodes;
                VertexFieldGraph<MeshType>::IndexNodes(currVIndex,Nodes);
                for (size_t k=0;k<Nodes.size();k++)
                {
                    size_t IndexN=Nodes[k];
                    Traceable[IndexN]=true;
                }
            }
            //std::cout<<std::endl;
        }
        //        if (DebugMsg)
        //            std::cout<<"selected traceable borders"<<std::endl;
        if (DebugMsg)
            std::cout<<"** DONE **"<<std::endl;

    }

    void InitTracer(ScalarType _Drift,bool _DebugMsg)
    {
        vcg::tri::InitFaceIMark(Mesh());
        vcg::tri::InitVertexIMark(Mesh());

        DebugMsg=_DebugMsg;
        Drift=_Drift;

        InitInternalStructures();

        InitEdgeL();

        ChoosenPaths.clear();
        //ChoosenIsLoop.clear();
        MaxNarrowWeight=sqrt(TotArea(Mesh()))*MAX_NARROW_CONST*Drift;
        InitAvEdge();
        InitEdgeDirTable();

        UpdatePartitionsFromChoosen(true);
        ColorByPartitions();

        Traceable=std::vector<bool>(VFGraph.NumNodes(),true);
        InitTraceableBorders();
        //std::cout<<"****Convex Size "<<EDirTable.ConvexV.size()<<std::endl;

        EdgeSel0.clear();
        EdgeSel0.resize(Mesh().face.size(),std::vector<bool>(3,false));
        EdgeSel1.clear();
        EdgeSel1.resize(Mesh().face.size(),std::vector<bool>(3,false));

        //        //            if (InitializePart)
        //        //            {
        //                        UpdatePartitionsFromChoosen(true);
        //                        ColorByPartitions();
        //            //        }
    }

    size_t CopyPathsFrom(MyTracerType &Ptr,
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

    template <class PatchTracerType>
    void CopyFrom(PatchTracerType &Ptr,
                  std::vector<size_t> &VertMap,
                  size_t IndexPatch)
    {
        Drift=Ptr.Drift;
        DebugMsg=Ptr.DebugMsg;
        MaxNarrowWeight=Ptr.MaxNarrowWeight;
        avgEdge=Ptr.avgEdge;

        //0. basic checks
        assert(VertMap.size()==Mesh().vert.size());
        VertType.resize(Mesh().vert.size(),TVNone);
        VerticesNeeds.resize(Mesh().vert.size(),0);

        //1. copy the original vert type
        for (size_t i=0;i<VertMap.size();i++)
        {
            size_t IndexV=VertMap[i];
            assert(IndexV<Ptr.Mesh().vert.size());
            assert(IndexV<Ptr.VertType.size());
            VertType[i]=Ptr.VertType[IndexV];
            VerticesNeeds[i]=Ptr.VerticesNeeds[IndexV];
        }

        //2. update coherently for the new patch
        //in case of concave/narrow
        for (size_t i=0;i<VertType.size();i++)
        {
            if ((VertType[i]!=TVNarrow)&&
                    (VertType[i]!=TVConcave)&&
                    Mesh().vert[i].IsB())
                VertType[i]=TVFlat;

            if ((VertType[i]==TVNarrow)&&(VerticesNeeds[i]==0))
                VertType[i]=TVFlat;

            if ((VertType[i]==TVConcave)&&(VerticesNeeds[i]==0))
                VertType[i]=TVFlat;
        }

        //3. then also set the corners as convex
        std::set<size_t> CornerSet(Ptr.PartitionCorners[IndexPatch].begin(),
                                   Ptr.PartitionCorners[IndexPatch].end());
        for (size_t i=0;i<VertType.size();i++)
        {
            size_t IndexV=VertMap[i];
            if (CornerSet.count(IndexV)>0)
                VertType[i]=TVConvex;
        }

        //4. get configuration on borders
        std::vector<size_t> FlatEmitters,FlatReceivers,ChosenEmitters,
                ChosenReceivers,FlatTangent,ChosenTangent;


        //        std::vector<std::vector<CoordType> > VertFlatDir,VertOrthoDir;
        //        VertexEmitter<MeshType>::GetOrthoFlatDirections(VFGraph,VertFlatDir,VertOrthoDir);

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

        InvalidateNonConcaveSing();
    }

    void GetPatchNodes(const size_t &IndexPatch,
                       std::vector<size_t> &PatchNodes)
    {
        //        for (size_t i=0;i<Partitions[IndexPatch].size();i++)
        //        {
        //            size_t FaceI=Partitions[IndexPatch][i];
        //            for (size_t j=0;j<3;j++)
        //            {
        //                size_t VertI=vcg::tri::Index(Mesh(),Mesh().face[FaceI].V(j));
        //                std::vector<size_t> NodeI;
        //                VertexFieldGraph<MeshType>::IndexNodes(VertI,NodeI);
        //                PatchNodes.insert(PatchNodes.end(),NodeI.begin(),NodeI.end());
        //            }
        //        }
        //        std::sort(PatchNodes.begin(),PatchNodes.end());
        //        auto last=std::unique(PatchNodes.begin(),PatchNodes.end());
        //        PatchNodes.erase(last, PatchNodes.end());

        for (size_t i=0;i<Mesh().vert.size();i++)
            vcg::tri::Mark(Mesh(),&Mesh().vert[i]);

        vcg::tri::UnMarkAll(Mesh());
        for (size_t i=0;i<Partitions[IndexPatch].size();i++)
        {
            size_t FaceI=Partitions[IndexPatch][i];
            for (size_t j=0;j<3;j++)
            {
                VertexType *v=Mesh().face[FaceI].V(j);
                size_t VertI=vcg::tri::Index(Mesh(),v);
                if (vcg::tri::IsMarked(Mesh(),v))continue;
                vcg::tri::Mark(Mesh(),v);
                std::vector<size_t> NodeI;
                VertexFieldGraph<MeshType>::IndexNodes(VertI,NodeI);
                PatchNodes.insert(PatchNodes.end(),NodeI.begin(),NodeI.end());
            }
        }
        //        std::sort(PatchNodes.begin(),PatchNodes.end());
        //        auto last=std::unique(PatchNodes.begin(),PatchNodes.end());
        //        PatchNodes.erase(last, PatchNodes.end());

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
        DiscardedCandidates.clear();

        if (DebugMsg)
            std::cout<<"Adding candidates (Trace From)"<<std::endl;

        assert(Traceable.size()==VFGraph.NumNodes());
        for (size_t i=0;i<VFGraph.NumNodes();i++)
        {
            //not the same kind
            if (!CanEmit[i])continue;
            if (!Traceable[i])continue;

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
        if (DebugMsg)
            std::cout<<"Aftern Expansion there are "<<Candidates.size()<<" candidates"<<std::endl;

        if (Candidates.size()==0)return false;

        int size0=ChoosenPaths.size();

        //        if (DebugMsg)
        //            std::cout<<"After Expansion there are "<<Candidates.size()<<" candidates"<<std::endl;


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

        //        size_t t0=clock();
        std::vector<size_t> PatchNodes;
        GetPatchNodes(IndexPatch,PatchNodes);

        //        size_t t1=clock();

        std::set<size_t> PatchNodesSet(PatchNodes.begin(),PatchNodes.end());

        //GetFlatEmitters(FlatEmitters);
        //GetFlatReceivers(FlatReceivers);

        //        size_t t2=clock();

        GetActiveEmittersType(TVFlat,FlatEmitters);
        GetActiveReceiversType(TVFlat,FlatReceivers);

        //        size_t t3=clock();

        //just filters the one in the patch
        FlatEmitters=FilterFromSet(FlatEmitters,PatchNodesSet);
        FlatReceivers=FilterFromSet(FlatReceivers,PatchNodesSet);

        GetFlatTangentNodes(FlatTangent);
        FlatTangent=FilterFromSet(FlatTangent,PatchNodesSet);

        //        size_t t4=clock();

        //then get the one from choosen
        std::vector<size_t> ChoosenNodes;
        //GetChoosenNodes(ChoosenNodes);
        GetCandidateNodes(ChoosenPaths,ChoosenNodes);
        ChoosenNodes=FilterFromSet(ChoosenNodes,PatchNodesSet);

        //        size_t t5=clock();

        //GetChoosenNodesAndTangent(ChosenTangent);
        GetCandidateNodesNodesAndTangent<MeshType>(ChoosenPaths,ChosenTangent);
        ChosenTangent=FilterFromSet(ChosenTangent,PatchNodesSet);

        //      size_t t6=clock();

        //get the ortho direction of the patch
        MeshType patchMesh;
        GetPatchMesh(IndexPatch,patchMesh,false);
        std::vector<std::vector<CoordType> > VertFlatDir,VertOrthoDir;
        VertexEmitter<MeshType>::GetOrthoFlatDirections(patchMesh,VertFlatDir,VertOrthoDir);

        //        size_t t7=clock();

        //then average ortho dirs
        std::vector<CoordType> AvOrthoDir(VertOrthoDir.size(),CoordType(0,0,0));
        for (size_t i=0;i<VertOrthoDir.size();i++)
        {
            for (size_t j=0;j<VertOrthoDir[i].size();j++)
                AvOrthoDir[i]+=VertOrthoDir[i][j];
            if (AvOrthoDir[i]!=CoordType(0,0,0))
                AvOrthoDir[i].Normalize();
        }

        //create a map to retrieve original vert
        std::map<size_t,size_t> VertMap;
        for (size_t i=0;i<patchMesh.vert.size();i++)
        {
            int IndexV=patchMesh.vert[i].Q();
            assert(IndexV>=0);
            assert(IndexV<(int)Mesh().vert.size());
            VertMap[IndexV]=i;
        }


        //then for each one get the ortho nodes
        for (size_t i=0;i<ChoosenNodes.size();i++)
        {
            size_t OrthoN0,OrthoN1;
            VertexFieldGraph<MeshType>::OrthoNode(ChoosenNodes[i],OrthoN0,OrthoN1);
            size_t IndexV=VertexFieldGraph<MeshType>::NodeVertI(ChoosenNodes[i]);
            assert(OrthoN0!=OrthoN1);
            CoordType Dir0=VFGraph.NodeDir(OrthoN0);
            CoordType Dir1=VFGraph.NodeDir(OrthoN1);
            //then get the target ortho one
            assert(VertMap.count(IndexV)>0);
            //remap to the submesh
            IndexV=VertMap[IndexV];
            assert(IndexV>=0);
            assert(IndexV<AvOrthoDir.size());
            CoordType TargetD=AvOrthoDir[IndexV];

            //in this case it is an internal edge that has been reinserted
            //then add both
            if (TargetD==CoordType(0,0,0))
            {
                ChosenEmitters.push_back(OrthoN0);
                ChosenReceivers.push_back(OrthoN0);
                ChosenEmitters.push_back(OrthoN1);
                ChosenReceivers.push_back(OrthoN1);
                continue;
            }
            ScalarType dot0=Dir0*TargetD;
            ScalarType dot1=Dir1*TargetD;
            if (dot0>dot1)
            {
                ChosenEmitters.push_back(OrthoN0);
                ChosenReceivers.push_back(OrthoN1);
            }
            else
            {
                ChosenEmitters.push_back(OrthoN1);
                ChosenReceivers.push_back(OrthoN0);
            }
            //            std::vector<size_t> NodeNeigh0,NodeNeigh1;

            //            VFGraph.GetNodeNeigh(OrthoN0,NodeNeigh0);
            //            NodeNeigh0=FilterFromSet(NodeNeigh0,PatchNodesSet);

            //            VFGraph.GetNodeNeigh(OrthoN1,NodeNeigh1);
            //            NodeNeigh1=FilterFromSet(NodeNeigh1,PatchNodesSet);

            //            if (NodeNeigh0.size()>NodeNeigh1.size())
            //            {
            //                ChosenEmitters.push_back(OrthoN0);
            //                ChosenReceivers.push_back(OrthoN1);
            //            }
            //            else
            //            {
            //                ChosenEmitters.push_back(OrthoN1);
            //                ChosenReceivers.push_back(OrthoN0);
            //            }
        }

        //        size_t t8=clock();

        //        Time_InitSubPatches2_0+=t1-t0;
        //        Time_InitSubPatches2_1+=t2-t1;
        //        Time_InitSubPatches2_2+=t3-t2;
        //        Time_InitSubPatches2_3+=t4-t3;
        //        Time_InitSubPatches2_4+=t5-t4;
        //        Time_InitSubPatches2_5+=t6-t5;
        //        Time_InitSubPatches2_6+=t7-t6;
        //        Time_InitSubPatches2_7+=t8-t7;


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
        PatchManager<MeshType>::ColorByVarianceLenght(Mesh(),Partitions,PartitionCorners,EdgeL);
    }

    void ColorByLenghtDistortion()
    {
        PatchManager<MeshType>::ColorByDistortionLenght(Mesh(),Partitions,PartitionCorners,EdgeL);
    }

    void ColorByArapDistortion()
    {
        PatchManager<MeshType>::ColorByUVDistortionFaces(Mesh(),Partitions,PartitionCorners,Arap,false,false);
    }

    void ColorByCClarkability()
    {
        InitEdgeL();
        PatchManager<MeshType>::ColorByCatmullClarkability(Mesh(),Partitions,PartitionCorners,
                                                           EdgeL,CClarkability,avgEdge,match_valence);
    }

    size_t UnsatisfiedNum()
    {
        size_t UnsatisfiedNum=0;
        for (size_t i=0;i<VerticesNeeds.size();i++)
            UnsatisfiedNum+=VerticesNeeds[i];
        return UnsatisfiedNum;
    }

    void JoinNarrowStep()//bool UpdatePartition=true)
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

        //        if(UpdatePartition)
        //        {
        //            UpdatePartitionsFromChoosen();
        //            ColorByPartitions();
        //        }
    }


    void JoinConcaveStep()//bool UpdatePartition=true)
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

            if (DebugMsg)
                std::cout<<"Still "<<UnsatisfiedNum()<<" Non Connected"<<std::endl;
        }
        while (Joined);
        //JoinConnection(Concave,Flat,DijkstraReceivers);

        size_t NumPath1=ChoosenPaths.size();
        if (NumPath1==NumPath0)return;

        //        if(UpdatePartition)
        //        {
        //            UpdatePartitionsFromChoosen();
        //            ColorByPartitions();
        //        }
    }


    void JoinBoundaries(bool UsePartitionNeeds)
    {
        if (DebugMsg)
            std::cout<<"**TRACING BORDERS ***"<<std::endl;

        bool Joined=true;
        do{
            Joined=false;
            if (DebugMsg)
                std::cout<<"*Trace Direct"<<std::endl;

            Joined|=JoinConnection(TVFlat,TVFlat,TraceDirect,UsePartitionNeeds);

            if (DebugMsg)
                std::cout<<"*Trace Dijkstra"<<std::endl;

            Joined|=JoinConnection(TVFlat,TVFlat,DijkstraReceivers,UsePartitionNeeds);
            //            InitCandidates(Flat,Flat,TraceDirect);
            //            InitCandidates(Flat,Flat,DijkstraReceivers);

            //            ChooseGreedyByDistance(false,UsePartitionNeeds);
        }while (Joined);
        //std::cout<<""
        //exit(0);
        //        if(UpdatePartition)
        //        {
        //            UpdatePartitionsFromChoosen();
        //            ColorByPartitions();
        //        }
    }

    void TraceBorderAndLoops(bool UsePartitionNeeds)
    {
        InitCandidates(TVInternal,TVInternal,TraceLoop);
        //std::cout<<"TEST 0 There are "<<Candidates.size()<<"candidates"<<std::endl;

        InitCandidates(TVFlat,TVFlat,TraceDirect);
        //std::cout<<"TEST 1 There are "<<Candidates.size()<<"candidates"<<std::endl;

        //UpdatePartitionsFromChoosen(true);

        ChooseGreedyByDistance(false,UsePartitionNeeds);

        //std::cout<<"Choosen "<<ChoosenPaths.size()<<"candidates"<<std::endl;

    }

    size_t TraceLoops(bool UsePartitionNeeds)
    {

        if (DebugMsg)
            std::cout<<"**TRACING LOOPS ***"<<std::endl;
        //ensure a minimal of loop to be traced, might happens when there is
        //recoursive trace
        //        std::vector<size_t> NodesSet;
        //        GetEmitterType(TVInternal,NodesSet);
        //        GetE
        //        if (NodesSet.size()<MIN_SAMPLES_HARD)
        //            SampleLoopEmitters();

        InitCandidates(TVInternal,TVInternal,TraceLoop);
        size_t Size0=ChoosenPaths.size();
        ChooseGreedyByDistance(false,UsePartitionNeeds);
        //        if(UpdatePartition)
        //        {
        //            UpdatePartitionsFromChoosen();
        //            ColorByPartitions();
        //        }
        size_t Size1=ChoosenPaths.size();
        assert(Size1>=Size0);
        return (Size1-Size0);
    }

    size_t BatchRemovalOnMesh(bool PreRemoveStep=true)//bool do_smooth=true)
    {
        //        size_t Size0=ChoosenPaths.size();
        //        RemovePaths();
        //        size_t Size1=ChoosenPaths.size();
        //        return (Size0-Size1);


        size_t Size0=ChoosenPaths.size();
        if (PreRemoveStep)
            RemovePaths();//false);
        size_t Size1=ChoosenPaths.size();

        if (split_on_removal)
            SplitIntoSubPaths();

        size_t Size2=ChoosenPaths.size();
        RemovePaths();

        if (!PreRemoveStep)
        {
            SplitIntoSubPaths();
            RemovePaths();// to be checked
        }

        if (DebugMsg)
            WriteInfo();

        size_t Size3=ChoosenPaths.size();
        assert(Size1<=Size0);
        assert(Size3<=Size2);
        return ((Size0-Size1)+(Size2-Size3));
    }

    void SetPriorityVect(const std::vector<ScalarType> &Vect)
    {
        assert((Vect.size()==0)||(Vect.size()==Mesh().vert.size()));
        PrioVect=Vect;
    }

    void ClearPriorityVect()
    {
        PrioVect.clear();
    }

//    void RemoveDarts0()//bool do_smooth=true)
//    {
//        std::cout<<"REMOVE DARTS"<<std::endl;
//        std::vector<std::vector<vcg::face::Pos<FaceType> > > PathPos;

//        //select pos
//        //std::cout<<"a"<<std::endl;
//        vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(Mesh());
//        GetPathPos(VFGraph,ChoosenPaths,PathPos);
//        Mesh().SelectPos(PathPos,true);

//        if (DebugMsg)
//            std::cout<<"Removing..."<<std::endl;

//        //only makes sense if quality is on (usually distortion)
//        assert(check_quality_functor);


//        size_t NumRem=0;
//        do{
//            MergeContiguousPaths();
//            SplitIntoSubPaths();
//            SplitIntoIntervals(ChoosenPaths);

//            int size0=ChoosenPaths.size();

//            RemovePaths();
//            RemoveEmptyPaths();

//            int size1=ChoosenPaths.size();
//            NumRem=size0-size1;
//        }while (NumRem>0);

//        std::cout<<"END REMOVE DARTS"<<std::endl;
//        //        if (HasRemoved)
//        //        {
//        //            RemoveEmptyPaths();
//        //        }

//    }

    void RemoveDarts()//bool do_smooth=true)
    {
        UpdatePartitionsFromChoosen(true);

        std::cout<<"REMOVE DARTS"<<std::endl;
        std::vector<std::vector<vcg::face::Pos<FaceType> > > PathPos;

        //select pos
        vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(Mesh());
        GetPathPos(VFGraph,ChoosenPaths,PathPos);
        Mesh().SelectPos(PathPos,true);

        //split in subpaths

        bool Has_changed=false;
        do{
            RemoveEmptyPaths();
            MergeContiguousPaths();
            SplitIntoSubPaths();

            //            vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(Mesh());
            //            GetPathPos(VFGraph,ChoosenPaths,PathPos);
            //            Mesh().SelectPos(PathPos,true);

            //sort by priority globally
            if (PrioVect.size()>0)
                SortPathByPriority(ChoosenPaths);
            if (DebugMsg)
                std::cout<<"Removing..."<<std::endl;

            //only makes sense if quality is on (usually distortion)
            assert(check_quality_functor);
            Has_changed=false;
            int InitialSize=ChoosenPaths.size();
            for (int i=InitialSize-1;i>=0;i--)
            {
                //get the current Path
                CandidateTrace CurrP=ChoosenPaths[i];

                //save the path on the mesh
                std::vector<vcg::face::Pos<FaceType> > CurrPos;
                VFGraph.GetNodesPos(CurrP.PathNodes,CurrP.IsLoop,CurrPos);

                CurrP.Unremovable=false;
                std::vector<CandidateTrace> To_split;
                //get the old patch index
                To_split.push_back(CurrP);

                //then split it
                SplitIntoIntervals(To_split);

                //sort by priority
                SortPathByPriority(To_split);

                //substitute paths
                ChoosenPaths[i].PathNodes.clear();

                //set the old one as non removeable
                for (size_t j=0;j<ChoosenPaths.size();j++)
                    ChoosenPaths[j].Unremovable=true;

                //add the splitted one the vector
                size_t num0=ChoosenPaths.size();
                for (size_t j=0;j<To_split.size();j++)
                {
                    ChoosenPaths.push_back(To_split[j]);
                    ChoosenPaths.back().Unremovable=false;
                }

                //then remove if possible the last substitution

                bool Has_Rem=false;
                bool Has_Removed_Once=false;
                ContinuosCheckUVInt=false;
                do{
                    Has_Rem=false;
                    for (int j=(int)ChoosenPaths.size()-1;j>=num0;j--)
                    {
                        std::cout<<"Remove Attempt"<<std::endl;
                        Has_Rem|=RemoveIfPossible(j);
                    }
                    Has_Removed_Once|=Has_Rem;

                }while (Has_Rem);

                if (!Has_Removed_Once)continue;

                if (!CheckUVIntersection) continue;

                MeshType SubMesh;
                PatchManager<MeshType>::GetSubMeshSurroundingUsingEdgeSel(Mesh(),CurrPos,SubMesh);

                //then apply the parametrizator
                PatchQualityFunctor PFunct;
                PFunct.ContinuousCheckSelfInt()=false;

                //then make last check
                ScalarType ValQ=PFunct(SubMesh);

                //check the self Intersection
                bool SelfInt=PatchManager<MeshType>::SelfOverlapUV(SubMesh);
                if (SelfInt)
                {
                    std::cout<<"RESTORED FOR SELF INTERSECTION"<<std::endl;
                    //restore the old one
                    ChoosenPaths[i]=CurrP;
                    Mesh().SelectPos(CurrPos,true);
                    //first remove the splitted

                    AddEdgeNodes<MeshType>(CurrP.PathNodes,CurrP.IsLoop,EDirTable);
                    //cancel the split ones
                    for (int j=((int)ChoosenPaths.size()-1);j>=num0;j--)
                        ChoosenPaths[j].PathNodes.clear();
                }
                else
                {
                    //in this case simply accept the modification
                    Has_changed=true;
                }
            }

        }while (Has_changed);
        ContinuosCheckUVInt=true;
        RemoveEmptyPaths();
        MergeContiguousPaths();
        std::cout<<"LAST REMOVAL OP"<<std::endl;
        UpdatePartitionsFromChoosen(true);
        WriteInfo();
        std::cout<<"END REMOVE DARTS"<<std::endl;


    }

    size_t BatchRemovalMetaMesh(bool PreRemoveStep=true)//bool do_smooth=true)
    {
        size_t t0,t1;//t2;
        size_t Size0=ChoosenPaths.size();
        if (PreRemoveStep)
        {
            if (DebugMsg)
                std::cout<<"**Pre Removal Step**"<<std::endl;
            t0=clock();
            InitMetaMesh();
            t1=clock();
            Time_InitMetaMesh+=t1-t0;

            RemoveMetaMeshStep();

        }
        size_t Size1=ChoosenPaths.size();

        if (split_on_removal)
            SplitIntoSubPaths();

        size_t Size2=ChoosenPaths.size();

        if (DebugMsg)
            std::cout<<"**Regular Removal Step**"<<std::endl;

        t0=clock();
        InitMetaMesh();
        t1=clock();
        Time_InitMetaMesh+=t1-t0;

        RemoveMetaMeshStep();

        if (!PreRemoveStep)
        {
            if (DebugMsg)
                std::cout<<"**Regular Removal Step**"<<std::endl;

            SplitIntoSubPaths();
            t0=clock();
            InitMetaMesh();
            t1=clock();
            Time_InitMetaMesh+=t1-t0;

            RemoveMetaMeshStep();

        }

        if (DebugMsg)
            WriteInfo();


        size_t Size3=ChoosenPaths.size();
        assert(Size1<=Size0);
        assert(Size3<=Size2);
        return ((Size0-Size1)+(Size2-Size3));
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
                angleSumH[IndexV] += vcg::face::WedgeAngleRad(TestMesh.face[i],j);
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

    void RemoveConvexCorners(std::vector<size_t> &Corners,size_t numRemoval)
    {
        assert(numRemoval>0);
        int to_remove=-1;
        size_t removed=0;
        do{
            to_remove=-1;
            ScalarType minAngle=0;
            for (size_t i=0;i<Corners.size();i++)
            {
                size_t IndexV=Corners[i];
                if (VertType[IndexV]!=TVConvex)continue;
                if (Mesh().vert[IndexV].Q()<minAngle)continue;
                to_remove=i;
                minAngle=Mesh().vert[IndexV].Q();
            }
            if (to_remove==-1)return;
            Corners.erase(Corners.begin()+to_remove);
            removed++;
        }while (removed<numRemoval);
    }

    size_t FixAmmittableValences()
    {
        if (!CheckQuadrangulationLimits) return 0;

        size_t NeedFix=0;
        for (size_t i=0;i<PartitionCorners.size();i++)
            if ((PartitionCorners[i].size()<MIN_ADMITTIBLE)||
                    (PartitionCorners[i].size()>MAX_ADMITTIBLE))
            {
                NeedFix++;
                FixCorners(i);
            }
        return NeedFix;
    }

    void FixValences()
    {

        size_t NeedFix=FixAmmittableValences();

        if (match_valence)
        {
            vcg::tri::UpdateQuality<MeshType>::VertexConstant(Mesh(),0);
            for(size_t i=0;i<Mesh().face.size();i++)
            {
                if(Mesh().face[i].IsD())continue;
                for(int j=0;j<3;j++)
                    Mesh().face[i].V(j)->Q()+= vcg::face::WedgeAngleRad(Mesh().face[i],j);
            }

            //solve valence 4 not valid because of corners
            for (size_t i=0;i<PartitionCorners.size();i++)
            {
                bool SingOnCorner=false;
                int ExpVal=PatchManager<MeshType>::ExpectedValence(Mesh(),Partitions[i],PartitionCorners[i],SingOnCorner);
                if (ExpVal!=4)continue;
                if (SingOnCorner)continue;
                if (ExpVal>(int)PartitionCorners[i].size())continue;
                if (PartitionCorners[i].size()!=4)
                {
                    NeedFix++;
                    RemoveConvexCorners(PartitionCorners[i],PartitionCorners[i].size()-4);
                }
            }
        }
        if (DebugMsg)
            std::cout<<"FINAL Fixed "<<NeedFix<<std::endl;
    }


    void BatchAddLoops(bool ForceReceivers,
                       bool AddOnlyNeeded,
                       bool OnlyNarrowConcave,
                       bool force_always)
    {
        //might need to resample
        //other disconnected components
        //std::vector<size_t> AllEmit;
        //GetAllEmitter(AllEmit);
        //        if (AllEmit.size()<MIN_SAMPLES_HARD)
        //        {
        ////            size_t genus=vcg::tri::Clean<MeshType>::MeshGenus(Mesh());
        //            size_t numB=vcg::tri::Clean<MeshType>::CountHoles(Mesh());
        //            if (numB==0)
        //                SampleLoopEmitters();
        //        }

        if (AddOnlyNeeded)
            UpdatePartitionsFromChoosen(true);

        if (DebugMsg)
            std::cout<<"**TRACING NARROW/CONCAVE ***"<<std::endl;

        size_t numInitial=ChoosenPaths.size();
        if (ForceReceivers)
        {
            AllReceivers=true;
            MaxNarrowWeight/=100;
            JoinNarrowStep();
            JoinConcaveStep();
            if (!force_always)
            {
                MaxNarrowWeight*=100;
                AllReceivers=false;
            }
        }
        JoinNarrowStep();
        JoinConcaveStep();

        if (DebugMsg)
            std::cout<<"done"<<std::endl;

        //then update the partitions
        size_t numAddedNarConc=ChoosenPaths.size();
        if (AddOnlyNeeded && (numInitial!=numAddedNarConc))
        {
            //            LazyUpdatePartitions();
            //            LazyUpdatePartitions();
            UpdatePartitionsFromChoosen(true);
        }

        Candidates.clear();
        if (OnlyNarrowConcave)
            return;

        //        if (DebugMsg)
        if (PrioMode==PrioModeLoop)
        {
            if (DebugMsg)
                std::cout<<"**TRACING LOOPS ***"<<std::endl;
            size_t numAddedLoops=TraceLoops(AddOnlyNeeded);

            if (DebugMsg)
                std::cout<<"done"<<std::endl;

            if (DebugMsg)
                std::cout<<"Num Chosen "<<ChoosenPaths.size()<<std::endl;

            if ((AddOnlyNeeded)&&(numAddedLoops>0))//||(numRemoved>0)))
                //THIS SHOULD SPEED UP
                LazyUpdatePartitions();
            //UpdatePartitionsFromChoosen(true);

            if (DebugMsg)
                std::cout<<"**TRACING BORDERS ***"<<std::endl;
            JoinBoundaries(AddOnlyNeeded);
            if (DebugMsg)
                std::cout<<"done"<<std::endl;
        }

        if (PrioMode==PrioModBorder)
        {

            if (DebugMsg)
                std::cout<<"**TRACING BORDERS ***"<<std::endl;
            JoinBoundaries(AddOnlyNeeded);
            if (DebugMsg)
                std::cout<<"done"<<std::endl;

            if (AddOnlyNeeded)//&&(numAddedLoops>0))//||(numRemoved>0)))
                UpdatePartitionsFromChoosen(true);

            if (DebugMsg)
                std::cout<<"**TRACING LOOPS ***"<<std::endl;
            TraceLoops(AddOnlyNeeded);

            if (DebugMsg)
                std::cout<<"done"<<std::endl;

            if (DebugMsg)
                std::cout<<"Num Chosen "<<ChoosenPaths.size()<<std::endl;

        }

        if (PrioMode==PrioModBlend)
        {
            TraceBorderAndLoops(AddOnlyNeeded);

            if (DebugMsg)
                std::cout<<"done"<<std::endl;

            if (DebugMsg)
                std::cout<<"Num Chosen "<<ChoosenPaths.size()<<std::endl;

        }

        //restore the setup
        if (ForceReceivers && force_always)
        {
            MaxNarrowWeight*=100;
            AllReceivers=false;
        }
        //        if (FinalRemoval)
        //            BatchRemoval();//SmoothOnRemoval);
    }



    void InitCandidates(TypeVert FromType,
                        TypeVert ToType,
                        TraceType TrType)
    {

        std::vector<bool> CanEmit,CanReceive,MustDisable;
        GetTracingConfiguration(FromType,ToType,TrType,CanEmit,CanReceive,MustDisable);

        //DiscardedCandidates.clear();

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

    void  GetCurrDiscarded(std::vector<std::vector<size_t> > &CurrDiscarded)
    {
        CurrDiscarded.clear();
        for (size_t i=0;i<DiscardedCandidates.size();i++)
            CurrDiscarded.push_back(DiscardedCandidates[i].PathNodes);
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
            std::vector<size_t> NodesI;
            VertexFieldGraph<MeshType>::IndexNodes(i,NodesI);
            for (size_t j=0;j<NodesI.size();j++)
            {
                //                if (!VFGraph.IsActive(NodesI[j]))continue;
                //                if((NodeEmitterTypes[NodesI[j]]==TVNarrow)||
                //                        (NodeEmitterTypes[NodesI[j]]==TVConcave))
                //                {
                Remaining.push_back(NodesI[j]);
                //                }
            }
        }
    }

    void GetUnsolvedPartitions(std::vector<std::vector<size_t> > &UnsolvedPartition,
                               std::vector<PatchType> &UnsolvedType,
                               bool UpdatePart)
    {
        UnsolvedPartition.clear();
        UnsolvedType.clear();

        if (UpdatePart)
            UpdatePartitionsFromChoosen(true);

        for (size_t i=0;i<Partitions.size();i++)
        {
            if (PartitionType[i]==IsOK)continue;
            UnsolvedPartition.push_back(Partitions[i]);
            UnsolvedType.push_back(PartitionType[i]);
        }
    }

    void GetUnsolvedMesh(MeshType &problem_mesh)
    {
        problem_mesh.Clear();
        std::vector<std::vector<size_t> > UnsolvedPartition;
        std::vector<PatchType> UnsolvedType;
        GetUnsolvedPartitions(UnsolvedPartition,UnsolvedType,false);
        PatchManager<MeshType>::SelectFaces(Mesh(),UnsolvedPartition);
        //        for (size_t i=0;i<UnsolvedType.size();i++)
        //        {
        //            if (UnsolvedType[i]!=NonDisk)continue;
        //            PatchManager<MeshType>::SelectFaces(Mesh(),UnsolvedPartition[i]);
        //        }
        vcg::tri::UpdateSelection<MeshType>::VertexFromFaceLoose(problem_mesh);

        vcg::tri::Append<MeshType,MeshType>::Mesh(problem_mesh,Mesh(),true);
        //       vcg::tri::io::ExporterPLY<MeshType>::Save(problem_mesh,"test_problem.ply");
        problem_mesh.UpdateAttributes();
    }

    void GetCurrChosenIsLoop(std::vector<bool> &ChosenIsLoop)
    {
        ChosenIsLoop.clear();
        GetIsLoop(ChoosenPaths,ChosenIsLoop);
    }

    void GetCurrDiscardedIsLoop(std::vector<bool> &DiscardedIsLoop)
    {
        DiscardedIsLoop.clear();
        GetIsLoop(DiscardedCandidates,DiscardedIsLoop);
    }

    //    void CutEarPath()
    //    {
    //        PatchManager<MeshType>::SelectMeshPatchBorders(Mesh(),FacePartitions);
    //        //NEED TO SELECT PATCH END POINTS !!!
    //        for (size_t i=0;i<ChoosenPaths.size();i++)
    //        {
    //            VertexFieldQuery<MeshType>::CutEarPath(VFGraph,ChoosenPaths[i].PathNodes,ChoosenPaths[i].IsLoop);
    //        }

    //        InitEdgeDirTable();
    //        InitEdgeL();
    //    }


    void MergeContiguousPaths(std::vector<CandidateTrace> &TraceSet)
    {
        TestedMerged=0;
        DoneMerged=0;
        while(MergeContiguousPathStep(TraceSet))
        {
            //std::cout<<"Size Merging Set"<<TraceSet.size()<<std::endl;
            RemoveEmptyPaths();
        }
        //        std::cout<<"Tested Merged "<<TestedMerged<<std::endl;
        //        std::cout<<"Done Merged "<<DoneMerged<<std::endl;
    }

    bool SingleRemoveStep(bool singleStep)//bool printNonRemoved=false)
    {
        std::vector<std::vector<vcg::face::Pos<FaceType> > > PathPos;

        //select pos
        vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(Mesh());
        GetPathPos(VFGraph,ChoosenPaths,PathPos);
        Mesh().SelectPos(PathPos,true);

        if (DebugMsg)
            std::cout<<"Removing..."<<std::endl;

        //int s=0;
        bool removed=RemoveIteration(singleStep);//true);

        if (removed)
            std::cout<<"** Successfully Removed, Num Path:"<<ChoosenPaths.size()<<std::endl;
        else
            std::cout<<"** NOT Removed, Num Path:"<<ChoosenPaths.size()<<std::endl;
        //CutEarPath();

        UpdatePartitionsFromChoosen(true);

        //        ColorByPartitions();

        if (DebugMsg)
            WriteInfo();

        if (!removed)
        {
            std::cout<<"** NON REMOVED BECAUSE **"<<std::endl;
            std::cout<<"* Zero Size:"<<RMZeroSize<<std::endl;
            std::cout<<"* Unremoveable:"<<RMUnremoveable<<std::endl;
            std::cout<<"* ConcaveNarrowGen:"<<RMHasConcaveNarrow<<std::endl;
            std::cout<<"* TJunctionCreation:"<<RMHasTJunction<<std::endl;
            std::cout<<"* HasDeadEnd:"<<RMHasDeadEnd<<std::endl;
            std::cout<<"* OutOfCornerNum:"<<RMOutOfCornerNum<<std::endl;
            std::cout<<"* NonProfitable:"<<RMNotProfitable<<std::endl;
        }
        return removed;
    }

    void GetEdgeNodes(const size_t &IndexV0,const size_t &IndexV1,
                      size_t &IndexN0,size_t &IndexN1)const
    {
        VFGraph.GetEdgeNodes(IndexV0,IndexV1,IndexN0,IndexN1);
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

    void GetTopologicallyNotOKPartitionsIndex(std::vector<size_t > &NonGenusPartitionIndex)
    {
        NonGenusPartitionIndex.clear();

        for (size_t i=0;i<Partitions.size();i++)
        {
            if (PatchManager<MeshType>::PatchGenus(Mesh(),Partitions[i],AllowDarts,AllowSelfGluedPatch)!=1)
                NonGenusPartitionIndex.push_back(i);
        }
    }

    void GetTopologicallyNotOKPartitionsIndex(std::set<size_t > &NonGenusPartitionIndexSet)
    {
        std::vector<size_t > NonGenusPartitionIndexV;
        GetTopologicallyNotOKPartitionsIndex(NonGenusPartitionIndexV);
        NonGenusPartitionIndexSet=std::set<size_t >(NonGenusPartitionIndexV.begin(),NonGenusPartitionIndexV.end());
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
        //UpdatePartitionsFromChoosen();
        //ColorByPartitions();
        InitEdgeL();
    }

    void ReinitPathFromEdgeSel()
    {
        ChoosenPaths.clear();
        std::vector<std::vector<PosType > > TestPosSeq;
        RetrievePosSeqFromSelEdges(Mesh(),TestPosSeq);
        for (size_t i=0;i<TestPosSeq.size();i++)
        {
            std::vector<size_t> NodeList;
            assert(TestPosSeq[i].size()>0);

            VertexType *v=TestPosSeq[i][0].VFlip();
            size_t IndexV=vcg::tri::Index(Mesh(),v);
            NodeList.push_back(IndexV);
            for (size_t j=0;j<TestPosSeq[i].size();j++)
            {
                v=TestPosSeq[i][j].V();
                IndexV=vcg::tri::Index(Mesh(),v);
                NodeList.push_back(IndexV);
            }

            bool IsLoop=false;
            IsLoop=(NodeList[0]==NodeList.back());
            if (IsLoop)NodeList.pop_back();

            assert(NodeList.size()>=2);
            //ChoosenPaths.resize(ChoosenPaths.size()+1);

            //set some dfault values
            CandidateTrace CTrace;
            CTrace.FromType=TVNone;
            CTrace.ToType=TVNone;
            CTrace.TracingMethod=TraceDirect;
            CTrace.IsLoop=IsLoop;
            CTrace.Updated=true;
            CTrace.Priority=0;
            CTrace.Unremovable=false;
            size_t IndexN0,IndexN1;
            for (size_t j=0;j<NodeList.size()-1;j++)
            {
                VFGraph.GetEdgeNodes(NodeList[j],NodeList[j+1],IndexN0,IndexN1);
                CTrace.PathNodes.push_back(IndexN0);
            }
            CTrace.PathNodes.push_back(IndexN1);
            CTrace.InitNode=CTrace.PathNodes[0];
            AddChoosen(CTrace);
        }
        UpdatePartitionsFromChoosen(true);

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
        PatchManager<MeshType>::ParametrizePatches(Mesh(),splittedUV,Partitions,PartitionCorners,Arap,
                                                   false,true,true,(AllowDarts||AllowSelfGluedPatch));
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

        PatchManager<MeshType>::DerivePerFacePartition(Mesh(),Partitions,FacePartitions);

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

    void Reset()
    {
        Partitions.clear();
        FacePartitions.clear();
        PartitionType.clear();
        PartitionCorners.clear();
        PatchInfos.clear();
        VertType.clear();
        NodeEmitterTypes.clear();
        NodeReceiverTypes.clear();
        CurrNodeDist.clear();
        Traceable.clear();
        Candidates.clear();
        ChoosenPaths.clear();
        VerticesNeeds.clear();
        EdgeL.clear();
    }

    //    void GLDraweMetaMesh()
    //    {
    //        MMesh.GLDraw();
    //    }

    void InitMetaMesh()
    {
        FixAmmittableValences();

        std::vector<std::vector<size_t> > CurrV,CurrDir;
        std::vector<bool> IsLoop;
        GetCurrVertDir(CurrV,CurrDir,IsLoop);

        std::vector<size_t> NarrowV,ConcaveV,ConvexV;
        GetVertexType(TVNarrow,NarrowV);
        GetVertexType(TVConcave,ConcaveV);
        GetVertexType(TVConvex,ConvexV);
        std::vector<size_t> FixedV=NarrowV;
        FixedV.insert(FixedV.end(),ConcaveV.begin(),ConcaveV.end());
        FixedV.insert(FixedV.end(),ConvexV.begin(),ConvexV.end());

        //std::vector<std::vector<ScalarType> >  SideLen;
        std::set<size_t > NonGenusPartitionIndexSet;
        GetTopologicallyNotOKPartitionsIndex(NonGenusPartitionIndexSet);
        MMesh.Init(&Mesh(),PartitionCorners,FacePartitions,CurrV,IsLoop,
                   FixedV,EdgeL,EDirTable,NonGenusPartitionIndexSet);
    }

    void RemoveMetaMeshStep()
    {
        size_t t0=clock();
        MMesh.MergeLoop(avgEdge*CClarkability,MinVal,MaxVal,CClarkability,avgEdge,match_valence);
        size_t t1=clock();
        std::vector<size_t> RemPaths;
        MMesh.GetRemainingPaths(RemPaths);
        std::set<size_t> RemPathsSet(RemPaths.begin(),RemPaths.end());
        for (size_t i=0;i<ChoosenPaths.size();i++)
        {
            if (RemPathsSet.count(i)>0)continue;
            ChoosenPaths[i].PathNodes.clear();
        }
        size_t t2=clock();
        //std::cout<<"0"<<std::endl;
        RemoveEmptyPaths();
        //std::cout<<"1"<<std::endl;
        MergeContiguousPaths(ChoosenPaths);
        //std::cout<<"2"<<std::endl;
        RemoveEmptyPaths();
        //std::cout<<"3"<<std::endl;
        size_t t3=clock();
        InitEdgeDirTable();
        size_t t4=clock();
        UpdatePartitionsFromChoosen(true);
        size_t t5=clock();

        Time_Collapse_Step0+=t1-t0;
        Time_Collapse_Step1+=t2-t1;
        Time_Collapse_Step2+=t3-t1;
        Time_Collapse_Step3+=t4-t3;
        Time_Collapse_Step4+=t5-t4;
    }

    PatchTracer(VertexFieldGraph<MeshType> &_VFGraph):VFGraph(_VFGraph)
    {
        split_on_removal=true;
        //avoid_increase_valence=true;
        //avoid_collapse_irregular=false;
        away_from_singular=true;
        match_valence=true;
        check_quality_functor=false;
        //        max_lenght_distortion=-1;//1.2;
        //        max_lenght_variance=-1;//2;
        CClarkability=1;
        sample_ratio=0.1;
        MinVal=3;
        MaxVal=5;
        AllReceivers=false;
        Concave_Need=1;
        AllowDarts=false;
        AllowSelfGluedPatch=false;
        CheckQuadrangulationLimits=true;
        //FirstBorder=false;
        AllowRemoveConcave=false;
        PrioMode=PrioModeLoop;
        CheckUVIntersection=false;
        ContinuosCheckUVInt=true;
        CheckTJunction=true;
        subInt=3;
        //max_patch_area=MeshArea(Mesh())*0.5;
        //TraceLoopsBorders=true;
    }
};


#endif
