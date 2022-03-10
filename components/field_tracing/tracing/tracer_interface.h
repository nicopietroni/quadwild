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

#ifndef PATCH_TRACER_INTERFACE
#define PATCH_TRACER_INTERFACE

#include "patch_tracer.h"

//class ComputationTimeStats
//{
//public:
//    //static size_t Time_InitSubP;
//    static size_t Time_InitSubPatches;
//    static size_t Time_TraceSubPatches;
//    static size_t Time_RetrieveSubPatches;
////    static size_t Time_AddingEmitters;
////    static size_t Time_AddingLoops;
////    static size_t Time_AddingBorder;
//};


// Basic subdivision class
template<class MESH_TYPE>
struct SplitLevTrace : public   std::unary_function<vcg::face::Pos<typename MESH_TYPE::FaceType> ,  typename MESH_TYPE::CoordType >
{
    typedef typename MESH_TYPE::VertexType VertexType;
    typedef typename MESH_TYPE::FaceType FaceType;
    typedef typename MESH_TYPE::CoordType CoordType;
    typedef typename MESH_TYPE::ScalarType ScalarType;

    typedef std::pair<CoordType,CoordType> EdgeCoordKey;

    std::map<EdgeCoordKey,CoordType> *SplitOps;

    void operator()(typename MESH_TYPE::VertexType &nv,vcg::face::Pos<typename MESH_TYPE::FaceType>  ep)
    {
        VertexType* v0=ep.f->V0(ep.z);
        VertexType* v1=ep.f->V1(ep.z);

        assert(v0!=v1);

        CoordType Pos0=v0->P();
        CoordType Pos1=v1->P();

        EdgeCoordKey CoordK(std::min(Pos0,Pos1),std::max(Pos0,Pos1));

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

    SplitLevTrace(std::map<EdgeCoordKey,CoordType> *_SplitOps){SplitOps=_SplitOps;}
    //SplitLevQ(){}
};

template <class MESH_TYPE>
class EdgePredTrace
{
    typedef typename MESH_TYPE::VertexType VertexType;
    typedef typename MESH_TYPE::FaceType FaceType;
    typedef typename MESH_TYPE::CoordType CoordType;
    typedef typename MESH_TYPE::ScalarType ScalarType;
    typedef std::pair<CoordType,CoordType> EdgeCoordKey;

    std::map<EdgeCoordKey,CoordType> *SplitOps;

public:

    bool operator()(vcg::face::Pos<typename MESH_TYPE::FaceType> ep) const
    {
        VertexType* v0=ep.f->V0(ep.z);
        VertexType* v1=ep.f->V1(ep.z);

        assert(v0!=v1);

        CoordType Pos0=v0->P();
        CoordType Pos1=v1->P();

        EdgeCoordKey CoordK(std::min(Pos0,Pos1),std::max(Pos0,Pos1));

        return (SplitOps->count(CoordK)>0);
    }

    EdgePredTrace(std::map<EdgeCoordKey,CoordType> *_SplitOps){SplitOps=_SplitOps;}
};



template <class MeshType>
void SplitAlongShap(MeshType &mesh)
{
    typedef typename MeshType::CoordType CoordType;
    typedef std::pair<CoordType,CoordType> EdgeCoordKey;

    //first mark the ones that must be splitted
    std::map<EdgeCoordKey,CoordType> SplitOps;

    mesh.SelectSharpFeatures();
    std::vector<std::pair<CoordType,CoordType> > SharpCoords;
    //for each face
    for (size_t i=0;i<mesh.face.size();i++)
    {
        //for each edge
        for (size_t j=0;j<3;j++)
        {
            CoordType Pos0=mesh.face[i].P0(j);
            CoordType Pos1=mesh.face[i].P1(j);
            EdgeCoordKey CoordK(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
            CoordType NewPos=(Pos0+Pos1)/2;
            if (mesh.face[i].IsFaceEdgeS(j))
            {
                SharpCoords.push_back(std::pair<CoordType,CoordType>(std::min(Pos0,NewPos),
                                                                     std::max(Pos0,NewPos)));
                SharpCoords.push_back(std::pair<CoordType,CoordType>(std::min(Pos1,NewPos),
                                                                     std::max(Pos1,NewPos)));
            }
            //see if it has been already marked as cut
            if (SplitOps.count(CoordK)>0)continue;

            SplitOps[CoordK]=NewPos;
        }
    }

    std::sort(SharpCoords.begin(),SharpCoords.end());
    auto last=std::unique(SharpCoords.begin(),SharpCoords.end());
    SharpCoords.erase(last, SharpCoords.end());

    //UpdateFromCoordPairs(const std::vector<std::pair<CoordType,CoordType> > &SharpCoords)

    //if (SplitOps.size()==0)return;

    std::cout<<"Performing "<<SplitOps.size()<< " split ops"<<std::endl;
    SplitLevTrace<MeshType> splMd(&SplitOps);
    EdgePredTrace<MeshType> eP(&SplitOps);

    //do the final split
    //bool done=
    vcg::tri::RefineE<MeshType,SplitLevTrace<MeshType>,EdgePredTrace<MeshType> >(mesh,splMd,eP);
    mesh.UpdateFromCoordPairs(SharpCoords);

    //return done;
}

//template <class MeshType>
//bool RemovePatches(const std::vector<size_t> &IndexPatches,
//                   PatchTracer<MeshType> &PTr)
//{
//    typedef typename MeshType::ScalarType ScalarType;

//    vcg::tri::io::ExporterPLY<MeshType>::Save(PTr.Mesh(),"test_before.ply");

//    //split to remove sub path portions
//    PTr.SplitIntoSubPaths();


//    //get the patches
//    std::vector<std::vector<size_t> > VertIdx;
//    std::vector<std::vector<size_t> > VertDir;
//    std::vector<bool> IsLoop;
//    PTr.GetCurrVertDir(VertIdx,VertDir,IsLoop);

//    //save original index of vert on quality
//    for (size_t i=0;i<PTr.Mesh().vert.size();i++)
//        PTr.Mesh().vert[i].Q()=i;

//    //select the one to remove
////    vcg::tri::UpdateSelection<MeshType>::FaceClear(PTr.Mesh());
////    vcg::tri::UpdateSelection<MeshType>::VertexClear(PTr.Mesh());
//    for (size_t i=0;i<IndexPatches.size();i++)
//        for (size_t j=0;j<PTr.Partitions[IndexPatches[i]].size();j++)
//        {
//            size_t IndexF=PTr.Partitions[IndexPatches[i]][j];
//            vcg::tri::Allocator<MeshType>::DeleteFace(PTr.Mesh(),PTr.Mesh().face[IndexF]);
//        }
//    vcg::tri::Clean<MeshType>::RemoveUnreferencedVertex(PTr.Mesh());
//    vcg::tri::Allocator<MeshType>::CompactEveryVector(PTr.Mesh());
////    vcg::tri::UpdateSelection<MeshType>::FaceClear(PTr.Mesh());
////    vcg::tri::UpdateSelection<MeshType>::VertexClear(PTr.Mesh());
//    PTr.Mesh().UpdateAttributes();

//    //get vertex map
//    std::map<size_t,size_t> VertMap;
//    for (size_t i=0;i<PTr.Mesh().vert.size();i++)
//    {
//        size_t OldIdx=PTr.Mesh().vert[i].Q();
//        VertMap[OldIdx]=i;
//    }

//    //remap vertices
//    for (size_t i=0;i<VertIdx.size();i++)
//        for (size_t j=0;j<VertIdx[i].size();j++)
//        {
//            size_t OldIdx=VertIdx[i][j];
//            //referring to a non existing vert
//            if (VertMap.count(OldIdx)==0)
//            {
//                VertIdx[i].clear();
//                break;
//            }
//            else
//            {
//                size_t NewIdx=VertMap[OldIdx];
//                VertIdx[i][j]=NewIdx;
//            }
//        }

//    //check new borders
//    std::set<std::pair<size_t,size_t> > BorderE;
//    //std::set<std::pair<size_t,size_t> > InternalE;
//    for (size_t i=0;i<PTr.Mesh().face.size();i++)
//        for (size_t j=0;j<3;j++)
//        {
//            size_t IndexV0=vcg::tri::Index(PTr.Mesh(),PTr.Mesh().face[i].V0(j));
//            size_t IndexV1=vcg::tri::Index(PTr.Mesh(),PTr.Mesh().face[i].V1(j));
//            std::pair<size_t,size_t>  key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
//            if (vcg::face::IsBorder(PTr.Mesh().face[i],j))
//                 BorderE.insert(key);
//                //                BorderE.insert(key);
//                //            else
//                //InternalE.insert(key);
//        }

//    //then filter the paths
//    std::vector<std::vector<size_t> > VertIdxSwap;
//    std::vector<std::vector<size_t> > VertDirSwap;
//    std::vector<bool> IsLoopSwap;
//    for (size_t i=0;i<VertIdx.size();i++)
//    {
//        bool to_delete=false;
//        if (VertIdx[i].size()==0)continue;

//        size_t Limit=VertIdx[i].size()-1;
//        if (IsLoop[i])Limit++;
//        size_t Size=VertIdx[i].size();
//        for (size_t j=0;j<Limit;j++)
//        {
//            size_t IndexV0=VertIdx[i][j];
//            size_t IndexV1=VertIdx[i][(j+1)%Size];
//            std::pair<size_t,size_t>  key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
//            if (BorderE.count(key)==0)continue;
//            to_delete=true;
//            std::cout<<"AAA"<<std::endl;
//            break;
//        }
//        if (!to_delete)
//        {
//            VertIdxSwap.push_back(VertIdx[i]);
//            VertDirSwap.push_back(VertDir[i]);
//            IsLoopSwap.push_back(IsLoop[i]);
//        }
//    }
//    VertIdx=VertIdxSwap;
//    VertDir=VertDirSwap;
//    IsLoop=IsLoopSwap;


//    std::set<std::pair<size_t,size_t> > AddedE;
//    for (size_t i=0;i<VertIdx.size();i++)
//    {
//        size_t Limit=VertIdx[i].size()-1;
//        if (IsLoop[i])Limit++;
//        size_t Size=VertIdx[i].size();
//        for (size_t j=0;j<Limit;j++)
//        {
//            size_t IndexV0=VertIdx[i][j];
//            size_t IndexV1=VertIdx[i][(j+1)%Size];
//            std::pair<size_t,size_t>  key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
//            assert(AddedE.count(key)==0);
//            assert(BorderE.count(key)==0);
//            AddedE.insert(key);
//            //std::cout<<"AAA"<<std::endl;
//            //break;
//        }
//    }


//    //update the graph
//    std::cout<<"A"<<std::endl;
//    ScalarType currDrift=PTr.Drift;
//    PTr.Reset();
//    PTr.InitTracer(currDrift,false);
//    //copy paths
//    std::cout<<"0"<<std::endl;
//    PTr.SetChoosenFromVertDir(VertIdx,VertDir,IsLoop);
//    std::cout<<"1"<<std::endl;
//    PTr.InitEdgeDirTable();
//    std::cout<<"2"<<std::endl;
//    PTr.InitEdgeL();
//    std::cout<<"3"<<std::endl;
//    vcg::tri::io::ExporterPLY<MeshType>::Save(PTr.Mesh(),"test_remain.ply");

//    PTr.UpdatePartitionsFromChoosen(true);
//    std::cout<<"4"<<std::endl;
//}

template <class TracerType>
bool TraceSubPatch(const size_t &IndexPatch,
                   TracerType &PTr,
                   std::vector<std::vector<size_t> > &VertIdx,
                   std::vector<std::vector<size_t> > &VertDir,
                   std::vector<bool> &IsLoop,
                   bool onlyneeded,
                   bool resample_loops,
                   bool DebugMsg,
                   bool force_always)
{
    typedef typename TracerType::MeshType MeshType;

    size_t t0=clock();
    
    //    //first copy the submesh
    //    for (size_t i=0;i<PTr.Mesh().vert.size();i++)
    //        PTr.Mesh().vert[i].Q()=i;

    MeshType SubMesh;
    PTr.GetPatchMesh(IndexPatch,SubMesh,false);
    SubMesh.UpdateAttributes();

    size_t t1=clock();

    std::vector<size_t> VertMap;
    for (size_t i=0;i<SubMesh.vert.size();i++)
        VertMap.push_back(SubMesh.vert[i].Q());

    //copy original normals
    for (size_t i=0;i<SubMesh.vert.size();i++)
        SubMesh.vert[i].N()=PTr.Mesh().vert[VertMap[i]].N();

    //make a subgraph
    VertexFieldGraph<MeshType> VFGraph(SubMesh);
    VFGraph.InitGraph();


    //initialize the tracer
    size_t t2=clock();
    //PatchTracer<MeshType> SubTr(VFGraph);
    TracerType SubTr(VFGraph);
    //    size_t t_2_0=clock();
    SubTr.CopyParametersFrom(PTr);
    //    size_t t_2_1=clock();
    SubTr.CopyFrom(PTr,VertMap,IndexPatch);
    //    size_t t_2_2=clock();
    SubTr.InitEdgeDirTable();
    size_t t3=clock();

    //    Time_InitSubPatches2_0+=t_2_0-t2;
    //    Time_InitSubPatches2_1+=t_2_1-t_2_0;
    //    Time_InitSubPatches2_2+=t_2_2-t_2_1;
    //    Time_InitSubPatches2_3+=t3-t_2_2;

    //SubTr.InitTraceableBorders();

    //SubTr.InitBorderDirMap();
    //this put 100 sample by default
    //SubTr.sample_ratio=-1;
    size_t Added_paths=SubTr.CopyPathsFrom(PTr,VertMap);
    SubTr.InitEdgeL();
    if (Added_paths>0)
    {
        if (DebugMsg)
            std::cout<<"ADDED "<<Added_paths<<" EXTRA PATH IN SUBDIVISION"<<std::endl;

        for (size_t i=0;i<SubTr.ChoosenPaths.size();i++)
            SubTr.ChoosenPaths[i].Unremovable=true;

        //added internal path, need to reupdate the table
        SubTr.InitEdgeDirTable();
        //SubTr.SampleLoopEmitters(false,MIN_SAMPLES_HARD);
    }
    if (resample_loops)
        SubTr.SampleLoopEmitters(false);
    
    size_t t4=clock();
    Time_InitSubPatches0+=t1-t0;
    Time_InitSubPatches1+=t2-t1;
    Time_InitSubPatches2+=t3-t2;
    Time_InitSubPatches3+=t4-t3;
    //SubTr.SampleLoopEmitters(false);

    //first update of the sub patch (fast and set needs)
    //SubTr.UpdatePartitionsFromChoosen(true);

    //then trace in the subpatch
    SubTr.DebugMsg=DebugMsg;
    SubTr.BatchAddLoops(true,onlyneeded,false,force_always);//inteleaveremoval,finalremoval);//,interleave_smooth);


    //copy back paths to the original
    SubTr.GetCurrVertDir(VertIdx,VertDir,IsLoop);

    //remove the first one in case they were already there
    if (Added_paths>0)
    {
        if (DebugMsg)
            std::cout<<"REMOVING EXTRA PATH"<<std::endl;

        VertIdx.erase(VertIdx.begin(), VertIdx.begin() + Added_paths);
        VertDir.erase(VertDir.begin(), VertDir.begin() + Added_paths);
        IsLoop.erase(IsLoop.begin(), IsLoop.begin() + Added_paths);
    }

    for (size_t i=0;i<VertIdx.size();i++)
        for (size_t j=0;j<VertIdx[i].size();j++)
            VertIdx[i][j]=VertMap[VertIdx[i][j]];

    size_t t5=clock();
    Time_TraceSubPatches+=t5-t4;

    return (VertIdx.size()>0);
}

void WriteUnsolvedStats(const std::vector<PatchType> &PatchTypes)
{
    size_t numLow=0;
    size_t numHigh=0;
    size_t numNonDisk=0;
    size_t numhasEmitter=0;
    size_t nummaxCCbility=0;
    size_t numwrongVal=0;
    size_t numOK=0;
    for (size_t i=0;i<PatchTypes.size();i++)
    {
        switch (PatchTypes[i])
        {
        case LowCorners: numLow++;break;
        case HighCorners: numHigh++;break;
        case NonDisk:numNonDisk++;break;
        case HasEmitter:numhasEmitter++;break;
        case MaxCClarkability:nummaxCCbility++;break;
        case NonMatchValence:numwrongVal++;break;
        default: numOK++;
        }
    }
    std::cout<<"** UNSATISFIED PATCHES **"<<std::endl;
    std::cout<<"*Num Patches:"<<numLow<<std::endl;
    std::cout<<"*Low Sides:"<<numLow<<std::endl;
    std::cout<<"*High Sides:"<<numHigh<<std::endl;
    std::cout<<"*Non Disks:"<<numNonDisk<<std::endl;
    std::cout<<"*Has Emit:"<<numhasEmitter<<std::endl;
    std::cout<<"*Max CCbility:"<<nummaxCCbility<<std::endl;
    std::cout<<"*Wrong Valence:"<<numwrongVal<<std::endl;
    std::cout<<"*IsOk:"<<numOK<<std::endl;
}

template <class TracerType>
bool FullTraced(TracerType &PTr,size_t &PartitionIndex)
{
    bool FullTraced=true;
    for (size_t i=0;i<PTr.Partitions[PartitionIndex].size();i++)
    {
        size_t IndexF=PTr.Partitions[PartitionIndex][i];
        FullTraced&=PTr.Mesh().face[IndexF].FullTraced;
        if (!FullTraced)return false;
    }
    return true;
}

template <class TracerType>
void FilterFullTracedPatches(TracerType &PTr,
                             std::vector<size_t> &PartitionIndexes,
                             bool DebugMsg=true)
{
    std::vector<size_t> PartitionIndexesSwap;
    for (size_t i=0;i<PartitionIndexes.size();i++)
    {
        if (FullTraced(PTr,PartitionIndexes[i]))continue;
        PartitionIndexesSwap.push_back(PartitionIndexes[i]);
    }
    if (DebugMsg)
        std::cout<<"Filtered out "<<PartitionIndexes.size()-PartitionIndexesSwap.size()<<" patches"<<std::endl;

    PartitionIndexes=PartitionIndexesSwap;
}

template <class TracerType>
void SolveSubPatches(TracerType &PTr,
                     bool onlyneeded,
                     bool only_non_disk,
                     bool force_always=false,
                     bool DebugMsg=true)
{
    std::vector<std::vector<size_t> > TotVertIdx;
    std::vector<std::vector<size_t> > TotVertDir;
    std::vector<bool> TotIsLoop;
    PTr.GetCurrVertDir(TotVertIdx,TotVertDir,TotIsLoop);

    bool solved=false;
    bool added_trace=false;
    std::vector<size_t> UnsolvedPartitionIndex;
    std::vector<PatchType> PatchTypes;
    bool augmented_max_narrow_dist=false;
    //std::vector<size_t> AddedPathOnStep;
    do{
        std::vector<std::vector<size_t> > OLDVertIdx=TotVertIdx;
        std::vector<std::vector<size_t> > OLDVertDir=TotVertDir;
        std::vector<bool> OLDIsLoop=TotIsLoop;

        PTr.LazyUpdatePartitions();
        if (DebugMsg)
            std::cout<<"**** RETRIEVING NON OK PATCHES ****"<<std::endl;

        if (!only_non_disk)
            PTr.GetUnsolvedPartitionsIndex(UnsolvedPartitionIndex,PatchTypes);
        else
            PTr.GetTopologicallyNotOKPartitionsIndex(UnsolvedPartitionIndex);

        if (UnsolvedPartitionIndex.size()==0)
            solved=true;

        if (!only_non_disk)
            FilterFullTracedPatches(PTr,UnsolvedPartitionIndex,DebugMsg);

        if (DebugMsg)
            std::cout<<"**** SUBPATCH TRACING - THERE ARE "<<UnsolvedPartitionIndex.size()<<" Unsolved Partitions ****"<<std::endl;

        if ((!only_non_disk)&&(DebugMsg))
            WriteUnsolvedStats(PatchTypes);
        //PTr.WriteInfo();

        //first copy the submesh
        for (size_t i=0;i<PTr.Mesh().vert.size();i++)
            PTr.Mesh().vert[i].Q()=i;

        size_t Has_traced=0;
        for(size_t i=0;i<UnsolvedPartitionIndex.size();i++)
        {

            std::vector<std::vector<size_t> > NewVertIdx;
            std::vector<std::vector<size_t> > NewVertDir;
            std::vector<bool> NewIsLoop;
            size_t currPartIndex=UnsolvedPartitionIndex[i];
            bool traced=TraceSubPatch<TracerType>(currPartIndex,PTr,NewVertIdx,NewVertDir,
                                                  NewIsLoop,onlyneeded,only_non_disk,false,force_always);
            if (!traced)
            {
                for (size_t j=0;j<PTr.Partitions[currPartIndex].size();j++)
                {
                    size_t IndxF=PTr.Partitions[currPartIndex][j];
                    assert(IndxF<PTr.Mesh().face.size());
                    PTr.Mesh().face[IndxF].FullTraced=true;
                }
            }
            else
                Has_traced++;

            TotVertIdx.insert(TotVertIdx.end(),NewVertIdx.begin(),NewVertIdx.end());
            TotVertDir.insert(TotVertDir.end(),NewVertDir.begin(),NewVertDir.end());
            TotIsLoop.insert(TotIsLoop.end(),NewIsLoop.begin(),NewIsLoop.end());

        }
        if (DebugMsg)
        {
            std::cout<<"Has traced into "<<Has_traced<<" patches"<<std::endl;
            std::cout<<"Updating Patches"<<std::endl;
        }
        PTr.SetChoosenFromVertDir(TotVertIdx,TotVertDir,TotIsLoop);
        PTr.InitEdgeDirTable();

        if (DebugMsg)
            std::cout<<"Done Updating Patches"<<std::endl;

        added_trace=false;
        added_trace|=(OLDVertIdx!=TotVertIdx);
        added_trace|=(OLDVertDir!=TotVertDir);
        added_trace|=(OLDIsLoop!=TotIsLoop);

        //add a last step withmore tolerance in case
        //it cannot close the concave/narrow
        if ((!only_non_disk)&&(!added_trace)&&
                (PTr.HasIncompleteEmitter())&&
                (!augmented_max_narrow_dist))
        {
            augmented_max_narrow_dist=true;
            PTr.MaxNarrowWeight*=100;
            added_trace=true;
        }
    }while(added_trace & (!solved));
}

//template <class MeshType>
//void RemoveNonTopologicallyOK(PatchTracer<MeshType> &PTr)
//{
//    std::vector<size_t> To_Remove;
//    for (size_t i=0;i<PTr.Partitions.size();i++)
//    {
//        if (PTr.PartitionType[i]!=NonDisk)continue;
//        To_Remove.push_back(i);
//    }
//    if (To_Remove.size()>0);
//    RemovePatches(To_Remove,PTr);
//}

//template <class MeshType>
//void ForceSolvePatch(const std::vector<size_t> &IndexPatches,
//                     PatchTracer<MeshType> &PTr)
//{
//    typedef typename MeshType::ScalarType ScalarType;

//    //save all existing edges
//    std::set<std::pair<size_t,size_t> > ExistingE;
//    std::vector<std::vector<size_t> > CurrV,CurrDir;
//    std::vector<bool> IsLoop;
//    PTr.GetCurrVertDir(CurrV,CurrDir,IsLoop);
//    for (size_t i=0;i<CurrV.size();i++)
//    {
//        size_t Limit=CurrV[i].size()-1;
//        if (IsLoop[i])Limit++;
//        size_t Size=CurrV[i].size();
//        for (size_t j=0;j<Limit;j++)
//        {
//            size_t IndexV0=CurrV[i][j];
//            size_t IndexV1=CurrV[i][(j+1)%Size];
//            std::pair<size_t,size_t>  key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
//            ExistingE.insert(key);
//        }
//    }
//    //    for (size_t i=0;i<PTr.Mesh().face.size();i++)
//    //        for (size_t j=0;j<3;j++)
//    //        {
//    //            size_t IndexV0=vcg::tri::Index(PTr.Mesh(),PTr.Mesh().face[i].V0(j));
//    //            size_t IndexV1=vcg::tri::Index(PTr.Mesh(),PTr.Mesh().face[i].V1(j));
//    //            std::pair<size_t,size_t>  key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
//    //            ExistingE.insert(key);
//    //        }

//    for (size_t i=0;i<IndexPatches.size();i++)
//        for (size_t j=0;j<PTr.Partitions[IndexPatches[i]].size();j++)
//        {
//            size_t IndexF=PTr.Partitions[IndexPatches[i]][j];
//            for (size_t k=0;k<3;k++)
//            {
//                size_t IndexV0=vcg::tri::Index(PTr.Mesh(),PTr.Mesh().face[IndexF].V0(k));
//                size_t IndexV1=vcg::tri::Index(PTr.Mesh(),PTr.Mesh().face[IndexF].V1(k));
//                std::pair<size_t,size_t>  key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
//                if (ExistingE.count(key)>0)continue;
//                size_t IndexN0,IndexN1;
//                PTr.GetEdgeNodes(IndexV0,IndexV1,IndexN0,IndexN1);
//                PTr.AddSingleEdgePath(IndexN0,IndexN1);
//            }
//        }

//    PTr.UpdatePartitionsFromChoosen(true);
//}


//template <class MeshType>
//void SolveNonTopologicallyOK(PatchTracer<MeshType> &PTr)
//{
//    std::vector<size_t> To_Solve;
//    for (size_t i=0;i<PTr.Partitions.size();i++)
//    {
//        if (PTr.PartitionType[i]!=NonDisk)continue;
//        To_Solve.push_back(i);
//    }
//    if (To_Solve.size()>0);
//    ForceSolvePatch(To_Solve,PTr);
//}

template <class TracerType>
void RecursiveProcess(TracerType &PTr,
                      const typename TracerType::ScalarType Drift,
                      bool onlyneeded,
                      bool finalremoval,
                      bool PreRemoveStep=true,
                      bool UseMetamesh=true,
                      bool ForceMultiSplit=false,
                      bool CheckSurfaceFolds=true,
                      bool SmoothBeforeRemove=false,
                      bool DebugMsg=true)
{
    typedef typename TracerType::ScalarType ScalarType;

    //cannot do metacollapse with the following conditions
    assert(!(PTr.AllowDarts && UseMetamesh));
    assert(!(PTr.AllowSelfGluedPatch && UseMetamesh));
    assert(!((!PTr.CheckQuadrangulationLimits) && UseMetamesh));
    assert(!((PTr.MinVal<2) && UseMetamesh));

    Time_FirstTrace=0;
    Time_InitSubPatches0=0;
    Time_InitSubPatches1=0;
    Time_InitSubPatches2=0;
    Time_InitSubPatches3=0;
    Time_TraceSubPatches=0;
    Time_UpdatePartitionsTotal=0;
    Time_UpdatePartitionsLazy=0;
    Time_Removal=0;

    int t0=clock();
    //do a first step of tracing
    //PTr.Init(Drift,true);



    if (DebugMsg)
        std::cout<<"**** FIRST TRACING STEP ****"<<std::endl;

    int NumE0=PTr.NumEmitterType(TVNarrow);
    int NumE1=PTr.NumEmitterType(TVConcave);

    //in this case need to trace loops
    if ((NumE0==0)||(NumE1==0))
    {
        //PTr.DebugMsg=true;
        PTr.BatchAddLoops(false,onlyneeded,false,false);
    }
    else
    {
        PTr.BatchAddLoops(false,onlyneeded,ForceMultiSplit,false);
    }
    //if no patch has been created than trace loops
    if ((ForceMultiSplit==true)&&(PTr.Partitions.size()==0))
        PTr.BatchAddLoops(false,onlyneeded,false,false);

    int t1=clock();
    Time_FirstTrace+=t1-t0;

    //then do a first partitions update
    if (DebugMsg)
        std::cout<<"Updating Patches"<<std::endl;

    //THIS SHOULD SPEED UP
    PTr.LazyUpdatePartitions();
    //PTr.UpdatePartitionsFromChoosen(true);

    if (DebugMsg)
    {
        std::cout<<"Updated"<<std::endl;
    }
    //solve sub patches normally
    if (DebugMsg)
        std::cout<<"**** FIRST SUBTRACING STEP ****"<<std::endl;

    SolveSubPatches(PTr,onlyneeded,false,false,DebugMsg);

    //then check if there is some non-disk-like patches
    if (DebugMsg)
        std::cout<<"**** CHECK NO DISK ONES ****"<<std::endl;

    SolveSubPatches(PTr,onlyneeded,true,false,DebugMsg);


    if (DebugMsg)
        std::cout<<"**** FORCE SOLVING NO DISK ONES ****"<<std::endl;


    SolveSubPatches(PTr,onlyneeded,true,true,DebugMsg);

    std::vector<size_t> UnsolvedPartitionIndex;
    std::vector<PatchType> PatchTypes;

    if (DebugMsg)
        std::cout<<"**** After All Insertion Steps ****"<<std::endl;

    int t2=clock();


    if (finalremoval)
    {
        if (SmoothBeforeRemove)
            PTr.SmoothPatches(3,0.5,CheckSurfaceFolds);

        if (DebugMsg)
            std::cout<<"**** FINAL REMOVAL ****"<<std::endl;

        //PTr.UpdatePartitionsFromChoosen(true);
        PTr.SetAllRemovable();

        bool HasNonManif=false;
        std::vector<size_t > NonGenusOK;
        //        PTr.UpdatePartitionsFromChoosen(true);
        PTr.GetTopologicallyNotOKPartitionsIndex(NonGenusOK);
        if (NonGenusOK.size()>0)
            HasNonManif=true;


        if (UseMetamesh)// && (!HasNonManif))
            PTr.BatchRemovalMetaMesh(PreRemoveStep);
        else
            PTr.BatchRemovalOnMesh(PreRemoveStep);

        if (DebugMsg)
            std::cout<<"**** After Last Removal Step ****"<<std::endl;

        //        PTr.CutEarPath();
        if (HasNonManif)
        {
            PTr.UpdatePartitionsFromChoosen(true);
            PTr.WriteInfo();
        }
    }
    else
    {
        //PTr.CutEarPath();
        PTr.UpdatePartitionsFromChoosen(true);
    }

    //PTr.CutEarPath();

    PTr.UpdatePartitionsFromChoosen(false);

    int t3=clock();
    Time_Removal+=t3-t2;

    if (DebugMsg)
        std::cout<<"Updating Patches"<<std::endl;
    PTr.GetUnsolvedPartitionsIndex(UnsolvedPartitionIndex,PatchTypes);

    if (DebugMsg)
    {
        std::cout<<"**** FINAL THERE ARE "<<UnsolvedPartitionIndex.size()<<" Unsolved Partitions ****"<<std::endl;
        std::cout<<"**** TOTAL "<<PTr.Partitions.size()<<" Partitions ****"<<std::endl;
        std::cout<<"Updated"<<std::endl;
    }

    if (DebugMsg)
        std::cout<<"Smoothing"<<std::endl;
    PTr.SmoothPatches(10,0.5,CheckSurfaceFolds);
    if (DebugMsg)
        std::cout<<"End Smoothing"<<std::endl;

    //    std::cout<<"Fix Valences"<<std::endl;
    PTr.FixValences();

    if (DebugMsg)
        PTr.WriteInfo();

    int t4=clock();
    ScalarType ElpsedSec=(ScalarType)(t4-t0)/CLOCKS_PER_SEC;

    if (DebugMsg)
        std::cout<<"**** FINAL ELAPSED TIME "<<ElpsedSec<<" Seconds ****"<<std::endl;

    //    //check if metamesh has not been initialized already
    //    if (!(finalremoval && UseMetamesh))
    //        PTr.InitMetaMesh();

    //    std::cout<<"Time First Trace "<<(ScalarType)Time_FirstTrace/CLOCKS_PER_SEC<<std::endl;
    //    std::cout<<"Time Init SubPatches 0 "<<(ScalarType)Time_InitSubPatches0/CLOCKS_PER_SEC<<std::endl;
    ////    std::cout<<"Time Init SubPatches 0 - 0 "<<(ScalarType)Time_InitSubPatches0_0/CLOCKS_PER_SEC<<std::endl;
    ////    std::cout<<"Time Init SubPatches 0 - 1 "<<(ScalarType)Time_InitSubPatches0_1/CLOCKS_PER_SEC<<std::endl;

    //    std::cout<<"Time Init SubPatches 1 "<<(ScalarType)Time_InitSubPatches1/CLOCKS_PER_SEC<<std::endl;
    //    std::cout<<"Time Init SubPatches 2 "<<(ScalarType)Time_InitSubPatches2/CLOCKS_PER_SEC<<std::endl;

    ////    std::cout<<"---time Init SubPatches 2 0 "<<(ScalarType)Time_InitSubPatches2_0/CLOCKS_PER_SEC<<std::endl;
    ////    std::cout<<"---time Init SubPatches 2 1 "<<(ScalarType)Time_InitSubPatches2_1/CLOCKS_PER_SEC<<std::endl;
    ////    std::cout<<"---time Init SubPatches 2 2 "<<(ScalarType)Time_InitSubPatches2_2/CLOCKS_PER_SEC<<std::endl;
    ////    std::cout<<"---time Init SubPatches 2 3 "<<(ScalarType)Time_InitSubPatches2_3/CLOCKS_PER_SEC<<std::endl;
    ////    std::cout<<"---time Init SubPatches 2 4 "<<(ScalarType)Time_InitSubPatches2_4/CLOCKS_PER_SEC<<std::endl;
    ////    std::cout<<"---time Init SubPatches 2 5 "<<(ScalarType)Time_InitSubPatches2_5/CLOCKS_PER_SEC<<std::endl;
    ////    std::cout<<"---time Init SubPatches 2 6 "<<(ScalarType)Time_InitSubPatches2_6/CLOCKS_PER_SEC<<std::endl;
    ////    std::cout<<"---time Init SubPatches 2 7 "<<(ScalarType)Time_InitSubPatches2_7/CLOCKS_PER_SEC<<std::endl;

    //    std::cout<<"Time Init SubPatches 3 "<<(ScalarType)Time_InitSubPatches3/CLOCKS_PER_SEC<<std::endl;
    //    std::cout<<"Time Trace SubPatches "<<(ScalarType)Time_TraceSubPatches/CLOCKS_PER_SEC<<std::endl;

    //    std::cout<<"Time Update Partition Total "<<(ScalarType)Time_UpdatePartitionsTotal/CLOCKS_PER_SEC<<std::endl;
    //    std::cout<<"Time Update Partition Lazy "<<(ScalarType)Time_UpdatePartitionsLazy/CLOCKS_PER_SEC<<std::endl;
    //    std::cout<<"Time Removal "<<(ScalarType)Time_Removal/CLOCKS_PER_SEC<<std::endl;
    //    std::cout<<"Time Init MMesh "<<(ScalarType)Time_InitMetaMesh/CLOCKS_PER_SEC<<std::endl;
    //    std::cout<<"Time Collapse MMesh 0 "<<(ScalarType)Time_Collapse_Step0/CLOCKS_PER_SEC<<std::endl;
    //    std::cout<<"Time Collapse MMesh 1 "<<(ScalarType)Time_Collapse_Step1/CLOCKS_PER_SEC<<std::endl;
    //    std::cout<<"Time Collapse MMesh 2 "<<(ScalarType)Time_Collapse_Step2/CLOCKS_PER_SEC<<std::endl;
    //    std::cout<<"Time Collapse MMesh 3 "<<(ScalarType)Time_Collapse_Step3/CLOCKS_PER_SEC<<std::endl;
    //    std::cout<<"Time Collapse MMesh 4 "<<(ScalarType)Time_Collapse_Step4/CLOCKS_PER_SEC<<std::endl;
    //    //PTr.InitMetaMesh();
}

template <class TracerType>
void RecursiveProcessWithDarts(TracerType &PTr,
                               const typename TracerType::ScalarType Drift,
                               bool onlyneeded,
                               bool finalremoval,
                               bool PreRemoveStep,
                               bool UseMetamesh,
                               bool ForceMultiSplit,
                               bool CheckSurfaceFolds,
                               const std::vector<typename TracerType::ScalarType> &DartPriority,
                               bool SmoothBeforeRemove,
                               bool DebugMsg)
{
    PTr.AllowDarts=true;
    RecursiveProcess(PTr,Drift,onlyneeded,finalremoval,PreRemoveStep,
                     UseMetamesh,ForceMultiSplit,CheckSurfaceFolds,SmoothBeforeRemove,DebugMsg);

    //    std::cout<<"DART STEP"<<std::endl;
    if (finalremoval)
    {
        //in this case either remove all or none
        if (!PTr.check_quality_functor)
        {
            PTr.SetAllRemovable();
            PTr.AllowDarts=true;
            PTr.MinVal=0;
            PTr.split_on_removal=true;
            PTr.CheckQuadrangulationLimits=false;
            //split at intersectons
            PTr.SetPriorityVect(DartPriority);
            PTr.BatchRemovalOnMesh(true);
            PTr.MergeContiguousPaths();
        }
        //in this case either remove gradually
        else
        {
            PTr.SetAllRemovable();
            PTr.AllowDarts=true;
            PTr.MinVal=0;
            PTr.split_on_removal=true;
            PTr.CheckQuadrangulationLimits=false;
            PTr.SetPriorityVect(DartPriority);
            PTr.RemoveDarts();
            PTr.ClearPriorityVect();
            PTr.MergeContiguousPaths();
        }
    }
}


template <class TracerType>
void RecursiveProcessForTexturing(TracerType &PTr,
                                  const typename TracerType::ScalarType Drift,
                                  bool onlyneeded,
                                  bool finalremoval,
                                  bool PreRemoveStep,
                                  bool UseMetamesh,
                                  bool ForceMultiSplit,
                                  bool CheckSurfaceFolds,
                                  bool SmoothBeforeRemove,
                                  bool DebugMsg)
{
    PTr.AllowDarts=false;
    PTr.AllowSelfGluedPatch=true;
    //    PTr.MinVal=0;
    //    PTr.CheckQuadrangulationLimits=false;

    RecursiveProcess(PTr,Drift,onlyneeded,finalremoval,PreRemoveStep,
                     UseMetamesh,ForceMultiSplit,CheckSurfaceFolds,SmoothBeforeRemove,DebugMsg);

    //then make a remove step
    //PTr.SplitIntoSubPaths();
    if (finalremoval)
    {
        PTr.SetAllRemovable();
        PTr.AllowDarts=false;
        PTr.AllowSelfGluedPatch=true;
        PTr.MinVal=0;
        PTr.split_on_removal=true;
        //PTr.match_valence=false;
        PTr.CheckQuadrangulationLimits=false;
        PTr.BatchRemovalOnMesh(true);
        PTr.MergeContiguousPaths();
    }
}

template <class TracerType>
void RecursiveProcessForTexturingWithDarts(TracerType &PTr,
                                           const typename TracerType::ScalarType Drift,
                                           bool onlyneeded,
                                           bool finalremoval,
                                           bool PreRemoveStep,
                                           bool UseMetamesh,
                                           bool ForceMultiSplit,
                                           bool CheckSurfaceFolds,
                                           const std::vector<typename TracerType::ScalarType> &DartPriority,
                                           bool SmoothBeforeRemove,
                                           bool DebugMsg)
{

    //remove darts later
    PTr.AllowDarts=true;
    //PTr.AllowDarts=false;
    PTr.AllowSelfGluedPatch=true;
    PTr.CheckQuadrangulationLimits=false;

    RecursiveProcess(PTr,Drift,onlyneeded,finalremoval,PreRemoveStep,
                     UseMetamesh,ForceMultiSplit,CheckSurfaceFolds,SmoothBeforeRemove,DebugMsg);

    if (finalremoval)
    {
        //size_t OldVal=PTr.MaxVal;
        //PTr.MaxVal=PTr.MaxVal*2;
        //PTr.MaxVal=PTr.MaxVal*2;
        //in this case either remove all or none
        if (!PTr.check_quality_functor)
        {
            PTr.SetAllRemovable();
            PTr.AllowDarts=true;
            PTr.MinVal=0;
            PTr.split_on_removal=true;

            //split at intersectons
            PTr.SetPriorityVect(DartPriority);
            PTr.BatchRemovalOnMesh(true);
            PTr.MergeContiguousPaths();
        }
        //in this case either remove gradually
        else
        {
            PTr.SetAllRemovable();
            PTr.AllowDarts=true;
            PTr.MinVal=0;
            PTr.split_on_removal=true;

            PTr.SetPriorityVect(DartPriority);
            PTr.RemoveDarts();
            PTr.ClearPriorityVect();
            PTr.MergeContiguousPaths();
        }
        //PTr.MaxVal=OldVal;
    }
}

template <class MeshType>
void SaveCSV(PatchTracer<MeshType> &PTr,
             const std::string &pathProject,
             const size_t CurrNum)
{
    typename PatchTracer<MeshType>::PatchInfoType PInfo;
    PTr.GetInfo(PInfo);

    std::string pathCSVFinal=pathProject;
    pathCSVFinal=pathCSVFinal+"_p"+std::to_string(CurrNum)+".csv";
    FILE *f=fopen(pathCSVFinal.c_str(),"wt");
    assert(f!=NULL);
    fprintf(f,"%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d \n",pathProject.c_str(),
            (int)PInfo.NumPatchs,
            (int)PInfo.HasEmit,
            (int)PInfo.HighC,
            (int)PInfo.LowC,
            (int)PInfo.NonDiskLike,
            (int)PInfo.SizePatches[0],
            (int)PInfo.SizePatches[1],
            (int)PInfo.SizePatches[2],
            (int)PInfo.SizePatches[3],
            (int)PInfo.SizePatches[4],
            (int)PInfo.SizePatches[5],
            (int)PInfo.SizePatches[6],
            (int)PInfo.SizePatches[7]);
    fclose(f);
}


//template <class MeshType>
//void SaveAllData(PatchTracer<MeshType> &PTr,
//                 const std::string &pathProject,
//                 const size_t CurrNum,
//                 bool subdivide_when_save,
//                 bool save_origFace)
//{
//    typedef typename MeshType::CoordType CoordType;

//    std::vector<std::pair<CoordType,CoordType> > SharpCoords;
//    PTr.Mesh().GetSharpCoordPairs(SharpCoords);

//    //copy the mesh
//    MeshType SaveM;
//    vcg::tri::Append<MeshType,MeshType>::Mesh(SaveM,PTr.Mesh());
//    std::vector<size_t> SharpCorners;
//    PTr.getCornerSharp(SharpCorners);
//    std::set<CoordType> SharpCornerPos;
//    for(size_t i=0;i<SharpCorners.size();i++)
//        SharpCornerPos.insert(PTr.Mesh().vert[SharpCorners[i]].P());

//    //merge vertices
//    vcg::tri::Clean<MeshType>::RemoveDuplicateVertex(SaveM);
//    vcg::tri::Clean<MeshType>::RemoveUnreferencedVertex(SaveM);
//    vcg::tri::Allocator<MeshType>::CompactEveryVector(SaveM);
//    SaveM.UpdateAttributes();
//    SaveM.UpdateFromCoordPairs(SharpCoords);

//    //save vert pos
//    std::map<CoordType,size_t> VertMap;
//    for (size_t i=0;i<SaveM.vert.size();i++)
//        VertMap[SaveM.vert[i].P()]=i;

//    for (size_t i=0;i<SaveM.face.size();i++)
//        SaveM.face[i].Q()=PTr.FacePartitions[i];

//    //vcg::tri::io::ExporterPLY<MeshType>::Save(SaveM,"test_final.ply");


//    if (subdivide_when_save)
//    {
//        SplitAlongShap(SaveM);
//    }

//    //update sharp vertices
//    SaveM.SharpCorners.clear();
//    for (size_t i=0;i<SaveM.vert.size();i++)
//        if (SharpCornerPos.count(SaveM.vert[i].P())>0)
//            SaveM.SharpCorners.push_back(i);

//    //save the mesh
//    int Mask=0;
//    std::string pathMeshFinal=pathProject;
//    pathMeshFinal=pathMeshFinal+"_p"+std::to_string(CurrNum)+".obj";
//    vcg::tri::io::ExporterOBJ<MeshType>::Save(SaveM,pathMeshFinal.c_str(),Mask);

//    std::string pathPartitions=pathProject;
//    //pathPartitions.append("_p.patch");
//    pathPartitions=pathPartitions+"_p"+std::to_string(CurrNum)+".patch";
//    FILE *F=fopen(pathPartitions.c_str(),"wt");
//    assert(F!=NULL);
//    //    assert(PTr.FacePartitions.size()==PTr.Mesh().face.size());
//    //    fprintf(F,"%d\n",PTr.FacePartitions.size());
//    //    for (size_t i=0;i<PTr.FacePartitions.size();i++)
//    //        fprintf(F,"%d\n",PTr.FacePartitions[i]);
//    //    fclose(F);
//    fprintf(F,"%d\n",SaveM.face.size());
//    for (size_t i=0;i<SaveM.face.size();i++)
//        fprintf(F,"%d\n",(int)SaveM.face[i].Q());
//    fclose(F);

//    std::string pathCorners=pathProject;
//    //pathCorners.append("_p.corners");
//    pathCorners=pathCorners+"_p"+std::to_string(CurrNum)+".corners";
//    F=fopen(pathCorners.c_str(),"wt");
//    assert(F!=NULL);
//    assert(PTr.PartitionCorners.size()==PTr.Partitions.size());
//    fprintf(F,"%d\n",PTr.PartitionCorners.size());

//    //std::cout<<"**** SAVING MESH ****"<<std::endl;
//    PTr.WriteInfo();
//    for (size_t i=0;i<PTr.PartitionCorners.size();i++)
//    {
//        //check uniqueness
//        std::set<size_t> TestCorn;
//        fprintf(F,"%d\n",PTr.PartitionCorners[i].size());

//        //        std::cout<<"Writing Patch "<<i<<std::endl;
//        //        std::cout<<"Size "<<PTr.PartitionCorners[i].size()<<std::endl;
//        for (size_t j=0;j<PTr.PartitionCorners[i].size();j++)
//        {
//            //            if (PTr.PartitionCorners[i].size()<MIN_ADMITTIBLE)
//            //            {
//            //                PTr.GetPatchMesh(i,mesh);
//            //                vcg::tri::io::ExporterPLY<MeshType>::Save(mesh,"")
//            //                assert(0);
//            //            }
//            assert(PTr.PartitionCorners[i].size()>=MIN_ADMITTIBLE);
//            assert(PTr.PartitionCorners[i].size()<=MAX_ADMITTIBLE);
//            size_t IndexV=PTr.PartitionCorners[i][j];
//            CoordType CornerPos=PTr.Mesh().vert[IndexV].P();
//            assert(VertMap.count(CornerPos)>0);
//            fprintf(F,"%d\n",VertMap[CornerPos]);
//            TestCorn.insert(VertMap[CornerPos]);
//            //std::cout<<VertMap[CornerPos]<<",";
//        }
//        assert(TestCorn.size()==PTr.PartitionCorners[i].size());

//        //        std::cout<<std::endl;
//        //        for (size_t j=0;j<PTr.PartitionCorners[i].size();j++)
//        //        {
//        //            size_t IndexV=PTr.PartitionCorners[i][j];
//        //            std::cout<<IndexV<<",";
//        //        }
//        //        std::cout<<std::endl;
//        //        std::cout<<std::endl;
//    }
//    fclose(F);

//    std::string featurePartitions=pathProject;
//    featurePartitions=featurePartitions+"_p"+std::to_string(CurrNum)+".feature";
//    SaveM.SaveFeatures(featurePartitions);

//    std::string featureCorners=pathProject;
//    featureCorners=featureCorners+"_p"+std::to_string(CurrNum)+".c_feature";
//    SaveM.SaveSharpCorners(featureCorners);

//    if (save_origFace)
//    {
//        std::string origfaceP=pathProject;
//        origfaceP=featureCorners+"_p"+std::to_string(CurrNum)+"_origf.txt";
//        SaveM.SaveOrigFace(origfaceP);
//    }
//}


template <class MeshType>
void SaveAllData(PatchTracer<MeshType> &PTr,
                 const std::string &pathProject,
                 const size_t CurrNum,
                 bool subdivide_when_save,
                 bool save_origFace)
{
    typedef typename MeshType::CoordType CoordType;

    std::vector<std::pair<CoordType,CoordType> > SharpCoords;
    PTr.Mesh().GetSharpCoordPairs(SharpCoords);

    //copy the mesh
    MeshType SaveM;
    vcg::tri::Append<MeshType,MeshType>::Mesh(SaveM,PTr.Mesh());
    std::vector<size_t> SharpCorners;
    PTr.getCornerSharp(SharpCorners);
    std::set<CoordType> SharpCornerPos;
    for(size_t i=0;i<SharpCorners.size();i++)
        SharpCornerPos.insert(PTr.Mesh().vert[SharpCorners[i]].P());

    for (size_t i=0;i<SaveM.face.size();i++)
        SaveM.face[i].Q()=PTr.FacePartitions[i];

    //merge across sharp features
    vcg::tri::Clean<MeshType>::RemoveDuplicateVertex(SaveM);
    vcg::tri::Clean<MeshType>::RemoveUnreferencedVertex(SaveM);
    vcg::tri::Allocator<MeshType>::CompactEveryVector(SaveM);
    SaveM.UpdateAttributes();

    //then remove the non disk-like
    std::set<size_t > ToErasePartitions;
    PTr.GetTopologicallyNotOKPartitionsIndex(ToErasePartitions);
    std::vector<std::vector<size_t> > PatchCorners;

    //add the patches of faces that has non-manifoldness (if occours)
    int num=vcg::tri::Clean<MeshType>::CountNonManifoldEdgeFF(SaveM,true);
    if (num>0)
    {
        for (size_t i=0;i<SaveM.face.size();i++)
        {
            if (!SaveM.face[i].IsS())continue;
            size_t CurrP=SaveM.face[i].Q();
            ToErasePartitions.insert(CurrP);
        }
    }
    bool has_erased=false;
    if (ToErasePartitions.size()>0)
    {
        std::cout<<ToErasePartitions.size()<<"- NON DISK OR NON MANIF PATCH ERASED- "<<std::endl;
        has_erased=true;
        //first remove the faces
        size_t NumPart=0;
        for (size_t i=0;i<SaveM.face.size();i++)
        {
            size_t IndexP=SaveM.face[i].Q();
            NumPart=std::max(NumPart,IndexP);
            if (ToErasePartitions.count(IndexP)==0)continue;
            vcg::tri::Allocator<MeshType>::DeleteFace(SaveM,SaveM.face[i]);
        }

        //then remap the partitions
        std::vector<int> PartitionMap(NumPart+1,-1);
        size_t currP=0;
        for (size_t i=0;i<PartitionMap.size();i++)
        {
            if (ToErasePartitions.count(i)>0)
                continue;//in this case do not consider that patch

            //otherwise add the patch corners and the index too
            PatchCorners.push_back(PTr.PartitionCorners[i]);
            PartitionMap[i]=currP;
            currP++;
        }

        //final check
        assert((currP+ToErasePartitions.size()-1)==NumPart);

        //remap faces partitions
        for (size_t i=0;i<SaveM.face.size();i++)
        {
            if (SaveM.face[i].IsD())continue;
            int currP=SaveM.face[i].Q();
            assert(currP>=0);
            assert(currP<PartitionMap.size());
            int newPIdx=PartitionMap[currP];
            assert(newPIdx>=0);
            SaveM.face[i].Q()=newPIdx;
        }
    }
    else
        PatchCorners=PTr.PartitionCorners;

    //merge vertices
    vcg::tri::Clean<MeshType>::RemoveDuplicateVertex(SaveM);
    vcg::tri::Clean<MeshType>::RemoveUnreferencedVertex(SaveM);
    vcg::tri::Allocator<MeshType>::CompactEveryVector(SaveM);
    SaveM.UpdateAttributes();
    SaveM.UpdateFromCoordPairs(SharpCoords,false);

    //save vert pos
    std::map<CoordType,size_t> VertMap;
    for (size_t i=0;i<SaveM.vert.size();i++)
        VertMap[SaveM.vert[i].P()]=i;

    //vcg::tri::io::ExporterPLY<MeshType>::Save(SaveM,"test_final.ply");


    if (subdivide_when_save)
    {
        SplitAlongShap(SaveM);
    }

    //update sharp vertices
    SaveM.SharpCorners.clear();
    for (size_t i=0;i<SaveM.vert.size();i++)
        if (SharpCornerPos.count(SaveM.vert[i].P())>0)
            SaveM.SharpCorners.push_back(i);

    //save the mesh
    int Mask=0;
    std::string pathMeshFinal=pathProject;
    pathMeshFinal=pathMeshFinal+"_p"+std::to_string(CurrNum)+".obj";
    vcg::tri::io::ExporterOBJ<MeshType>::Save(SaveM,pathMeshFinal.c_str(),Mask);

    std::string pathPartitions=pathProject;
    //pathPartitions.append("_p.patch");
    pathPartitions=pathPartitions+"_p"+std::to_string(CurrNum)+".patch";
    FILE *F=fopen(pathPartitions.c_str(),"wt");
    assert(F!=NULL);
    //    assert(PTr.FacePartitions.size()==PTr.Mesh().face.size());
    //    fprintf(F,"%d\n",PTr.FacePartitions.size());
    //    for (size_t i=0;i<PTr.FacePartitions.size();i++)
    //        fprintf(F,"%d\n",PTr.FacePartitions[i]);
    //    fclose(F);
    fprintf(F,"%d\n",SaveM.face.size());
    for (size_t i=0;i<SaveM.face.size();i++)
        fprintf(F,"%d\n",(int)SaveM.face[i].Q());
    fclose(F);

    std::string pathCorners=pathProject;
    //pathCorners.append("_p.corners");
    pathCorners=pathCorners+"_p"+std::to_string(CurrNum)+".corners";
    F=fopen(pathCorners.c_str(),"wt");
    assert(F!=NULL);
    //assert(PTr.PartitionCorners.size()==PTr.Partitions.size());
    fprintf(F,"%d\n",PatchCorners.size());


    PTr.WriteInfo();

    std::cout<<"**** SAVING MESH ****"<<std::endl;
    for (size_t i=0;i<PatchCorners.size();i++)
    {
        //check uniqueness
        std::set<size_t> TestCorn;
        fprintf(F,"%d\n",PatchCorners[i].size());
        for (size_t j=0;j<PatchCorners[i].size();j++)
        {
            if (PTr.CheckQuadrangulationLimits)
            {
                assert(PatchCorners[i].size()>=MIN_ADMITTIBLE);
                assert(PatchCorners[i].size()<=MAX_ADMITTIBLE);
            }
            int IndexV=PatchCorners[i][j];//PTr.PartitionCorners[i][j];
            assert(IndexV<PTr.Mesh().vert.size());
            assert(IndexV>=0);
            CoordType CornerPos=PTr.Mesh().vert[IndexV].P();

            assert(VertMap.count(CornerPos)>0);
            fprintf(F,"%d\n",VertMap[CornerPos]);
            TestCorn.insert(VertMap[CornerPos]);
            //std::cout<<VertMap[CornerPos]<<",";
        }
        if (TestCorn.size()!=PatchCorners[i].size())
        {
            std::cout<<"WARNING DOUBLE VERT: "<<TestCorn.size()<<"!="<<PatchCorners[i].size()<<std::endl;
            //assert(0);
        }
        //assert(TestCorn.size()==PatchCorners[i].size());

        //        std::cout<<std::endl;
        //        for (size_t j=0;j<PTr.PartitionCorners[i].size();j++)
        //        {
        //            size_t IndexV=PTr.PartitionCorners[i][j];
        //            std::cout<<IndexV<<",";
        //        }
        //        std::cout<<std::endl;
        //        std::cout<<std::endl;
    }
    fclose(F);

    std::string featurePartitions=pathProject;
    featurePartitions=featurePartitions+"_p"+std::to_string(CurrNum)+".feature";
    SaveM.SaveFeatures(featurePartitions);

    std::string featureCorners=pathProject;
    featureCorners=featureCorners+"_p"+std::to_string(CurrNum)+".c_feature";
    SaveM.SaveSharpCorners(featureCorners);

    if (save_origFace)
    {
        std::string origfaceP=pathProject;
        origfaceP=featureCorners+"_p"+std::to_string(CurrNum)+"_origf.txt";
        SaveM.SaveOrigFace(origfaceP);
    }
}
#endif
