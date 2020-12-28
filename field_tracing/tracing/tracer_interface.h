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
struct SplitLev : public   std::unary_function<vcg::face::Pos<typename MESH_TYPE::FaceType> ,  typename MESH_TYPE::CoordType >
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

    SplitLev(std::map<EdgeCoordKey,CoordType> *_SplitOps){SplitOps=_SplitOps;}
    //SplitLevQ(){}
};

template <class MESH_TYPE>
class EdgePred
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

    EdgePred(std::map<EdgeCoordKey,CoordType> *_SplitOps){SplitOps=_SplitOps;}
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
    SplitLev<MeshType> splMd(&SplitOps);
    EdgePred<MeshType> eP(&SplitOps);

    //do the final split
    bool done=vcg::tri::RefineE<MeshType,SplitLev<MeshType>,EdgePred<MeshType> >(mesh,splMd,eP);
    mesh.UpdateFromCoordPairs(SharpCoords);

    //return done;
}

template <class MeshType>
bool TraceSubPatch(const size_t &IndexPatch,
                   PatchTracer<MeshType> &PTr,
                   std::vector<std::vector<size_t> > &VertIdx,
                   std::vector<std::vector<size_t> > &VertDir,
                   std::vector<bool> &IsLoop,
                   bool onlyneeded,
                   bool resample_loops,
                   bool DebugMsg)
{
    size_t t0=clock();
    
//    //first copy the submesh
//    for (size_t i=0;i<PTr.Mesh().vert.size();i++)
//        PTr.Mesh().vert[i].Q()=i;

    MeshType SubMesh;
    PTr.GetPatchMesh(IndexPatch,SubMesh);
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
    VFGraph.Init();


    //initialize the tracer
    size_t t2=clock();
    PatchTracer<MeshType> SubTr(VFGraph);
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
    SubTr.BatchAddLoops(true,onlyneeded);//inteleaveremoval,finalremoval);//,interleave_smooth);


    //copy back paths to the original
    SubTr.GetCurrVertDir(VertIdx,VertDir,IsLoop);

    //remove the first one in case they were already there
    if (Added_paths>0)
    {
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
    std::cout<<"*Low Sides:"<<numLow<<std::endl;
    std::cout<<"*High Sides:"<<numHigh<<std::endl;
    std::cout<<"*Non Disks:"<<numNonDisk<<std::endl;
    std::cout<<"*Has Emit:"<<numhasEmitter<<std::endl;
    std::cout<<"*Max CCbility:"<<nummaxCCbility<<std::endl;
    std::cout<<"*Wrong Valence:"<<numwrongVal<<std::endl;
    std::cout<<"*IsOk:"<<numOK<<std::endl;
}

template <class MeshType>
bool FullTraced(PatchTracer<MeshType> &PTr,size_t &PartitionIndex)
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

template <class MeshType>
void FilterFullTracedPatches(PatchTracer<MeshType> &PTr,std::vector<size_t> &PartitionIndexes)
{
    std::vector<size_t> PartitionIndexesSwap;
    for (size_t i=0;i<PartitionIndexes.size();i++)
    {
     if (FullTraced(PTr,PartitionIndexes[i]))continue;
     PartitionIndexesSwap.push_back(PartitionIndexes[i]);
    }
     std::cout<<"Filtered out "<<PartitionIndexes.size()-PartitionIndexesSwap.size()<<" patches"<<std::endl;
    PartitionIndexes=PartitionIndexesSwap;  
}

template <class MeshType>
void SolveSubPatches(PatchTracer<MeshType> &PTr,
                     bool onlyneeded,
                     bool only_non_disk)
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
        std::cout<<"**** RETRIEVING NON OK PATCHES ****"<<std::endl;
        if (!only_non_disk)
            PTr.GetUnsolvedPartitionsIndex(UnsolvedPartitionIndex,PatchTypes);
        else
            PTr.GetTopologicallyOKPartitionsIndex(UnsolvedPartitionIndex);

        if (UnsolvedPartitionIndex.size()==0)
            solved=true;

        FilterFullTracedPatches(PTr,UnsolvedPartitionIndex);

        std::cout<<"**** SUBPATCH TRACING - THERE ARE "<<UnsolvedPartitionIndex.size()<<" Unsolved Partitions ****"<<std::endl;
        WriteUnsolvedStats(PatchTypes);
        //PTr.WriteInfo();

        //first copy the submesh
        for (size_t i=0;i<PTr.Mesh().vert.size();i++)
            PTr.Mesh().vert[i].Q()=i;

        for(size_t i=0;i<UnsolvedPartitionIndex.size();i++)
        {

            std::vector<std::vector<size_t> > NewVertIdx;
            std::vector<std::vector<size_t> > NewVertDir;
            std::vector<bool> NewIsLoop;
            size_t currPartIndex=UnsolvedPartitionIndex[i];
            bool traced=TraceSubPatch<MeshType>(currPartIndex,PTr,
                                                NewVertIdx,NewVertDir,
                                                NewIsLoop,onlyneeded,only_non_disk,false);
            if (!traced)
            {
                for (size_t j=0;j<PTr.Partitions[currPartIndex].size();j++)
                {
                    size_t IndxF=PTr.Partitions[currPartIndex][j];
                    assert(IndxF<PTr.Mesh().face.size());
                    PTr.Mesh().face[j].FullTraced=true;
                }
            }
            TotVertIdx.insert(TotVertIdx.end(),NewVertIdx.begin(),NewVertIdx.end());
            TotVertDir.insert(TotVertDir.end(),NewVertDir.begin(),NewVertDir.end());
            TotIsLoop.insert(TotIsLoop.end(),NewIsLoop.begin(),NewIsLoop.end());

        }

        std::cout<<"Updating Patches"<<std::endl;
        PTr.SetChoosenFromVertDir(TotVertIdx,TotVertDir,TotIsLoop);
        PTr.InitEdgeDirTable();
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

template <class MeshType>
void RecursiveProcess(PatchTracer<MeshType> &PTr,
                      const typename MeshType::ScalarType Drift,
                      bool onlyneeded,
                      bool finalremoval,
                      bool PreRemoveStep=true,
                      bool UseMetamesh=true)
{
    typedef typename MeshType::ScalarType ScalarType;

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

    std::cout<<"**** FIRST TRACING STEP ****"<<std::endl;
    PTr.BatchAddLoops(false,onlyneeded);
    int t1=clock();
    Time_FirstTrace+=t1-t0;

    //then do a first partitions update
    std::cout<<"Updating Patches"<<std::endl;
    PTr.UpdatePartitionsFromChoosen(true);
    std::cout<<"Updated"<<std::endl;

    //solve sub patches normally
    std::cout<<"**** FIRST SUBTRACING STEP ****"<<std::endl;
    SolveSubPatches(PTr,onlyneeded,false);

    //then check if there is some non-disk-like patches
    std::cout<<"**** CHECK NO DISK ONES ****"<<std::endl;
    SolveSubPatches(PTr,onlyneeded,true);

    std::vector<size_t> UnsolvedPartitionIndex;
    std::vector<PatchType> PatchTypes;
    std::cout<<"**** After All Insertion Steps ****"<<std::endl;

    std::cout<<"Updating Patches"<<std::endl;
    PTr.LazyUpdatePartitions();
    PTr.GetUnsolvedPartitionsIndex(UnsolvedPartitionIndex,PatchTypes);
    WriteUnsolvedStats(PatchTypes);
    std::cout<<"Updated"<<std::endl;

    int t2=clock();
    if (finalremoval)
    {
        std::cout<<"**** FINAL REMOVAL ****"<<std::endl;
        //PTr.UpdatePartitionsFromChoosen(true);
        PTr.SetAllRemovable();
        if (UseMetamesh)
            PTr.BatchRemovalMetaMesh(PreRemoveStep);
        else
            PTr.BatchRemovalOnMesh(PreRemoveStep);

        std::cout<<"**** After Last Removal Step ****"<<std::endl;
        PTr.GetUnsolvedPartitionsIndex(UnsolvedPartitionIndex,PatchTypes);
        WriteUnsolvedStats(PatchTypes);
    }
    else
        PTr.UpdatePartitionsFromChoosen(true);
    int t3=clock();
    Time_Removal+=t3-t2;

    std::cout<<"Updating Patches"<<std::endl;
    PTr.GetUnsolvedPartitionsIndex(UnsolvedPartitionIndex,PatchTypes);
    std::cout<<"**** FINAL THERE ARE "<<UnsolvedPartitionIndex.size()<<" Unsolved Partitions ****"<<std::endl;
    std::cout<<"**** TOTAL "<<PTr.Partitions.size()<<" Partitions ****"<<std::endl;
    std::cout<<"Updated"<<std::endl;

    std::cout<<"Smoothing"<<std::endl;
    PTr.SmoothPatches(10);
    std::cout<<"Fix Valences"<<std::endl;
    PTr.FixValences();
    PTr.WriteInfo();
    int t4=clock();
    ScalarType ElpsedSec=(ScalarType)(t4-t0)/CLOCKS_PER_SEC;
    std::cout<<"**** FINAL ELAPSED TIME "<<ElpsedSec<<" Seconds ****"<<std::endl;

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
            PInfo.NumPatchs,
            PInfo.HasEmit,
            PInfo.HighC,
            PInfo.LowC,
            PInfo.NonDiskLike,
            PInfo.SizePatches[0],
            PInfo.SizePatches[1],
            PInfo.SizePatches[2],
            PInfo.SizePatches[3],
            PInfo.SizePatches[4],
            PInfo.SizePatches[5],
            PInfo.SizePatches[6],
            PInfo.SizePatches[7]);
    fclose(f);
}


template <class MeshType>
void SaveAllData(PatchTracer<MeshType> &PTr,
                 const std::string &pathProject,
                 const size_t CurrNum)
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

    //merge vertices
    vcg::tri::Clean<MeshType>::RemoveDuplicateVertex(SaveM);
    vcg::tri::Clean<MeshType>::RemoveUnreferencedVertex(SaveM);
    vcg::tri::Allocator<MeshType>::CompactEveryVector(SaveM);
    SaveM.UpdateAttributes();
    SaveM.UpdateFromCoordPairs(SharpCoords);

    //save vert pos
    std::map<CoordType,size_t> VertMap;
    for (size_t i=0;i<SaveM.vert.size();i++)
        VertMap[SaveM.vert[i].P()]=i;

    for (size_t i=0;i<SaveM.face.size();i++)
        SaveM.face[i].Q()=PTr.FacePartitions[i];

    //vcg::tri::io::ExporterPLY<MeshType>::Save(SaveM,"test_final.ply");


    SplitAlongShap(SaveM);

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
    assert(PTr.PartitionCorners.size()==PTr.Partitions.size());
    fprintf(F,"%d\n",PTr.PartitionCorners.size());

    //std::cout<<"**** SAVING MESH ****"<<std::endl;
    PTr.WriteInfo();
    for (size_t i=0;i<PTr.PartitionCorners.size();i++)
    {
        fprintf(F,"%d\n",PTr.PartitionCorners[i].size());
        for (size_t j=0;j<PTr.PartitionCorners[i].size();j++)
        {
//            if (PTr.PartitionCorners[i].size()<MIN_ADMITTIBLE)
//            {
//                PTr.GetPatchMesh(i,mesh);
//                vcg::tri::io::ExporterPLY<MeshType>::Save(mesh,"")
//                assert(0);
//            }
            assert(PTr.PartitionCorners[i].size()>=MIN_ADMITTIBLE);
            assert(PTr.PartitionCorners[i].size()<=MAX_ADMITTIBLE);
            size_t IndexV=PTr.PartitionCorners[i][j];
            CoordType CornerPos=PTr.Mesh().vert[IndexV].P();
            assert(VertMap.count(CornerPos)>0);
            fprintf(F,"%d\n",VertMap[CornerPos]);
        }
    }
    fclose(F);

    std::string featurePartitions=pathProject;
    featurePartitions=featurePartitions+"_p"+std::to_string(CurrNum)+".feature";
    SaveM.SaveFeatures(featurePartitions);

    std::string featureCorners=pathProject;
    featureCorners=featureCorners+"_p"+std::to_string(CurrNum)+".c_feature";
    SaveM.SaveSharpCorners(featureCorners);
}


#endif
