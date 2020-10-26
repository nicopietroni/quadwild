#ifndef PATCH_TRACER_INTERFACE
#define PATCH_TRACER_INTERFACE


#include "patch_tracer.h"

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
                   bool inteleaveremoval,
                   bool finalremoval,
                   bool interleave_smooth)
{
    //first copy the submesh
    for (size_t i=0;i<PTr.Mesh().vert.size();i++)
        PTr.Mesh().vert[i].Q()=i;

    MeshType SubMesh;
    PTr.GetPatchMesh(IndexPatch,SubMesh);
    SubMesh.UpdateAttributes();

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
    PatchTracer<MeshType> SubTr(VFGraph);
    //SubTr.Init(PTr.Drift);
    SubTr.CopyParametersFrom(PTr);
    SubTr.CopyFrom(PTr,VertMap,IndexPatch);
    size_t Added_paths=SubTr.CopyPathsFrom(PTr,VertMap);
    if (Added_paths>0)
    {
        std::cout<<"ADDED ONE EXTRA PATH IN SUBDIVISION"<<std::endl;
    }
    //SubTr.UpdatePartitionsFromChoosen(true);
    //then trace in the subpatch
    //SubTr.BatchProcess(false,true);
    SubTr.DebugMsg=true;
    SubTr.BatchAddLoops(true,onlyneeded,inteleaveremoval,finalremoval,interleave_smooth);

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

    return (VertIdx.size()>0);
}


//template <class MeshType>
//bool TraceSubPatch3(const size_t &IndexPatch,
//                    PatchTracer<MeshType> &PTr,
//                    std::vector<std::vector<size_t> > &VertIdx,
//                    std::vector<std::vector<size_t> > &VertDir,
//                    std::vector<bool> &IsLoop)
//{
//    //first copy the submesh
//    for (size_t i=0;i<PTr.Mesh().vert.size();i++)
//        PTr.Mesh().vert[i].Q()=i;

//    MeshType SubMesh;
//    PTr.GetPatchMesh(IndexPatch,SubMesh);
//    SubMesh.UpdateAttributes();

//    std::vector<size_t> VertMap;
//    for (size_t i=0;i<SubMesh.vert.size();i++)
//        VertMap.push_back(SubMesh.vert[i].Q());

//    //copy original normals
//    for (size_t i=0;i<SubMesh.vert.size();i++)
//        SubMesh.vert[i].N()=PTr.Mesh().vert[VertMap[i]].N();

//    //make a subgraph
//    VertexFieldGraph<MeshType> VFGraph(SubMesh);
//    VFGraph.Init();

//    //initialize the tracer
//    PatchTracer<MeshType> SubTr(VFGraph);
//    SubTr.CopyParametersFrom(PTr);
//    //SubTr.Init(PTr.Drift);
//    SubTr.CopyFrom(PTr,VertMap,IndexPatch);

//    //then trace in the subpatch
//    //SubTr.BatchIterativeInsertion(true);
//    SubTr.BatchAddLoops(true,true,true,true);

//    //copy back paths to the original
//    SubTr.GetCurrVertDir(VertIdx,VertDir,IsLoop);
//    for (size_t i=0;i<VertIdx.size();i++)
//        for (size_t j=0;j<VertIdx[i].size();j++)
//            VertIdx[i][j]=VertMap[VertIdx[i][j]];

//    return (VertIdx.size()>0);
//}

//template <class MeshType>
//void RecursiveProcess3(PatchTracer<MeshType> &PTr,
//                       const typename MeshType::ScalarType Drift)
//{
//    //do a first step of tracing
//    PTr.Init(Drift);
//    PTr.BatchAddLoops(false,true,true,true);

//    PTr.UpdatePartitionsFromChoosen();
//    PTr.ColorByPartitions();
//    PTr.WriteInfo();


//    PTr.SmoothPatches(20);
//    //    return;

//    std::vector<std::vector<size_t> > TotVertIdx;
//    std::vector<std::vector<size_t> > TotVertDir;
//    std::vector<bool> TotIsLoop;
//    PTr.GetCurrVertDir(TotVertIdx,TotVertDir,TotIsLoop);

//    bool solved=false;
//    bool traced=false;

//    do{
//        std::vector<size_t> UnsolvedPartitionIndex;
//        PTr.GetUnsolvedPartitionsIndex(UnsolvedPartitionIndex);
//        if (UnsolvedPartitionIndex.size()==0)
//            solved=true;

//        traced=false;
//        std::cout<<"**** THERE ARE "<<UnsolvedPartitionIndex.size()<<" Unsolved Partitions ****"<<std::endl;
//        //        for(size_t i=8;i<9;i++)
//        //        {
//        //            for (size_t j=0;j<PTr.Partitions[UnsolvedPartitionIndex[i]].size();j++)
//        //            {
//        //               size_t IndexF=PTr.Partitions[UnsolvedPartitionIndex[i]][j];
//        //               PTr.Mesh().face[IndexF].C()=vcg::Color4b::Red;
//        //            }
//        //        }
//        for(size_t i=0;i<UnsolvedPartitionIndex.size();i++)
//        {
//            std::vector<std::vector<size_t> > NewVertIdx;
//            std::vector<std::vector<size_t> > NewVertDir;
//            std::vector<bool> NewIsLoop;

//            std::cout<<"**** SUB PATCH STEP ****"<<std::endl;
//            traced|=TraceSubPatch<MeshType>(UnsolvedPartitionIndex[i],PTr,NewVertIdx,NewVertDir,NewIsLoop,
//                                            true,true);

//            if (traced)
//                std::cout<<"TRACED"<<std::endl;
//            else
//                std::cout<<"NON TRACED"<<std::endl;

//            TotVertIdx.insert(TotVertIdx.end(),NewVertIdx.begin(),NewVertIdx.end());
//            TotVertDir.insert(TotVertDir.end(),NewVertDir.begin(),NewVertDir.end());
//            TotIsLoop.insert(TotIsLoop.end(),NewIsLoop.begin(),NewIsLoop.end());


//        }
//        if (traced)
//        {
//            std::cout<<"Updating Patches"<<std::endl;
//            PTr.SetChoosenFromVertDir(TotVertIdx,TotVertDir,TotIsLoop);
//            std::cout<<"Done Updating Patches"<<std::endl;
//        }

//    }while(traced & (!solved));

//    std::vector<size_t> UnsolvedPartitionIndex;
//    PTr.GetUnsolvedPartitionsIndex(UnsolvedPartitionIndex);
//    std::cout<<"**** FINAL THERE ARE "<<UnsolvedPartitionIndex.size()<<" Unsolved Partitions ****"<<std::endl;

//    PTr.SmoothPatches(20);
//    PTr.WriteInfo();
//    PTr.FixValences();
//    PTr.WriteInfo();
//}

void WriteUnsolvedStats(const std::vector<PatchType> &PatchTypes)
{
    size_t numLow=0;
    size_t numHigh=0;
    size_t numNonDisk=0;
    size_t numhasEmitter=0;
    size_t nummaxCCbility=0;
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
            default: numOK++;
        }
    }
    std::cout<<"** UNSATISFIED PATCHES **"<<std::endl;
    std::cout<<"*Low Sides:"<<numLow<<std::endl;
    std::cout<<"*High Sides:"<<numHigh<<std::endl;
    std::cout<<"*Non Disks:"<<numNonDisk<<std::endl;
    std::cout<<"*Has Emit:"<<numhasEmitter<<std::endl;
    std::cout<<"*Max CCbility:"<<nummaxCCbility<<std::endl;
    std::cout<<"*IsOk:"<<numOK<<std::endl;
}

template <class MeshType>
void RecursiveProcess(PatchTracer<MeshType> &PTr,
                       const typename MeshType::ScalarType Drift,
                       bool onlyneeded,
                       bool inteleaveremoval,
                       bool finalremoval,
                       bool interleave_smooth=false)
{
    //do it at the very end
//    bool InterleaveRemove=PTr.split_on_removal;
//    PTr.split_on_removal=false;
    //do a first step of tracing
    PTr.Init(Drift,true);
    //PTr.BatchProcess();
    PTr.BatchAddLoops(false,onlyneeded,inteleaveremoval,finalremoval,interleave_smooth);
    //PTr.BatchAddLoops(false,onlyneeded,inteleaveremoval,false,interleave_smooth);
    PTr.UpdatePartitionsFromChoosen(true);

    if (interleave_smooth)
        PTr.SmoothPatches(20);
    //return;
    std::vector<std::vector<size_t> > TotVertIdx;
    std::vector<std::vector<size_t> > TotVertDir;
    std::vector<bool> TotIsLoop;
    PTr.GetCurrVertDir(TotVertIdx,TotVertDir,TotIsLoop);

    bool solved=false;
    bool traced=false;
    do{
        std::vector<size_t> UnsolvedPartitionIndex;
        std::vector<PatchType> PatchTypes;
        //std::vector<PatchType> UnsolvedType;
        //PTr.GetUnsolvedPartitions(UnsolvedPartitions,UnsolvedType);
        PTr.UpdatePartitionsFromChoosen(true);
        PTr.GetUnsolvedPartitionsIndex(UnsolvedPartitionIndex,PatchTypes);
        //PTr.GetUnsolvedPartitions(UnsolvedPartitions,UnsolvedType);
        if (UnsolvedPartitionIndex.size()==0)
            solved=true;

        traced=false;
        //std::cout<<"**** THERE ARE "<<UnsolvedPartitionIndex.size()<<" Unsolved Partitions ****"<<std::endl;
        WriteUnsolvedStats(PatchTypes);

        for(size_t i=0;i<UnsolvedPartitionIndex.size();i++)
        {

            std::vector<std::vector<size_t> > NewVertIdx;
            std::vector<std::vector<size_t> > NewVertDir;
            std::vector<bool> NewIsLoop;

//            traced|=TraceSubPatch<MeshType>(UnsolvedPartitionIndex[i],PTr,
//                                            NewVertIdx,NewVertDir,
//                                            NewIsLoop,onlyneeded,
//                                            InterleaveRemove,false,
//                                            interleave_smooth);
            traced|=TraceSubPatch<MeshType>(UnsolvedPartitionIndex[i],PTr,
                                            NewVertIdx,NewVertDir,
                                            NewIsLoop,onlyneeded,
                                            inteleaveremoval,finalremoval,
                                            interleave_smooth);

            TotVertIdx.insert(TotVertIdx.end(),NewVertIdx.begin(),NewVertIdx.end());
            TotVertDir.insert(TotVertDir.end(),NewVertDir.begin(),NewVertDir.end());
            TotIsLoop.insert(TotIsLoop.end(),NewIsLoop.begin(),NewIsLoop.end());

        }

        if (traced)
        {
            std::cout<<"Updating Patches"<<std::endl;
            PTr.SetChoosenFromVertDir(TotVertIdx,TotVertDir,TotIsLoop);
            std::cout<<"Done Updating Patches"<<std::endl;
        }

    }while(traced & (!solved));

    if (finalremoval)
    {
        //PTr.split_on_removal=InterleaveRemove;
        PTr.UpdatePartitionsFromChoosen(true);
        PTr.BatchRemoval(interleave_smooth);
    }

    //PTr.FixValences();
    //PTr.WriteInfo();
    //std::cout<<"STEP 2"<<std::endl;
    std::vector<size_t> UnsolvedPartitionIndex;
    std::vector<PatchType> PatchTypes;
    PTr.UpdatePartitionsFromChoosen(true);
    PTr.GetUnsolvedPartitionsIndex(UnsolvedPartitionIndex,PatchTypes);
    std::cout<<"**** FINAL THERE ARE "<<UnsolvedPartitionIndex.size()<<" Unsolved Partitions ****"<<std::endl;
    std::cout<<"**** TOTAL "<<PTr.Partitions.size()<<" Partitions ****"<<std::endl;
    //PTr.WriteInfo();

    PTr.SmoothPatches(20);
    PTr.FixValences();
    PTr.WriteInfo();
//    PTr.split_on_removal=true;
//    PTr.BatchRemoval(false);
    //PTr.SmoothPatches(20);
    //PTr.WriteInfo();
}


//template <class MeshType>
//void UpdateTangentDirections(MeshType &mesh)
//{
//    typedef typename MeshType::CoordType CoordType;
//    typedef typename MeshType::ScalarType ScalarType;

//    //update field
//    for (size_t i=0;i<mesh.face.size();i++)
//    {
//        CoordType OldNorm=mesh.face[i].PD1()^mesh.face[i].PD2();
//        OldNorm.Normalize();
//        CoordType NewNorm=mesh.face[i].N();
//        vcg::Matrix33<ScalarType> M=vcg::RotationMatrix(OldNorm,NewNorm);
//        mesh.face[i].PD1()=M*mesh.face[i].PD1();
//        mesh.face[i].PD2()=mesh.face[i].N()^mesh.face[i].PD1();
//    }
//    for (size_t i=0;i<mesh.vert.size();i++)
//    {
//        CoordType OldNorm=mesh.vert[i].PD1()^mesh.vert[i].PD2();
//        OldNorm.Normalize();
//        CoordType NewNorm=mesh.vert[i].N();
//        vcg::Matrix33<ScalarType> M=vcg::RotationMatrix(OldNorm,NewNorm);
//        mesh.vert[i].PD1()=M*mesh.vert[i].PD1();
//        mesh.vert[i].PD1().Normalize();
//        mesh.vert[i].PD2()=mesh.vert[i].N()^mesh.vert[i].PD1();
//        mesh.vert[i].PD2().Normalize();
//    }
//}

//template <class MeshType>
//bool FindPaths(MeshType &mesh,
//               std::vector<size_t> &ValidFaces,
//               const typename MeshType::ScalarType Drift,
//               const typename MeshType::ScalarType Sample_Rate,
//               std::vector<std::vector<size_t> > &VertIdx,
//               std::vector<std::vector<size_t> > &VertDir,
//               std::vector<bool> &IsLoop,
//               std::vector<std::vector<size_t> > &UnsolvedPartitions,
//               std::vector<PatchType> &UnsolvedType)
//{
//    VertIdx.clear();
//    VertDir.clear();
//    UnsolvedPartitions.clear();
//    UnsolvedType.clear();

//    vcg::tri::UpdateSelection<MeshType>::Clear(mesh);
//    for (size_t i=0;i<ValidFaces.size();i++)
//        mesh.face[ValidFaces[i]].SetS();

//    for (size_t i=0;i<mesh.vert.size();i++)
//        mesh.vert[i].Q()=i;
//    for (size_t i=0;i<mesh.face.size();i++)
//        mesh.face[i].Q()=i;

//    vcg::tri::UpdateSelection<MeshType>::VertexFromFaceLoose(mesh);

//    MeshType curr_mesh;
//    vcg::tri::Append<MeshType,MeshType>::Mesh(curr_mesh,mesh,true);
//    curr_mesh.UpdateAttributes();
//    int num_splitted=vcg::tri::Clean<MeshType>::SplitNonManifoldVertex(curr_mesh,0);
//    if (num_splitted>0)
//        curr_mesh.UpdateAttributes();

//    //save indexes
//    std::vector<size_t> OriginalV,OriginalF;
//    for (size_t i=0;i<curr_mesh.vert.size();i++)
//        OriginalV.push_back(curr_mesh.vert[i].Q());

//    for (size_t i=0;i<curr_mesh.face.size();i++)
//        OriginalF.push_back(curr_mesh.face[i].Q());

//    //copy old normals
//    for (size_t i=0;i<curr_mesh.vert.size();i++)
//    {
//        size_t OrIndex=OriginalV[i];
//        curr_mesh.vert[i].N()=mesh.vert[OrIndex].N();
//    }

//    //then update the tangent
//    VertexFieldGraph<MeshType> VGraph(curr_mesh);
//    VGraph.Init();
//    PatchTracer<MeshType> PTr(VGraph);
//    PTr.sample_ratio=Sample_Rate;
//    PTr.Init(Drift);
//    PTr.BatchProcess();

//    //smooth
//    PTr.GetCurrVertDir(VertIdx,VertDir,IsLoop);
//    for (size_t i=0;i<VertIdx.size();i++)
//        for (size_t j=0;j<VertIdx[i].size();j++)
//            VertIdx[i][j]=OriginalV[VertIdx[i][j]];

//    //    if (VertIdx.size()>0)
//    //        PTr.SmoothPatches(10);

//    //copy smooth on the original
//    for (size_t i=0;i<curr_mesh.vert.size();i++)
//    {
//        size_t OrIndex=OriginalV[i];
//        mesh.vert[OrIndex].P()=curr_mesh.vert[i].P();
//        mesh.vert[OrIndex].PD1()=curr_mesh.vert[i].PD1();
//        mesh.vert[OrIndex].PD2()=curr_mesh.vert[i].PD2();
//    }

//    //finally recover the unsolved faces
//    PTr.GetUnsolvedPartitions(UnsolvedPartitions,UnsolvedType);

//    //convert to original index
//    for (size_t i=0;i<UnsolvedPartitions.size();i++)
//        for (size_t j=0;j<UnsolvedPartitions[i].size();j++)
//            UnsolvedPartitions[i][j]=OriginalF[UnsolvedPartitions[i][j]];

//    return (VertIdx.size()>0);
//}

//#define MAX_ITERATION 8

//template <class MeshType>
//size_t RecursiveProcessGeoSmooth(PatchTracer<MeshType> &PTr,const typename MeshType::ScalarType Drift)
//{
//    std::vector<size_t> ValidFaces(PTr.Mesh().face.size(),0);
//    for (size_t i=0;i<ValidFaces.size();i++)
//        ValidFaces[i]=i;

//    std::vector<std::vector<size_t> > TotVertIdx,TotVertDir;
//    std::vector<bool> TotLoop;
//    size_t performed_step=0;
//    bool HasAdded=false;
//    std::vector<std::vector<size_t> > TotUnsolvedPartitions;
//    bool has_completed=false;
//    do
//    {
//        std::vector<std::vector<size_t> > VertIdx,VertDir;
//        std::vector<bool> IsLoop;
//        std::vector<std::vector<size_t> > UnsolvedPartitions;
//        std::vector<PatchType> UnsolvedType;
//        FindPaths(PTr.Mesh(),ValidFaces,Drift,PTr.sample_ratio,
//                  VertIdx,VertDir,IsLoop,UnsolvedPartitions,
//                  UnsolvedType);

//        TotVertIdx.insert(TotVertIdx.end(),VertIdx.begin(),VertIdx.end());
//        TotVertDir.insert(TotVertDir.end(),VertDir.begin(),VertDir.end());
//        TotLoop.insert(TotLoop.end(),IsLoop.begin(),IsLoop.end());

//        HasAdded=(VertIdx.size()>0);
//        //if added something then update
//        if (HasAdded)
//            TotUnsolvedPartitions.insert(TotUnsolvedPartitions.end(),
//                                         UnsolvedPartitions.begin(),
//                                         UnsolvedPartitions.end());

//        performed_step++;
//        std::cout<<"****** PERFORMED "<<performed_step<<" ******"<<std::endl;
//        std::cout<<"* Added "<<VertIdx.size()<<" PATHS"<<std::endl;
//        if (VertIdx.size()==1)
//        {
//            std::cout<<"* Size "<<VertIdx[0].size()<<" nodes"<<std::endl;
//            std::cout<<"* V0 "<<VertIdx[0][0]<<std::endl;
//            std::cout<<"* V1 "<<VertIdx[0][1]<<std::endl;
//            size_t VIndex0=VertIdx[0][0];
//            size_t VIndex1=VertIdx[0][1];
//            std::cout<<"* Coord0 "<<PTr.Mesh().vert[VIndex0].P().X()
//                    <<","<<PTr.Mesh().vert[VIndex0].P().Y()
//                   <<","<<PTr.Mesh().vert[VIndex0].P().Z()<<std::endl;
//            std::cout<<"* Coord1 "<<PTr.Mesh().vert[VIndex1].P().X()
//                    <<","<<PTr.Mesh().vert[VIndex1].P().Y()
//                   <<","<<PTr.Mesh().vert[VIndex1].P().Z()<<std::endl;
//        }
//        std::cout<<"* Unsolved "<<UnsolvedPartitions.size()<<" partitions "<<std::endl;
//        if ((UnsolvedPartitions.size()==1)&&(VertIdx.size()==1))
//        {
//            std::cout<<"* Size "<<UnsolvedPartitions[0].size()<<" faces"<<std::endl;
//            vcg::tri::UpdateSelection<MeshType>::FaceClear(PTr.Mesh());
//            vcg::tri::UpdateSelection<MeshType>::VertexClear(PTr.Mesh());
//            for(size_t i=0;i<UnsolvedPartitions[0].size();i++)
//                PTr.Mesh().face[UnsolvedPartitions[0][i]].SetS();

//            size_t VIndex0=VertIdx[0][0];
//            size_t VIndex1=VertIdx[0][1];
//            PTr.Mesh().vert[VIndex0].SetS();
//            PTr.Mesh().vert[VIndex1].SetS();
//            vcg::tri::io::ExporterPLY<MeshType>::Save(PTr.Mesh(),"test.ply",vcg::tri::io::Mask::IOM_FACEFLAGS|
//                                                      vcg::tri::io::Mask::IOM_VERTFLAGS);
//        }

//        if (TotUnsolvedPartitions.size()>0)
//        {
//            ValidFaces=TotUnsolvedPartitions.back();
//            TotUnsolvedPartitions.pop_back();
//        }
//        else has_completed=true;

//    }while ((!has_completed)&&(performed_step<MAX_ITERATION));

//    //finally reasseble the result in a single tracer
//    PTr.Init(Drift);
//    //PTr.WriteInfo();
//    std::cout<<"****** REASSEMBLING ******"<<std::endl;
//    PTr.SetChoosenFromVertDir(TotVertIdx,TotVertDir,TotLoop);
//    PTr.WriteInfo();
//    //PTr.FixValences();
//    PTr.WriteInfo();
//    return performed_step;
//}

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
    fprintf(f,"%s %d %d %d %d %d %d %d %d %d %d %d %d %d \n",pathProject.c_str(),
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


//template <class MeshType>
//void SaveAllData(PatchTracer<MeshType> &PTr,
//                 const std::string &pathProject,
//                 const size_t CurrNum)
//{
//    typedef typename MeshType::CoordType CoordType;

////    std::vector<std::pair<CoordType,CoordType> > SharpCoords;
////    PTr.Mesh().GetSharpCoordPairs(SharpCoords);

////    //copy the mesh
////    MeshType SaveM;
////    vcg::tri::Append<MeshType,MeshType>::Mesh(SaveM,PTr.Mesh());
////    std::vector<size_t> SharpCorners;
////    PTr.getCornerSharp(SharpCorners);
////    std::set<CoordType> SharpCornerPos;
////    for(size_t i=0;i<SharpCorners.size();i++)
////        SharpCornerPos.insert(PTr.Mesh().vert[SharpCorners[i]].P());

////    //merge vertices
////    vcg::tri::Clean<MeshType>::RemoveDuplicateVertex(SaveM);
////    vcg::tri::Clean<MeshType>::RemoveUnreferencedVertex(SaveM);
////    vcg::tri::Allocator<MeshType>::CompactEveryVector(SaveM);
////    SaveM.UpdateAttributes();
////    SaveM.UpdateFromCoordPairs(SharpCoords);

////    //save vert pos
////    std::map<CoordType,size_t> VertMap;
////    for (size_t i=0;i<SaveM.vert.size();i++)
////        VertMap[SaveM.vert[i].P()]=i;

////    for (size_t i=0;i<SaveM.face.size();i++)
////        SaveM.face[i].Q()=PTr.FacePartitions[i];

////    SplitAlongShap(SaveM);

////    //update sharp vertices
////    SaveM.SharpCorners.clear();
////    for (size_t i=0;i<SaveM.vert.size();i++)
////        if (SharpCornerPos.count(SaveM.vert[i].P())>0)
////            SaveM.SharpCorners.push_back(i);

//    //save the mesh
//    int Mask=0;
//    std::string pathMeshFinal=pathProject;
//    pathMeshFinal=pathMeshFinal+"_p"+std::to_string(CurrNum)+".obj";
//    vcg::tri::io::ExporterOBJ<MeshType>::Save(PTr.Mesh(),pathMeshFinal.c_str(),Mask);

//    std::string pathPartitions=pathProject;
//    //pathPartitions.append("_p.patch");
//    pathPartitions=pathPartitions+"_p"+std::to_string(CurrNum)+".patch";
//    FILE *F=fopen(pathPartitions.c_str(),"wt");
//    assert(F!=NULL);
//        assert(PTr.FacePartitions.size()==PTr.Mesh().face.size());
//        fprintf(F,"%d\n",PTr.FacePartitions.size());
//        for (size_t i=0;i<PTr.FacePartitions.size();i++)
//            fprintf(F,"%d\n",PTr.FacePartitions[i]);
//        fclose(F);
////    fprintf(F,"%d\n",SaveM.face.size());
////    for (size_t i=0;i<SaveM.face.size();i++)
////        fprintf(F,"%d\n",(int)SaveM.face[i].Q());
////    fclose(F);

//    std::string pathCorners=pathProject;
//    //pathCorners.append("_p.corners");
//    pathCorners=pathCorners+"_p"+std::to_string(CurrNum)+".corners";
//    F=fopen(pathCorners.c_str(),"wt");
//    assert(F!=NULL);
//    assert(PTr.PartitionCorners.size()==PTr.Partitions.size());
//    fprintf(F,"%d\n",PTr.PartitionCorners.size());
//    for (size_t i=0;i<PTr.PartitionCorners.size();i++)
//    {
//        fprintf(F,"%d\n",PTr.PartitionCorners[i].size());
//        for (size_t j=0;j<PTr.PartitionCorners[i].size();j++)
//        {
//            size_t IndexV=PTr.PartitionCorners[i][j];
//            //CoordType CornerPos=PTr.Mesh().vert[IndexV].P();
//            //assert(VertMap.count(CornerPos)>0);
//            fprintf(F,"%d\n",IndexV);
//        }
//    }
//    fclose(F);

//    std::string featurePartitions=pathProject;
//    featurePartitions=featurePartitions+"_p"+std::to_string(CurrNum)+".feature";
//    PTr.Mesh().SaveFeatures(featurePartitions);

//    std::string featureCorners=pathProject;
//    featureCorners=featureCorners+"_p"+std::to_string(CurrNum)+".c_feature";
//    PTr.Mesh().SaveSharpCorners(featureCorners);
//}

#endif
