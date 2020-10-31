#ifndef QR_PATCH_ASSEMBLER_H
#define QR_PATCH_ASSEMBLER_H

#include <vcg/complex/complex.h>
#include <vcg/math/matrix33.h>
#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>
#include <vcg/complex/algorithms/isotropic_remeshing.h>
#include <vcg/complex/algorithms/geodesic.h>
#include <wrap/io_trimesh/export.h>
#include <vcg/complex/algorithms/polygonal_algorithms.h>
#include <vcg/complex/algorithms/hole.h>
#include <vcg/complex/algorithms/update/curvature.h>

#include "qr_field_tracer.h"
#include "qr_utils.h"
#include "qr_field_smoother.h"

#define CONVEX 8.0
#define CONCAVE 8.0

namespace QuadRetopology {
namespace internal {


template <class MeshType>
class PatchSplitter
{
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;

    MeshType &patchMesh;

    bool printMsg;
    VertexFieldTracer<MeshType> &VFTracer;
    std::vector<std::vector<size_t> > TraceVertCandidates,TraceDirCandidates;
    std::vector<std::pair<ScalarType,size_t> > CandidatesPathLenghts;

    void GetConvexVertices(const std::set<size_t> &ConvexOriginal,
                           std::set<size_t> &ConvexMesh)
    {
        for (size_t  i=0;i<patchMesh.vert.size();i++)
        {
            size_t IndexOriginal=patchMesh.vert[i].Q();
            if(ConvexOriginal.count(IndexOriginal)==0)continue;
            ConvexMesh.insert(i);
        }
    }

    void GetConcaveVertices(const std::set<size_t> &ConcaveOriginal,
                            std::set<size_t> &ConcaveMesh)
    {
        for (size_t  i=0;i<patchMesh.vert.size();i++)
        {
            size_t IndexOriginal=patchMesh.vert[i].Q();
            if(ConcaveOriginal.count(IndexOriginal)==0)continue;
            ConcaveMesh.insert(i);
        }
    }

    void InitTracerConnections(const std::set<size_t> &ConvexOriginal)
    {
        std::set<size_t> ConvexMesh;
        GetConvexVertices(ConvexOriginal,ConvexMesh);
        VFTracer.InitConnections(ConvexMesh);
    }

    void InitCandidatesPathLenghts()
    {
        CandidatesPathLenghts.clear();
        for (size_t i=0;i<TraceVertCandidates.size();i++)
        {
            ScalarType currL=VFTracer.TraceLenght(TraceVertCandidates[i],TraceDirCandidates[i]);
            CandidatesPathLenghts.push_back(std::pair<ScalarType,size_t>(currL,i));
        }
        std::sort(CandidatesPathLenghts.begin(),CandidatesPathLenghts.end());
    }

    void PruneConcaveCandidates(bool addAll=false)
    {
        std::vector<std::vector<size_t > > TraceVertSwap;
        std::vector<std::vector<size_t > > TraceDirSwap;

        //count how many traces has been done for each concave
        std::vector<size_t> VertAllocator(patchMesh.vert.size(),1);
        InitCandidatesPathLenghts();

        //first round
        int Num0=TraceVertCandidates.size();
        for (size_t i=0;i<CandidatesPathLenghts.size();i++)
        {
            //get the current trace from the sorted ones
            size_t TraceIndex0=CandidatesPathLenghts[i].second;
            std::vector<size_t > TraceV0=TraceVertCandidates[TraceIndex0];
            std::vector<size_t > TraceD0=TraceDirCandidates[TraceIndex0];

            //get the first vertex to check if it has been already traced or not
            size_t IndexV0=TraceV0[0];

            //if it has already a trace then go on
            if (VertAllocator[IndexV0]==0)continue;

            bool collide=false;
            for (size_t i=0;i<TraceVertSwap.size();i++)
            {
                std::vector<size_t > TraceV1=TraceVertSwap[i];
                std::vector<size_t > TraceD1=TraceDirSwap[i];

                collide|=VFTracer.CollideTraces(TraceV0,TraceD0,TraceV1,TraceD1);
                if (collide)break;
            }
            if (!collide)
            {
                TraceVertSwap.push_back(TraceV0);
                TraceDirSwap.push_back(TraceD0);
                VertAllocator[IndexV0]--;
            }
        }

        //second round
        if (addAll)
        {
            VertAllocator.clear();
            VertAllocator.resize(patchMesh.vert.size(),1);
            for (size_t i=0;i<TraceVertCandidates.size();i++)
            {
                size_t TraceIndex0=CandidatesPathLenghts[i].second;
                std::vector<size_t > TraceV0=TraceVertCandidates[TraceIndex0];
                std::vector<size_t > TraceD0=TraceDirCandidates[TraceIndex0];
                size_t IndexV0=TraceV0[0];
                if (VertAllocator[IndexV0]==0)continue;

                bool collide=false;
                for (size_t i=0;i<TraceVertSwap.size();i++)
                {
                    std::vector<size_t > TraceV1=TraceVertSwap[i];
                    std::vector<size_t > TraceD1=TraceDirSwap[i];

                    collide|=VFTracer.CollideTraces(TraceV0,TraceD0,TraceV1,TraceD1);
                    if (collide)break;
                }
                if (!collide)
                {
                    TraceVertSwap.push_back(TraceV0);
                    TraceDirSwap.push_back(TraceD0);
                    VertAllocator[IndexV0]--;
                }
            }
        }
        TraceVertCandidates=TraceVertSwap;
        TraceDirCandidates=TraceDirSwap;
        int Num1=TraceVertCandidates.size();

        if (printMsg) {
            std::cout<<"Pruned "<<Num0-Num1<<" out ouf "<<Num0<<std::endl;
        }
        //InitPathPriority();
    }

    void TraceConcaveCandidates(const std::set<size_t> &ConvexOriginal,
                                const std::set<size_t> &ConcaveOriginal,
                                const std::vector<size_t> &TracingVert,
                                const std::vector<std::vector<size_t> > &TracingDir)
    {
        if (printMsg) {
            std::cout<<"Tracing Concave Corners"<<std::endl;
        }

        assert(TracingVert.size()==TracingDir.size());
        InitTracerConnections(ConvexOriginal);

        std::set<size_t> ConcaveMesh;
        GetConcaveVertices(ConcaveOriginal,ConcaveMesh);

        size_t tested=0;
        size_t non_expanded=0;
        size_t non_traced=0;

        TraceVertCandidates.clear();
        TraceDirCandidates.clear();
        //vcg::tri::UpdateSelection<MeshType>::VertexClear(patchMesh);
        for (size_t i=0;i<TracingVert.size();i++)
        {
            size_t TraceV=TracingVert[i];
            if (ConcaveMesh.count(TraceV)==0)continue;
            //std::cout<<"*** Size Start "<<TracingDir[i].size()<<std::endl;

            //patchMesh.vert[TraceV].SetS();
            for (size_t j=0;j<TracingDir[i].size();j++)
            {
                tested++;
                std::vector<size_t> IndexV,IndexDir;
                bool traced=VFTracer.TraceFrom(TraceV,
                                               TracingDir[i][j],
                                               IndexV,IndexDir);
                if (!traced)
                {
                    non_traced++;

                    if (printMsg) {
                        std::cout<<"WARNING: Non Traced "<<TraceV<<" "<<TracingDir[i][j]<<std::endl;
                    }
                    continue;
                }

                bool expanded=VFTracer.ExpandPath(IndexV,IndexDir);
                if (!expanded)
                {
                    non_expanded++;
                    if (printMsg) {
                        std::cout<<"WARNING: Non Expanded "<<TraceV<<" "<<TracingDir[i][j]<<std::endl;
                    }
                    continue;
                }

                TraceVertCandidates.push_back(IndexV);
                TraceDirCandidates.push_back(IndexDir);
            }
        }
        if (printMsg)
        {
            std::cout<<"Non Traced "<<non_traced<<std::endl;
            std::cout<<"Non Expanded "<<non_expanded<<std::endl;
            std::cout<<"Tested "<<tested<<std::endl;
        }

#ifdef SAVE_MESHES_FOR_DEBUG
        vcg::tri::io::ExporterOBJ<MeshType>::Save(patchMesh,"results/tracer_mesh_trace.obj",vcg::tri::io::Mask::IOM_VERTFLAGS);
#endif
    }

    void TraceAllBordersCandidates(const std::set<size_t> &ConvexOriginal,
                                   const std::vector<size_t> &TracingVert,
                                   const std::vector<std::vector<size_t> > &TracingDir)
    {
        if (printMsg) {
            std::cout<<"Tracing Regular Paths"<<std::endl;
        }

        assert(TracingVert.size()==TracingDir.size());
        InitTracerConnections(ConvexOriginal);

        size_t tested=0;
        size_t non_expanded=0;
        size_t non_traced=0;

        TraceVertCandidates.clear();
        TraceDirCandidates.clear();
        //vcg::tri::UpdateSelection<MeshType>::VertexClear(patchMesh);
        for (size_t i=0;i<TracingVert.size();i++)
        {
            //std::cout<<"Tracing From Vertex "<<TracingVert[i]<<std::endl;

            size_t TraceV=TracingVert[i];
            if (TracingDir[i].size()==0)continue;
            //if it is flat only one direction is feasible
            assert(TracingDir[i].size()<2);

            tested++;
            std::vector<size_t> IndexV,IndexDir;

            //std::cout<<"Test Trace "<<TracingVert[i]<<std::endl;

            bool traced=VFTracer.TraceFrom(TraceV,TracingDir[i][0],
                    IndexV,IndexDir);
            //std::cout<<"End Trace "<<TracingVert[i]<<std::endl;

            if (!traced)
            {
                non_traced++;
                if (printMsg) {
                    std::cout<<"WARNING: Non Traced "<<TraceV<<" "<<TracingDir[i][0]<<std::endl;
                }
                continue;
            }

            bool expanded=VFTracer.ExpandPath(IndexV,IndexDir);
            if (!expanded)
            {
                non_expanded++;
                if (printMsg) {
                    std::cout<<"WARNING: Non Expanded "<<TraceV<<" "<<TracingDir[i][0]<<std::endl;
                }
                continue;
            }

            TraceVertCandidates.push_back(IndexV);
            TraceDirCandidates.push_back(IndexDir);
        }
        if (printMsg)
        {
            std::cout<<"Non Traced "<<non_traced<<std::endl;
            std::cout<<"Non Expanded "<<non_expanded<<std::endl;
            std::cout<<"Tested "<<tested<<std::endl;
        }
    }

    void RetrievePartitioningFrom(const std::set<std::pair<CoordType,CoordType> > BorderEdges,
                                  const size_t &IndexF,std::vector<size_t> &partition)
    {
        partition.clear();

        std::vector<size_t> stack;
        std::set<size_t> explored;

        stack.push_back(IndexF);
        explored.insert(IndexF);
        do
        {
            size_t currF=stack.back();
            stack.pop_back();

            partition.push_back(currF);
            for (size_t i=0;i<patchMesh.face[currF].VN();i++)
            {
                if (vcg::face::IsBorder(patchMesh.face[currF],i))continue;

                CoordType Pos0=patchMesh.face[currF].P0(i);
                CoordType Pos1=patchMesh.face[currF].P1(i);

                std::pair<CoordType,CoordType> key(std::min(Pos0,Pos1),
                                                   std::max(Pos0,Pos1));

                if (BorderEdges.count(key)==1)continue;

                int NextFIndex=vcg::tri::Index(patchMesh,patchMesh.face[currF].FFp(i));

                if (explored.count(NextFIndex)>0)continue;

                explored.insert(NextFIndex);
                stack.push_back(NextFIndex);
            }
        }while (!stack.empty());
    }

    void GetBorderSet(std::set<std::pair<CoordType,CoordType> > &BorderEdges)
    {
        BorderEdges.clear();
        for (size_t i=0;i<TraceVertCandidates.size();i++)
        {
            if (TraceVertCandidates[i].size()==0)continue;
            for (size_t j=0;j<TraceVertCandidates[i].size()-1;j++)
            {
                size_t IndexV0=TraceVertCandidates[i][j];
                size_t IndexV1=TraceVertCandidates[i][j+1];

                CoordType PosV0=patchMesh.vert[IndexV0].P();
                CoordType PosV1=patchMesh.vert[IndexV1].P();

                std::pair<CoordType,CoordType> Key(std::min(PosV0,PosV1),
                                                   std::max(PosV0,PosV1));
                BorderEdges.insert(Key);
            }
        }
    }

    void RetrievePatchesFromPaths(std::vector<std::vector<size_t> > &Partitions)
    {
        std::set<std::pair<CoordType,CoordType> > BorderEdges;
        GetBorderSet(BorderEdges);

        Partitions.clear();
        vcg::tri::UpdateFlags<MeshType>::FaceClearV(patchMesh);
        for (size_t i=0;i<patchMesh.face.size();i++)
        {
            if (patchMesh.face[i].IsV())continue;

            //std::cout<<"test"<<std::endl;
            std::vector<size_t> partition;
            RetrievePartitioningFrom(BorderEdges,i,partition);

            for (size_t j=0;j<partition.size();j++)
                patchMesh.face[partition[j]].SetV();

            Partitions.push_back(partition);
        }
    }

    void SavePatchMesh(std::vector<std::vector<size_t> > &FacePatches,std::set<size_t> &SelV)
    {
        for (size_t i=0;i<FacePatches.size();i++)
        {
            for (size_t j=0;j<FacePatches[i].size();j++)
            {
                size_t IndexF=FacePatches[i][j];
                vcg::Color4b PatchCol=vcg::Color4b::Scatter(FacePatches.size(),i);
                patchMesh.face[IndexF].C()=PatchCol;
            }
        }

        vcg::tri::UpdateSelection<MeshType>::VertexClear(patchMesh);
        for (size_t i=0;i<patchMesh.vert.size();i++)
        {
            if (SelV.count(i)==0)continue;
            patchMesh.vert[i].SetS();
        }

#ifdef SAVE_MESHES_FOR_DEBUG
        vcg::tri::io::ExporterOBJ<MeshType>::Save(patchMesh,"results/tracer_mesh.obj",
                                                  vcg::tri::io::Mask::IOM_FACECOLOR|
                                                  vcg::tri::io::Mask::IOM_VERTFLAGS);
#endif

    }

    void SavePatchMeshConcaveQ(const std::set<size_t> &ConcaveOriginal)
    {
        vcg::tri::UpdateSelection<MeshType>::VertexClear(patchMesh);
        for (size_t i=0;i<patchMesh.vert.size();i++)
        {
            size_t OriginalQ=patchMesh.vert[i].Q();
            if (ConcaveOriginal.count(OriginalQ)==0)continue;
            patchMesh.vert[i].SetS();
        }

#ifdef SAVE_MESHES_FOR_DEBUG
        vcg::tri::io::ExporterOBJ<MeshType>::Save(patchMesh,"results/tracer_mesh_flag1.obj",vcg::tri::io::Mask::IOM_VERTFLAGS);
#endif
    }

    void SavePatchMeshConcave(const std::set<size_t> &ConcaveOriginal)
    {
        std::set<size_t> ConcaveMesh;
        GetConcaveVertices(ConcaveOriginal,ConcaveMesh);

        vcg::tri::UpdateSelection<MeshType>::VertexClear(patchMesh);
        for (size_t i=0;i<patchMesh.vert.size();i++)
        {
            if (ConcaveMesh.count(i)==0)continue;
            patchMesh.vert[i].SetS();
        }

#ifdef SAVE_MESHES_FOR_DEBUG
        vcg::tri::io::ExporterOBJ<MeshType>::Save(patchMesh,"results/tracer_mesh_flag2.obj",vcg::tri::io::Mask::IOM_VERTFLAGS);
#endif
    }

    //    CoordType MakeOrthogonalTo(const CoordType &Dir,
    //                               const CoordType &Norm)
    //    {
    //        CoordType TestDir=Dir;
    //        //check if orthogonal
    //        if(fabs(TestDir*Norm)<0.01)return Dir;

    //        CoordType RotAxis=Norm^TestDir;
    //        CoordType TargetN=TestDir^RotAxis;
    //        TargetN.Normalize();
    //        vcg::Matrix33<ScalarType> rot=vcg::RotationMatrix(TargetN,Norm);
    //        TestDir=rot*TestDir;
    //        TestDir.Normalize();
    //        assert(fabs(TestDir*Norm)<0.01);
    //        return TestDir;
    //    }

    void PruneCollidingCandidates()
    {
        std::vector<std::vector<size_t > > TraceVertSwap;
        std::vector<std::vector<size_t > > TraceDirSwap;

        InitCandidatesPathLenghts();
        for (size_t i=0;i<TraceVertCandidates.size();i++)
        {
            size_t TraceIndex0=CandidatesPathLenghts[i].second;
            std::vector<size_t > TraceV0=TraceVertCandidates[TraceIndex0];
            std::vector<size_t > TraceD0=TraceDirCandidates[TraceIndex0];
            bool collide=false;
            for (size_t i=0;i<TraceVertSwap.size();i++)
            {
                std::vector<size_t > TraceV1=TraceVertSwap[i];
                std::vector<size_t > TraceD1=TraceDirSwap[i];

                collide|=VFTracer.CollideTraces(TraceV0,TraceD0,TraceV1,TraceD1);
                if (collide)break;
            }
            if (!collide)
            {
                TraceVertSwap.push_back(TraceV0);
                TraceDirSwap.push_back(TraceD0);
            }
        }
        TraceVertCandidates=TraceVertSwap;
        TraceDirCandidates=TraceDirSwap;
    }

    bool IsOkPatchHoles(const std::vector<size_t> &FaceIndexes)
    {
        //get the mesh
        MeshType testMesh;
        vcg::tri::UpdateSelection<MeshType>::Clear(patchMesh);
        for (size_t i=0;i<FaceIndexes.size();i++)
            patchMesh.face[FaceIndexes[i]].SetS();
        vcg::tri::UpdateSelection<MeshType>::VertexFromFaceLoose(patchMesh);
        vcg::tri::Append<MeshType,MeshType>::Mesh(testMesh,patchMesh,true);
        updateAllMeshAttributes(testMesh);
        vcg::tri::UpdateSelection<MeshType>::Clear(patchMesh);
        vcg::tri::UpdateSelection<MeshType>::Clear(testMesh);

        return (vcg::tri::Clean<MeshType>::CountHoles(testMesh)==1);
    }

    bool IsOkPatchCorner(const std::vector<size_t> &FaceIndexes,
                         const std::set<size_t> &CurrVertices)
    {
        //retrieve the vertices
        std::set<size_t> Vertices;
        for (size_t i=0;i<FaceIndexes.size();i++)
        {
            size_t IndexF=FaceIndexes[i];
            for (size_t j=0;j<3;j++)
            {
                size_t IndexV=vcg::tri::Index(patchMesh,patchMesh.face[IndexF].V(j));
                Vertices.insert(IndexV);
            }
        }
        std::set<int> intersect;
        std::set_intersection(Vertices.begin(),Vertices.end(),CurrVertices.begin(),CurrVertices.end(),
                              std::inserter(intersect,intersect.begin()));

        //std::cout<<"Num Sides "<<intersect.size()<<std::endl;

        return ((intersect.size()>=3)&&(intersect.size()<=6));
    }

    void CornersFromCurrentPaths(std::set<size_t> &CurrV)
    {
        std::vector<size_t> NumTraces(patchMesh.vert.size(),0);
        for (size_t i=0;i<TraceVertCandidates.size();i++)
            for (size_t j=0;j<TraceVertCandidates[i].size();j++)
            {
                size_t IndexV=TraceVertCandidates[i][j];
                NumTraces[IndexV]++;
            }
        CurrV.clear();
        for (size_t i=0;i<NumTraces.size();i++)
        {
            if (patchMesh.vert[i].IsB() && NumTraces[i]>0)
            {
                CurrV.insert(i);
                continue;
            }
            if (NumTraces[i]>=2)
                CurrV.insert(i);
        }
        CurrV.insert(ConvexPatchMesh.begin(),ConvexPatchMesh.end());
    }

    bool CanRemove(size_t IndexCandidate,size_t MinPatches)
    {
        assert(IndexCandidate<TraceVertCandidates.size());

        std::vector<size_t > TraceVertSwap=TraceVertCandidates[IndexCandidate];
        std::vector<size_t > TraceDirSwap=TraceDirCandidates[IndexCandidate];

        //get current patches before removal
        std::vector<std::vector<size_t> > CurrPatches;
        RetrievePatchesFromPaths(CurrPatches);
        if (CurrPatches.size()<=MinPatches)return false;

        //then update the new corners
        std::set<size_t> CurrV;
        CornersFromCurrentPaths(CurrV);
        //std::cout<<"there are "<<CurrV.size()<<" vertices"<<std::endl;

        //        SavePatchMesh(CurrPatches,CurrV);
        //        exit(0);
        //then count the ones non ok before the remove
        size_t NonOkPatches=0;
        for (size_t i=0;i<CurrPatches.size();i++)
        {
            if (!IsOkPatchCorner(CurrPatches[i],CurrV)){NonOkPatches+=CurrPatches[i].size();continue;}
            if (!IsOkPatchHoles(CurrPatches[i])){NonOkPatches+=CurrPatches[i].size();continue;}
        }

        //remove the current path
        TraceVertCandidates[IndexCandidate].clear();
        TraceDirCandidates[IndexCandidate].clear();

        //get new patches after removal
        std::vector<std::vector<size_t> > NewPatches;
        RetrievePatchesFromPaths(NewPatches);

        //then update the new corners
        std::set<size_t> CurrVNew;
        CornersFromCurrentPaths(CurrVNew);

        //then count the ones non ok before the remove
        size_t NonOkPatchesNew=0;
        for (size_t i=0;i<NewPatches.size();i++)
        {
            if (!IsOkPatchCorner(NewPatches[i],CurrVNew)){NonOkPatchesNew+=NewPatches[i].size();continue;}
            if (!IsOkPatchHoles(NewPatches[i])){NonOkPatchesNew+=NewPatches[i].size();continue;}
        }
        TraceVertCandidates[IndexCandidate]=TraceVertSwap;
        TraceDirCandidates[IndexCandidate]=TraceDirSwap;

        //        std::cout<<"Before: "<<NonOkPatches<<std::endl;
        //        std::cout<<"After: "<<NonOkPatchesNew<<std::endl;

        //exit(0);

        return (NonOkPatchesNew<=NonOkPatches);
    }

public:

    void GetVertexTopology(const std::vector<std::vector<size_t> > &NewPatches,
                           const std::set<size_t> &ConvexOriginal,
                           const std::set<size_t> &ConcaveOriginal,
                           std::vector<std::vector<size_t> > &ConvexPatchVert,
                           std::vector<std::vector<size_t> > &ConcavePatchVert)
    {
        //count the number of pathes per edge
        std::vector<std::vector<size_t> > VertPatchCount;
        std::vector<std::vector<ScalarType> > VertPatchAngle;

        VertPatchCount.resize(patchMesh.vert.size());
        VertPatchAngle.resize(patchMesh.vert.size());
        for (size_t i=0;i<NewPatches.size();i++)
        {
            for (size_t j=0;j<NewPatches[i].size();j++)
            {
                size_t IndexF=NewPatches[i][j];
                assert(IndexF>=0);
                assert(IndexF<patchMesh.face.size());
                for (size_t k=0;k<patchMesh.face[IndexF].VN();k++)
                {
                    //compute the angle for the patch
                    // VertexType *VCurr=patchMesh.face[IndexF].V(k);
                    // VertexType *V0=patchMesh.face[IndexF].V1(k);
                    // VertexType *V1=patchMesh.face[IndexF].V2(k);
                    // CoordType Dir0=V0->P()-VCurr->P();
                    // CoordType Dir1=V1->P()-VCurr->P();
                    // Dir0.Normalize();
                    // Dir1.Normalize();
                    // Dir0=MakeOrthogonalTo(Dir0,VCurr->N());
                    // Dir1=MakeOrthogonalTo(Dir1,VCurr->N());
                    // ScalarType Angle=vcg::Angle(Dir0,Dir1);
                    ScalarType Angle=vcg::face::WedgeAngleRad(patchMesh.face[IndexF],k);//vcg::Angle(Dir0,Dir1);
                    //assert(Angle>0);

                    //then add the angle per vert per patch
                    size_t IndexV=vcg::tri::Index(patchMesh,patchMesh.face[IndexF].V(k));

                    //see if there is the patch already
                    std::vector<size_t>::iterator it;
                    it = find (VertPatchCount[IndexV].begin(),VertPatchCount[IndexV].end(),i);
                    if (it == VertPatchCount[IndexV].end())
                    {
                        VertPatchCount[IndexV].push_back(i);
                        VertPatchAngle[IndexV].push_back(Angle);
                    }
                    else
                    {
                        size_t Index=std::distance(VertPatchCount[IndexV].begin(),it);
                        assert(VertPatchCount[IndexV][Index]==i);
                        VertPatchAngle[IndexV][Index]+=Angle;
                    }
                }
            }
        }


        ConvexPatchVert=std::vector<std::vector<size_t> >(NewPatches.size(),std::vector<size_t>());
        ConcavePatchVert=std::vector<std::vector<size_t> >(NewPatches.size(),std::vector<size_t>());

        for (size_t i=0;i<VertPatchCount.size();i++)
        {
            //at least one patch per vertex
            assert(VertPatchCount[i].size()>0);
            //no more than 4
            //assert(VertPatchCount[i].size()<=4);

            //retrieve the original vert
            size_t IndexOriginal=patchMesh.vert[i].Q();

            //if was already originally convex then it will remain convex
            if (ConvexOriginal.count(IndexOriginal)>0)
            {
                //only one patch for the convex ones
                assert(patchMesh.vert[i].IsB());
                //assert(VertPatchCount[i].size()==1);
                if (VertPatchCount[i].size()!=1)
                {
                    if (printMsg) {
                        std::cout<<"WARNING: VPatch Size"<<VertPatchCount[i].size()<<std::endl;
                    }
                    for (size_t j=0;j<VertPatchCount[i].size();j++)
                    {
                        size_t IndexPatch=VertPatchCount[i][j];
                        ConvexPatchVert[IndexPatch].push_back(IndexOriginal);
                    }
                }
                else
                {
                    size_t IndexPatch=VertPatchCount[i][0];
                    ConvexPatchVert[IndexPatch].push_back(IndexOriginal);
                }
                continue;
            }
            if (ConcaveOriginal.count(IndexOriginal)>0)
            {
                assert(patchMesh.vert[i].IsB());
                //assert(VertPatchCount[i].size()<=3);
                //it remain concave
                if (VertPatchCount[i].size()==1)
                {
                    size_t IndexPatch=VertPatchCount[i][0];
                    ConcavePatchVert[IndexPatch].push_back(IndexOriginal);
                    continue;
                }
                //the one with bigger angle becomes flat, the other convex
                if (VertPatchCount[i].size()==2)
                {
                    size_t IndexPatch0=VertPatchCount[i][0];
                    size_t IndexPatch1=VertPatchCount[i][1];
                    ScalarType angle0=VertPatchAngle[i][0];
                    ScalarType angle1=VertPatchAngle[i][1];
                    if (angle0<angle1)
                        ConvexPatchVert[IndexPatch0].push_back(IndexOriginal);
                    else
                        ConvexPatchVert[IndexPatch1].push_back(IndexOriginal);
                    continue;
                }
                //all convex in this case
                if (VertPatchCount[i].size()==3)
                {
                    size_t IndexPatch0=VertPatchCount[i][0];
                    size_t IndexPatch1=VertPatchCount[i][1];
                    size_t IndexPatch2=VertPatchCount[i][2];

                    ConvexPatchVert[IndexPatch0].push_back(IndexOriginal);
                    ConvexPatchVert[IndexPatch1].push_back(IndexOriginal);
                    ConvexPatchVert[IndexPatch2].push_back(IndexOriginal);
                }
            }
            //if valence >=2 and border it becomes convex in both pathces
            if ((VertPatchCount[i].size()>=2)&&(patchMesh.vert[i].IsB()))
            {
                //                size_t IndexPatch0=VertPatchCount[i][0];
                //                size_t IndexPatch1=VertPatchCount[i][1];
                //                ConvexPatchVert[IndexPatch0].push_back(IndexOriginal);
                //                ConvexPatchVert[IndexPatch1].push_back(IndexOriginal);
                for (size_t j=0;j<VertPatchCount[i].size();j++)
                {
                    size_t IndexPatch=VertPatchCount[i][j];
                    ConvexPatchVert[IndexPatch].push_back(IndexOriginal);
                }
                continue;
            }
            //cross-intersections, all convex
            if (VertPatchCount[i].size()>=4)
            {

                //                size_t IndexPatch0=VertPatchCount[i][0];
                //                size_t IndexPatch1=VertPatchCount[i][1];
                //                size_t IndexPatch2=VertPatchCount[i][2];
                //                size_t IndexPatch3=VertPatchCount[i][3];
                //                ConvexPatchVert[IndexPatch0].push_back(IndexOriginal);
                //                ConvexPatchVert[IndexPatch1].push_back(IndexOriginal);
                //                ConvexPatchVert[IndexPatch2].push_back(IndexOriginal);
                //                ConvexPatchVert[IndexPatch3].push_back(IndexOriginal);
                for (size_t j=0;j<VertPatchCount[i].size();j++)
                {
                    size_t IndexPatch=VertPatchCount[i][j];
                    ConvexPatchVert[IndexPatch].push_back(IndexOriginal);
                }
                continue;
            }
        }
        //std::cout<<"de2"<<std::endl;

    }

    void CopyPositionsFrom(const MeshType &Source)
    {
        for (size_t i=0;i<patchMesh.vert.size();i++)
        {
            size_t IndexOriginal=patchMesh.vert[i].Q();
            patchMesh.vert[i].P()=Source.vert[IndexOriginal].cP();
        }
        updateAllMeshAttributes(patchMesh);
    }

    bool SplitByConcave(const std::set<size_t> &ConvexOriginal,
                        const std::set<size_t> &ConcaveOriginal,
                        const std::vector<size_t> &TracingVert,
                        const std::vector<std::vector<size_t> > &TracingDir,
                        std::vector<std::vector<size_t> > &LocalPatches,
                        std::vector<std::vector<size_t> > &NewPatches,
                        bool splitAll)
    //                        std::vector<std::vector<size_t> > &NewConvexPatchOriginalVert,
    //                        std::vector<std::vector<size_t> > &NewConcavePatchOriginalVert)
    {
        //SavePatchMeshConcaveQ(ConcaveOriginal);
        //SavePatchMeshConcave(ConcaveOriginal);

        TraceConcaveCandidates(ConvexOriginal,ConcaveOriginal,TracingVert,TracingDir);

        PruneConcaveCandidates(splitAll);

        if (TraceVertCandidates.size()==0)return false;
        RetrievePatchesFromPaths(LocalPatches);

        //set the index of the original mesh
        NewPatches=LocalPatches;
        for (size_t i=0;i<NewPatches.size();i++)
            for (size_t j=0;j<NewPatches[i].size();j++)
            {
                size_t IndexF=NewPatches[i][j];
                IndexF=patchMesh.face[IndexF].Q();
                NewPatches[i][j]=IndexF;
            }
        return true;
    }

    std::set<size_t> ConvexPatchMesh;

    void MapToOriginalIndex(const std::vector<std::vector<size_t> > &LocalIndex,
                            std::vector<std::vector<size_t> > &OriginalIndex)
    {
        OriginalIndex=LocalIndex;
        for (size_t i=0;i<OriginalIndex.size();i++)
            for (size_t j=0;j<OriginalIndex[i].size();j++)
            {
                size_t IndexF=OriginalIndex[i][j];
                assert(IndexF>=0);
                assert(IndexF<patchMesh.face.size());
                IndexF=patchMesh.face[IndexF].Q();
                OriginalIndex[i][j]=IndexF;
            }
    }


    bool SplitByBorders(const std::set<size_t> &ConvexOriginal,
                        const std::vector<size_t> &TracingVert,
                        const std::vector<std::vector<size_t> > &TracingDir,
                        std::vector<std::vector<size_t> > &LocalPatches,
                        std::vector<std::vector<size_t> > &NewPatches,
                        size_t MinTraces=0)
    {
        //trace all
        TraceAllBordersCandidates(ConvexOriginal,TracingVert,TracingDir);

        //then prune the colliding ones
        PruneCollidingCandidates();

        if (TraceVertCandidates.size()==0)
        {
            return false;
        }

        //remove one by one till the conditions are satisfied
        InitCandidatesPathLenghts();
        std::reverse(CandidatesPathLenghts.begin(),CandidatesPathLenghts.end());
        size_t num_removed=0;
        if (printMsg) {
            std::cout<<"*** Initial "<<CandidatesPathLenghts.size()<<" Traces"<<std::endl;
        }

        GetConvexVertices(ConvexOriginal,ConvexPatchMesh);

        //int KeptCandidates=CandidatesPathLenghts.size();
        for (size_t i=0;i<CandidatesPathLenghts.size();i++)
        {
            //if (KeptCandidates<=MinTraces)break;

            size_t CurrTrace=CandidatesPathLenghts[i].second;
            if (!CanRemove(CurrTrace,MinTraces))continue;
            num_removed++;

            //KeptCandidates--;

            TraceVertCandidates[CurrTrace].clear();
            TraceDirCandidates[CurrTrace].clear();
        }
        if (printMsg) {
            std::cout<<"*** Removed "<<num_removed<<" Traces"<<std::endl;
        }
        //exit(0);
        //erase the empty ones
        std::vector<std::vector<size_t > > TraceVertSwap;
        std::vector<std::vector<size_t > > TraceDirSwap;
        for (size_t i=0;i<TraceVertCandidates.size();i++)
        {
            if (TraceVertCandidates[i].size()==0)continue;
            TraceVertSwap.push_back(TraceVertCandidates[i]);
            TraceDirSwap.push_back(TraceDirCandidates[i]);
        }
        TraceVertCandidates=TraceVertSwap;
        TraceDirCandidates=TraceDirSwap;
        if (TraceVertCandidates.size()==0)
        {
            return false;
        }

        //finally retrieve the patches
        RetrievePatchesFromPaths(LocalPatches);
        MapToOriginalIndex(LocalPatches,NewPatches);
        return true;
    }

    PatchSplitter(MeshType &_patchMesh,
                  VertexFieldTracer<MeshType> &_VFTracer,
                  bool _printMsg):patchMesh(_patchMesh),VFTracer(_VFTracer)
    {printMsg=_printMsg;}
};

template <class MeshType,class QuadMeshType>
class PatchAssembler
{
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;

public:
    struct Parameters
    {
        bool InitialRemesh;
        bool PrintDebug;
        ScalarType EdgeSizeFactor;
        size_t PropagationSteps;
        bool Reproject;
        bool SplitAllConcave;
        bool CleanFolds;
        bool InterleaveSmooth;
        bool FinalSmooth;

        Parameters()
        {
            EdgeSizeFactor=1;
            PropagationSteps=1;
            PrintDebug=false;
            Reproject=true;
            InitialRemesh=true;
            SplitAllConcave=false;
            CleanFolds=false;
            InterleaveSmooth=false;
            FinalSmooth=true;
        }
    };

private:

    MeshType &mesh;
    QuadMeshType &quadMesh;

    //the original mesh to reproject to
    MeshType OriginalMesh;
    vcg::GridStaticPtr<FaceType,typename FaceType::ScalarType> Grid;

    Parameters Param;
    //size_t currPatch=0;

    struct PatchInfo
    {
        std::set<size_t> ConvexVertIndex;
        std::set<size_t> ConcaveVertIndex;
        std::set<size_t> FlatVertIndex;
        std::vector<size_t> IndexF;
        std::vector<size_t> TracingVert;
        std::vector<std::vector<size_t> > TracingDir;
        bool Active;
        MeshType *PatchMesh;
        VertexFieldTracer<MeshType> *PatchFieldTracer;

        PatchInfo()
        {
            PatchMesh=NULL;
            PatchFieldTracer=NULL;
        }
    };

    std::vector<PatchInfo> Patches;

    void GetBorderOrthoDirections(const MeshType &testMesh,
                                  std::vector<std::vector<CoordType> > &OrthoDirs)
    {
        OrthoDirs.clear();
        OrthoDirs.resize(testMesh.vert.size());

        for (size_t i=0;i<testMesh.face.size();i++)
            for (size_t j=0;j<testMesh.face[i].VN();j++)
            {
                if(!vcg::face::IsBorder(testMesh.face[i],j))continue;

                size_t IndexV0=vcg::tri::Index(testMesh,testMesh.face[i].cV0(j));
                size_t IndexV1=vcg::tri::Index(testMesh,testMesh.face[i].cV1(j));
                CoordType P0=testMesh.vert[IndexV0].P();
                CoordType P1=testMesh.vert[IndexV1].P();
                CoordType Dir=P1-P0;
                Dir.Normalize();

                CoordType OrthoDir=Dir;
                OrthoDir=testMesh.face[i].cN()^OrthoDir;

                vcg::Matrix33<ScalarType> Rot0=vcg::RotationMatrix(testMesh.face[i].cN(),testMesh.vert[IndexV0].N());
                vcg::Matrix33<ScalarType> Rot1=vcg::RotationMatrix(testMesh.face[i].cN(),testMesh.vert[IndexV1].N());

                CoordType OrthoDir0=Rot0*OrthoDir;
                CoordType OrthoDir1=Rot1*OrthoDir;

                OrthoDir0.Normalize();
                OrthoDir1.Normalize();


                OrthoDirs[IndexV0].push_back(OrthoDir0);
                OrthoDirs[IndexV1].push_back(OrthoDir1);
            }
    }



    void GetPatchMesh(int IndexPatch,MeshType &partition_mesh)
    {
        partition_mesh.Clear();
        assert(IndexPatch>=0);
        assert(IndexPatch<Patches.size());
        assert(Patches[IndexPatch].Active);

        vcg::tri::UpdateFlags<MeshType>::VertexClearS(mesh);
        vcg::tri::UpdateFlags<MeshType>::FaceClearS(mesh);
        for (size_t i=0;i<Patches[IndexPatch].IndexF.size();i++)
        {
            size_t IndexF=Patches[IndexPatch].IndexF[i];
            mesh.face[IndexF].SetS();
        }
        vcg::tri::UpdateSelection<MeshType>::VertexFromFaceLoose(mesh);
        vcg::tri::Append<MeshType,MeshType>::Mesh(partition_mesh,mesh,true);
        updateAllMeshAttributes(partition_mesh);
        vcg::tri::UpdateFlags<MeshType>::VertexClearS(partition_mesh);
        vcg::tri::UpdateFlags<MeshType>::FaceClearS(partition_mesh);
        vcg::tri::UpdateFlags<MeshType>::VertexClearS(mesh);
        vcg::tri::UpdateFlags<MeshType>::FaceClearS(mesh);
    }

    void SetPatchInitialDirections(int IndexPatch)
    {
        MeshType *PatchMesh=Patches[IndexPatch].PatchMesh;
        VertexFieldTracer<MeshType> *PatchFieldTracer=Patches[IndexPatch].PatchFieldTracer;

        std::vector<std::vector<CoordType> > OrthoDirs;
        GetBorderOrthoDirections(*PatchMesh,OrthoDirs);

        //then go around the borders
        //        for (size_t i=0;i<PatchMesh->face.size();i++)
        //            for (size_t j=0;j<PatchMesh->face[i].VN();j++)
        //            {
        //                if (vcg::face::IsBorder(PatchMesh->face[i],j))
        //                {
        for (size_t i=0;i<PatchMesh->vert.size();i++)
        {
            if (!PatchMesh->vert[i].IsB())continue;
            //TypeBorderVert TypeV=PatchBorderVertexType(IndexPatch,PatchMesh->face[i].V(j));
            TypeBorderVert TypeV=PatchBorderVertexType(IndexPatch,&PatchMesh->vert[i]);
            //int IndexOriginal=PatchMesh->face[i].V(j)->Q();
            size_t IndexV=i;//vcg::tri::Index(*PatchMesh,PatchMesh->face[i].V(j));
            if (TypeV==VBConvex)
            {
                //push an empty vector
                Patches[IndexPatch].TracingVert.push_back(IndexV);
                Patches[IndexPatch].TracingDir.push_back(std::vector<size_t>());
                continue;
            }
            if (TypeV==VBConcave)
            {
                Patches[IndexPatch].TracingVert.push_back(IndexV);
                CoordType Ortho0=OrthoDirs[IndexV][0];
                CoordType Ortho1=OrthoDirs[IndexV][1];
                size_t BestDir0=PatchFieldTracer->GetClosestDirTo(IndexV,Ortho0);
                size_t BestDir1=PatchFieldTracer->GetClosestDirTo(IndexV,Ortho1);
                Patches[IndexPatch].TracingDir.push_back(std::vector<size_t>());
                Patches[IndexPatch].TracingDir.back().push_back(BestDir0);
                if (BestDir1!=BestDir0)
                    Patches[IndexPatch].TracingDir.back().push_back(BestDir1);
                continue;
            }
            if (TypeV==VBFlat)
            {
                Patches[IndexPatch].TracingVert.push_back(IndexV);
                CoordType Ortho0=OrthoDirs[IndexV][0];
                CoordType Ortho1=OrthoDirs[IndexV][1];
                CoordType TargetD=Ortho0+Ortho1;
                TargetD.Normalize();
                size_t BestDir=PatchFieldTracer->GetClosestDirTo(IndexV,TargetD);
                Patches[IndexPatch].TracingDir.push_back(std::vector<size_t>(1,BestDir));
                continue;
            }
        }
        //       }
    }

    void SelectFacesWithBorderEdge()
    {
        vcg::tri::UpdateSelection<MeshType>::Clear(mesh);
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<mesh.face[i].VN();j++)
            {
                if (!vcg::face::IsBorder(mesh.face[i],j))continue;
                mesh.face[i].SetS();
            }
    }

    void SelectFacesWithBorderVertex()
    {
        vcg::tri::UpdateSelection<MeshType>::Clear(mesh);
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<mesh.face[i].VN();j++)
            {
                if (mesh.face[i].IsB(j))
                    mesh.face[i].SetS();
            }
    }

    void InvertFaceSelection()
    {
        for (size_t i=0;i<mesh.face.size();i++)
        {
            if (mesh.face[i].IsS())
                mesh.face[i].ClearS();
            else
                mesh.face[i].SetS();
        }
    }

    void InvertVertexSelection()
    {
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (mesh.vert[i].IsS())
                mesh.vert[i].ClearS();
            else
                mesh.vert[i].SetS();
        }
    }

    ScalarType AverageEdgeSize(bool OnlyBorder=true)
    {
        ScalarType AVGEdge=0;
        size_t Num=0;
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<mesh.face[i].VN();j++)
            {
                if (OnlyBorder && (!vcg::face::IsBorder(mesh.face[i],j)))continue;
                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
                CoordType P0=mesh.vert[IndexV0].P();
                CoordType P1=mesh.vert[IndexV1].P();
                AVGEdge+=(P0-P1).Norm();
                Num++;
            }
        AVGEdge/=Num;
        return (AVGEdge);
    }

    ScalarType SmallestBorderEdge()
    {
        ScalarType SmallE=std::numeric_limits<ScalarType>::max();
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<mesh.face[i].VN();j++)
            {
                if(!vcg::face::IsBorder(mesh.face[i],j))continue;
                CoordType P0=mesh.face[i].P0(j);
                CoordType P1=mesh.face[i].P1(j);
                SmallE=std::min(SmallE,(P0-P1).Norm());
            }
        return (SmallE);
    }

    void SmoothMesh()
    {
        vcg::tri::UpdateSelection<MeshType>::Clear(mesh);
        vcg::tri::UpdateSelection<MeshType>::VertexFromBorderFlag(mesh);

        InvertVertexSelection();
        vcg::PolygonalAlgorithm<MeshType>::LaplacianReproject(mesh,3,0.5,true);

        updateAllMeshAttributes(mesh);
        vcg::tri::UpdateSelection<MeshType>::Clear(mesh);
    }

    void RemeshMesh()
    {

        ScalarType EdgeStep=AverageEdgeSize()/2;

        typename vcg::tri::IsotropicRemeshing<MeshType>::Params Par;

        Par.swapFlag     = true;
        Par.collapseFlag = true;
        Par.smoothFlag=true;
        Par.projectFlag=Param.Reproject;
        Par.SetFeatureAngleDeg(100);
        Par.SetTargetLen(EdgeStep/Param.EdgeSizeFactor);
        Par.selectedOnly=true;
        Par.adapt=false;
        Par.iter=20;

        updateAllMeshAttributes(mesh);
        //vcg::tri::Smooth<MeshType>::VertexCoordLaplacian(mesh,3);
        //vcg::tri::UpdateCurvature<MeshType>::MeanAndGaussian(mesh);
//        vcg::tri::UpdateCurvature<MeshType>::PerVertex(mesh);
//        vcg::tri::FieldSmoother<MeshType>::InitByCurvature(mesh,2,true);
//        for (size_t i=0;i<mesh.vert.size();i++)
//        {
//           mesh.vert[i].Q()=i;
//        }
//        vcg::tri::UpdateColor<MeshType>::PerVertexQualityRamp(mesh);
//        for (size_t i=0;i<mesh.vert.size();i++)
//        {
//           vcg::Color4b col=mesh.vert[i].C();
//           std::cout<<"Col "<<col[0]<<" "<<col[1]<<" "<<col[2]<<std::endl;
//        }

//#ifdef SAVE_MESHES_FOR_DEBUG
//        vcg::tri::io::ExporterOBJ<MeshType>::Save(mesh,"results/colored.obj",vcg::tri::io::Mask::IOM_VERTCOLOR);
//#endif
        //exit(0);

        SelectFacesWithBorderEdge();
        //SelectFacesWithBorderVertex();
        InvertFaceSelection();

        if (Param.PrintDebug) {
            std::cout<<"Remeshing step"<<std::endl;
            std::cout<<"Edge Size " <<EdgeStep<<std::endl;
        }

        MeshType Reproject;
        vcg::tri::Append<MeshType,MeshType>::Mesh(Reproject,mesh);
        updateAllMeshAttributes(Reproject);

        vcg::tri::IsotropicRemeshing<MeshType>::Do(mesh,Reproject,Par);
        //vcg::tri::UpdateColor<MeshType>::PerFaceConstant(mesh,vcg::Color4b::LightGray);

        vcg::tri::UpdateFlags<MeshType>::Clear(mesh);
        updateAllMeshAttributes(mesh);


//#ifdef SAVE_MESHES_FOR_DEBUG
//        vcg::tri::io::ExporterOBJ<MeshType>::Save(mesh,"remeshed.obj", vcg::tri::io::Mask::IOM_NONE);
//#endif
//        exit(0);
//        std::cout<<std::flush;
    }

    void GetGeometricCorners(std::vector<CoordType> &ConvexV,
                             std::vector<CoordType> &ConcaveV)
    {
        ConvexV.clear();
        ConcaveV.clear();
        vcg::tri::UpdateSelection<MeshType>::VertexCornerBorder(mesh,M_PI-M_PI/CONVEX);
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (!mesh.vert[i].IsS())continue;
            ConvexV.push_back(mesh.vert[i].P());
        }

        //        //add the non manifold ones
        //        vcg::tri::Clean<MeshType>::CountNonManifoldVertexFF(mesh);
        //        for (size_t i=0;i<mesh.vert.size();i++)
        //        {
        //            if (!mesh.vert[i].IsS())continue;
        //            ConvexV.push_back(mesh.vert[i].P());
        //        }

        //then add the concave
        vcg::tri::UpdateSelection<MeshType>::VertexCornerBorder(mesh,M_PI+M_PI/CONVEX);
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (!mesh.vert[i].IsB())continue;
            if (mesh.vert[i].IsS())continue;
            //            typename std::vector<CoordType>::iterator It;
            //            It=std::find(ConvexV.begin(),ConvexV.end(),mesh.vert[i].P());
            //            if (It!=ConvexV.end())continue;
            ConcaveV.push_back(mesh.vert[i].P());
        }
    }



    void ColorPatch(const size_t &IndexPatch,
                    const vcg::Color4b &PatchCol)
    {
        for (size_t j=0;j<Patches[IndexPatch].IndexF.size();j++)
        {
            int CurrF=Patches[IndexPatch].IndexF[j];
            mesh.face[CurrF].C()=PatchCol;
        }
    }

    void ColorPatches()
    {
        for (size_t i=0;i<Patches.size();i++)
        {
            if (!Patches[i].Active)continue;
            vcg::Color4b PatchCol=vcg::Color4b::Scatter(Patches.size(),i);
            ColorPatch(i,PatchCol);
        }
    }

    enum TypeBorderVert{VBConcave,VBConvex,VBFlat};

    TypeBorderVert PatchBorderVertexType(const size_t &IndexPatch,
                                         const size_t &IndexOriginal)
    {
        if (Patches[IndexPatch].FlatVertIndex.count(IndexOriginal)>0)return VBFlat;
        if (Patches[IndexPatch].ConcaveVertIndex.count(IndexOriginal)>0)return VBConcave;
        if (Patches[IndexPatch].ConvexVertIndex.count(IndexOriginal)>0)return VBConvex;
        assert(0);
    }

    TypeBorderVert PatchBorderVertexType(const size_t &IndexPatch,
                                         const VertexType *v)
    {
        size_t IndexOriginal=v->Q();
        return  PatchBorderVertexType(IndexPatch,IndexOriginal);
    }

#ifdef DRAWTRACE

    void GLDrawPatch(const size_t &IndexPatch)
    {
        MeshType *PatchMesh=Patches[IndexPatch].PatchMesh;

        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.99998);
        glPointSize(10);

        glBegin(GL_POINTS);
        for (size_t i=0;i<PatchMesh->face.size();i++)
        {
            int CurrF=i;//Patches[IndexPatch].IndexF[i];
            CoordType Bary=(PatchMesh->face[CurrF].P(0)+
                            PatchMesh->face[CurrF].P(1)+
                            PatchMesh->face[CurrF].P(2))/3;
            for (size_t j=0;j<3;j++)
            {
                vcg::Color4b Col;
                VertexType *currV=PatchMesh->face[CurrF].V(j);
                if (!currV->IsB())continue;
                TypeBorderVert TypeV=PatchBorderVertexType(IndexPatch,currV);
                if (TypeV==VBConvex)
                {
                    Col=vcg::Color4b::Red;
                }
                if (TypeV==VBConcave)
                {
                    Col=vcg::Color4b::Blue;
                }
                if (TypeV==VBFlat)
                {
                    Col=vcg::Color4b::Green;
                }
                CoordType Pos=currV->P()*0.8+Bary*0.2;
                vcg::glColor(Col);
                vcg::glVertex(Pos);
                //std::cout<<"Test 1"<<std::endl;
            }
        }
        glEnd();
        glPopAttrib();
    }

    void GLDrawPatchDir(const size_t &IndexPatch)
    {
        assert(IndexPatch>=0);
        assert(IndexPatch<Patches.size());
        assert(Patches[IndexPatch].Active);

        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.998);

        MeshType *PatchMesh=Patches[IndexPatch].PatchMesh;
        VertexFieldTracer<MeshType> *PatchFieldTracer=Patches[IndexPatch].PatchFieldTracer;

        for (size_t i=0;i<Patches[IndexPatch].TracingVert.size();++i)
        {
            size_t IndexV=Patches[IndexPatch].TracingVert[i];
            assert(IndexV>=0);
            assert(IndexV<PatchMesh->vert.size());

            CoordType Pos0=PatchMesh->vert[IndexV].P();

            vcg::glColor(vcg::Color4b(0,255,0,255));
            glLineWidth(5);
            if (Patches[IndexPatch].TracingDir[i].size()>1)
            {
                vcg::glColor(vcg::Color4b(255,0,0,255));
                glLineWidth(10);
            }

            for (size_t j=0;j<Patches[IndexPatch].TracingDir[i].size();++j)
            {
                size_t currDirection=Patches[IndexPatch].TracingDir[i][j];
                CoordType Dir=PatchFieldTracer->GetDirection(IndexV,currDirection);

                CoordType Pos1=Pos0+Dir*mesh.bbox.Diag()*0.01;

                glBegin(GL_LINES);
                vcg::glVertex(Pos0);
                vcg::glVertex(Pos1);
                glEnd();

            }
        }
        glPopAttrib();
    }
#endif

    void CopyPatchCornerByOriginalIndex(const size_t IndexPatch,
                                        const std::vector<size_t> &ConvexOrigIndex,
                                        const std::vector<size_t> &ConcaveOrigIndex)
    {
        std::set<size_t> ConvexVSet(ConvexOrigIndex.begin(),ConvexOrigIndex.end());
        std::set<size_t> ConcaveVSet(ConcaveOrigIndex.begin(),ConcaveOrigIndex.end());

        Patches[IndexPatch].ConvexVertIndex.clear();
        Patches[IndexPatch].ConcaveVertIndex.clear();
        Patches[IndexPatch].FlatVertIndex.clear();

        for (size_t i=0;i<Patches[IndexPatch].PatchMesh->vert.size();i++)
        {
            size_t IndexOriginal=Patches[IndexPatch].PatchMesh->vert[i].Q();

            if (ConvexVSet.count(IndexOriginal)>0)
            {
                Patches[IndexPatch].ConvexVertIndex.insert(IndexOriginal);
                continue;
            }
            if (ConcaveVSet.count(IndexOriginal)>0)
            {
                Patches[IndexPatch].ConcaveVertIndex.insert(IndexOriginal);
                continue;
            }
            if (Patches[IndexPatch].PatchMesh->vert[i].IsB())
            {
                Patches[IndexPatch].FlatVertIndex.insert(IndexOriginal);
                continue;
            }
        }
    }

    void CopyPatchCornerByPos(size_t IndexPatch,
                              const std::vector<CoordType> &ConvexV,
                              const std::vector<CoordType> &ConcaveV)
    {
        std::set<CoordType> ConvexVSet(ConvexV.begin(),ConvexV.end());
        std::set<CoordType> ConcaveVSet(ConcaveV.begin(),ConcaveV.end());

        Patches[IndexPatch].ConvexVertIndex.clear();
        Patches[IndexPatch].ConcaveVertIndex.clear();
        Patches[IndexPatch].FlatVertIndex.clear();

        for (size_t i=0;i<Patches[IndexPatch].PatchMesh->vert.size();i++)
        {
            CoordType Pos=Patches[IndexPatch].PatchMesh->vert[i].P();
            size_t IndexOriginal=Patches[IndexPatch].PatchMesh->vert[i].Q();

            if (ConvexVSet.count(Pos)>0)
            {
                Patches[IndexPatch].ConvexVertIndex.insert(IndexOriginal);
                continue;
            }
            if (ConcaveVSet.count(Pos)>0)
            {
                Patches[IndexPatch].ConcaveVertIndex.insert(IndexOriginal);
                continue;
            }
            if (Patches[IndexPatch].PatchMesh->vert[i].IsB())
            {
                Patches[IndexPatch].FlatVertIndex.insert(IndexOriginal);
                continue;
            }
        }
    }

    void InitPatchMesh(size_t IndexPatch)
    {
        assert(IndexPatch>=0);
        assert(IndexPatch<Patches.size());
        assert(Patches[IndexPatch].Active);

        //create the temporaty mesh
        Patches[IndexPatch].PatchMesh=new MeshType();
        GetPatchMesh(IndexPatch,*Patches[IndexPatch].PatchMesh);
        //update the field
        vcg::tri::CrossField<MeshType>::SetVertCrossVectorFromFace(*Patches[IndexPatch].PatchMesh);
        //create the tracer
        Patches[IndexPatch].PatchFieldTracer=new VertexFieldTracer<MeshType>(*Patches[IndexPatch].PatchMesh,
                                                                             Param.PropagationSteps);
    }

    void SetInitialPatch(const std::vector<CoordType> &ConvexV,
                         const std::vector<CoordType> &ConcaveV)
    {
        Patches.clear();

        Patches.push_back(PatchInfo());

        for (size_t i=0;i<mesh.face.size();i++)
            Patches.back().IndexF.push_back(i);

        Patches.back().Active=true;

        InitPatchMesh(0);

        //copy corners
        CopyPatchCornerByPos(0,ConvexV,ConcaveV);

        //then initialize the tracing directions
        SetPatchInitialDirections(0);
    }

    void InitField()
    {
        if (Param.PrintDebug) {
            std::cout<<"Initializig Field"<<std::endl;
        }
        typename FieldSmoother<MeshType>::SmoothParam SParam;

        SParam.align_borders=true;
        SParam.curvRing=2;
        SParam.alpha_curv=0;
        SParam.Ndir=4;
        SParam.sharp_thr=0;
        SParam.curv_thr=0;
        SParam.SmoothM=SMNPoly;

#ifdef SAVE_MESHES_FOR_DEBUG
        //vcg::tri::io::ExporterOBJ<MeshType>::Save(mesh,"results/tracer_field.obj", vcg::tri::io::Mask::IOM_NONE);
#endif

        FieldSmoother<MeshType>::SmoothDirections(mesh,SParam);
        if (Param.PrintDebug) {
            std::cout<<"Done! Now propagating on Vertices"<<std::endl;
        }

        vcg::tri::CrossField<MeshType>::SetVertCrossVectorFromFace(mesh);

        if (Param.PrintDebug)
            std::cout<<"Done!"<<std::endl;

        //InitIndexOnQ();
    }

    void CheckIndexOnQ(int TestIndex)
    {        
        if (Param.PrintDebug) {
            std::cout<<"*** TESTING INDEX "<<TestIndex<<std::endl;
        }
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            size_t Index0=i;
            size_t Index1=mesh.vert[i].Q();
            if (Index0!=Index1)
            {
                std::cout<<"ERROR: i: "<<i<<" Q: "<<mesh.vert[i].Q()<<std::endl;
                assert(0);
            }
        }

        for (size_t i=0;i<mesh.face.size();i++)
        {
            size_t Index0=i;
            size_t Index1=mesh.face[i].Q();
            if (Index0!=Index1)
            {
                std::cout<<"ERROR i: "<<i<<" Q: "<<mesh.face[i].Q()<<std::endl;
                assert(0);
            }
        }
    }

    void CheckPatches()
    {
        //CheckIndexOnQ(10);

        for (size_t i=0;i<Patches.size();i++)
        {
            if (!Patches[i].Active)continue;
            //std::cout<<"Testing Patch "<<i<<std::endl;
            assert(Patches[i].PatchMesh!=NULL);
            for (size_t j=0;j<Patches[i].PatchMesh->vert.size();j++)
            {
                VertexType *currV=&Patches[i].PatchMesh->vert[j];
                if (!currV->IsB())continue;
                TypeBorderVert TypeV=PatchBorderVertexType(i,currV);
            }
            for (size_t j=0;j<Patches[i].IndexF.size();j++)
            {
                int CurrF=Patches[i].IndexF[j];
                for (size_t k=0;k<3;k++)
                {
                    VertexType *currV=mesh.face[CurrF].V(k);
                    size_t Index0=vcg::tri::Index(mesh,currV);
                    size_t Index1=currV->Q();
                    assert(Index0==Index1);
                    if (!currV->IsB())continue;
                    TypeBorderVert TypeV=PatchBorderVertexType(i,currV);
                }
            }
            //std::cout<<"Done Testing"<<std::endl;
        }
    }

    void SplitPatchIntoSubPatches(size_t IndexPatch,
                                  const std::vector<std::vector<size_t> > &NewPatches,
                                  std::vector<size_t> &NewIndexes)
    {
        if(Param.PrintDebug)
            std::cout<<"De-activated patch: "<<IndexPatch<<std::endl;

        Patches[IndexPatch].Active=false;
        NewIndexes.clear();
        //std::cout<<"Re-initializing Num Patches: "<<NewPatches.size()<<std::endl;
        for (size_t i=0;i<NewPatches.size();i++)
        {
            //std::cout<<"Re-initializing patch: "<<i<<std::endl;
            Patches.push_back(PatchInfo());
            Patches.back().IndexF=NewPatches[i];
            Patches.back().Active=true;
            NewIndexes.push_back(Patches.size()-1);

            if(Param.PrintDebug)
                std::cout<<"Created new patch: "<<NewIndexes.back()<<std::endl;

        }
    }

    void UpdateSubPatchesVertexType(const std::vector<size_t> &ToUpdate,
                                    const std::vector<std::vector<size_t> > &NewConvexPatchOriginalVert,
                                    const std::vector<std::vector<size_t> > &NewConcavePatchOriginalVert)
    {
        assert(ToUpdate.size()==NewConvexPatchOriginalVert.size());
        assert(ToUpdate.size()==NewConcavePatchOriginalVert.size());

        if (Param.PrintDebug)
            std::cout<<"Re-initializing Num Patches: "<<ToUpdate.size()<<std::endl;
        for (size_t i=0;i<ToUpdate.size();i++)
        {
            size_t IndexPatch=ToUpdate[i];
            InitPatchMesh(IndexPatch);
            if (Param.PrintDebug)
                std::cout<<"updatig patch vertex type: "<<IndexPatch<<std::endl;

            CopyPatchCornerByOriginalIndex(IndexPatch,
                                           NewConvexPatchOriginalVert[i],
                                           NewConcavePatchOriginalVert[i]);
            if (Param.PrintDebug)
                std::cout<<"setting new tracing directions "<<IndexPatch<<std::endl;

            SetPatchInitialDirections(IndexPatch);
        }
        if (Param.PrintDebug)
            std::cout<<"done "<<std::endl;

    }

    std::vector<int> FacePatches;
    void UpdateFacePatches()
    {
        FacePatches.clear();
        FacePatches.resize(mesh.face.size(),-1);

        for (size_t i=0;i<Patches.size();i++)
        {
            if (!Patches[i].Active)continue;
            for (size_t j=0;j<Patches[i].IndexF.size();j++)
            {
                size_t IndexF=Patches[i].IndexF[j];
                assert(IndexF>=0);
                assert(IndexF<mesh.face.size());
                FacePatches[IndexF]=i;
            }
        }
    }

    bool OnBorderPatch(const FaceType &f,
                       size_t IndexE)
    {
        if (vcg::face::IsBorder(f,IndexE))return false;
        FaceType *f_opp=f.cFFp(IndexE);
        size_t Index0=vcg::tri::Index(mesh,f);
        size_t Index1=vcg::tri::Index(mesh,f_opp);
        int Partition0=FacePatches[Index0];
        int Partition1=FacePatches[Index1];
        return (Partition0!=Partition1);
    }

    std::set<std::pair<size_t,size_t> > BorderPatches;
    void UpdateBorderPatches()
    {
        BorderPatches.clear();
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<mesh.face[i].VN();j++)
            {
                if (!OnBorderPatch(mesh.face[i],j))continue;
                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
                BorderPatches.insert(std::pair<size_t,size_t>(std::min(IndexV0,IndexV1),
                                                              std::max(IndexV0,IndexV1)));
            }
    }

    //    void UpdateVertexPatchCount()
    //    {

    //    }

    void ReProjectMesh()
    {
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (mesh.vert[i].IsB())continue;
            if (FixedVert.count(i)==0)continue;
            ScalarType maxD=mesh.bbox.Diag();
            ScalarType minD=0;
            CoordType closestPT;
            FaceType *f=vcg::tri::GetClosestFaceBase(OriginalMesh,Grid,mesh.vert[i].P(),maxD,minD,closestPT);
            mesh.vert[i].P()=closestPT;
        }
    }


    void SmoothPatchStep(ScalarType Damp)
    {
        //save old positio
        updateAllMeshAttributes(mesh);

        //        std::vector<CoordType> OldNorm;
        //        for (size_t i=0;i<mesh.face.size();i++)
        //            OldNorm.push_back(mesh.face[i].N());

        std::vector<CoordType> AvPos(mesh.vert.size(),CoordType(0,0,0));
        std::vector<size_t> NumDiv(mesh.vert.size(),0);

        vcg::tri::UpdateFlags<MeshType>::VertexClearS(mesh);
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<mesh.face[i].VN();j++)
            {
                CoordType Pos0=mesh.face[i].P0(j);
                CoordType Pos1=mesh.face[i].P1(j);

                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));

                std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),
                                             std::max(IndexV0,IndexV1));

                if (vcg::face::IsBorder(mesh.face[i],j))continue;

                if (BorderPatches.count(key)==0)continue;

                mesh.vert[IndexV0].SetS();
                mesh.vert[IndexV1].SetS();
                AvPos[IndexV0]+=Pos1;
                AvPos[IndexV1]+=Pos0;
                NumDiv[IndexV0]++;
                NumDiv[IndexV1]++;
            }
        for (size_t i=0;i<AvPos.size();i++)
        {
            if (mesh.vert[i].IsB()){
                mesh.vert[i].SetS();
                continue;
            }
            if (NumDiv[i]==0)continue;
            AvPos[i]/=NumDiv[i];

            if (FixedVert.count(i)==0)
                mesh.vert[i].P()=mesh.vert[i].P()*Damp+AvPos[i]*(1-Damp);
            else
                mesh.vert[i].SetS();
        }


        AvPos=std::vector<CoordType>(mesh.vert.size(),CoordType(0,0,0));
        NumDiv=std::vector<size_t>(mesh.vert.size(),0);
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<mesh.face[i].VN();j++)
            {
                CoordType Pos0=mesh.face[i].P0(j);
                CoordType Pos1=mesh.face[i].P1(j);

                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));

                AvPos[IndexV0]+=Pos1;
                AvPos[IndexV1]+=Pos0;
                NumDiv[IndexV0]++;
                NumDiv[IndexV1]++;
            }

        for (size_t i=0;i<AvPos.size();i++)
        {
            if (mesh.vert[i].IsB())continue;
            if (mesh.vert[i].IsS())continue;
            if (NumDiv[i]==0)continue;
            AvPos[i]/=NumDiv[i];

            mesh.vert[i].P()=mesh.vert[i].P()*Damp+AvPos[i]*(1-Damp);
        }

        //vcg::tri::UpdateSelection<MeshType>::VertexInvert(mesh);
        //vcg::tri::Smooth<MeshType>::VertexCoordLaplacian(mesh,1,true);

        vcg::tri::UpdateFlags<MeshType>::VertexClearS(mesh);

        //MakeFieldTangential(OldNorm);
    }

    std::set<size_t> FixedVert;

    void CleanFoldsNeeded()
    {
        vcg::tri::UpdateSelection<MeshType>::Clear(mesh);

        vcg::tri::Clean<MeshType>::SelectFoldedFaceFromOneRingFaces(mesh,-0.58);
        size_t NumSel=vcg::tri::UpdateSelection<MeshType>::FaceCount(mesh);
        if (NumSel==0)return;

        vcg::tri::UpdateSelection<MeshType>::VertexFromFaceLoose(mesh);

        //check that are not border
        for (size_t i=0;i<mesh.vert.size();i++)
            if (mesh.vert[i].IsB())mesh.vert[i].ClearS();

        vcg::tri::Smooth<MeshType>::VertexCoordLaplacian(mesh,1,true);

        vcg::tri::UpdateSelection<MeshType>::Clear(mesh);

        updateAllMeshAttributes<MeshType>(mesh);
    }

    void SmoothPatches(ScalarType Damp,size_t numSteps,bool Reproject)
    {
        for (size_t i=0;i<numSteps;i++)
        {
            SmoothPatchStep(Damp);
            if (Reproject)
                ReProjectMesh();
        }

        vcg::tri::UpdateFlags<MeshType>::Clear(mesh);
        updateAllMeshAttributes(mesh);
    }

    void UpdatePatchPosition(size_t IndexPatch)
    {
        //std::cout<<"Patch "<<IndexPatch<<std::endl;

        assert(IndexPatch>=0);
        assert(IndexPatch<Patches.size());
        assert(Patches[IndexPatch].Active);
        MeshType *PatchMesh=Patches[IndexPatch].PatchMesh;
        for (size_t i=0;i<PatchMesh->vert.size();i++)
        {
            //std::cout<<"Vertex "<<i<<std::endl;
            size_t IndexOriginal=PatchMesh->vert[i].Q();
            assert(IndexOriginal>=0);
            assert(IndexOriginal<mesh.vert.size());
            PatchMesh->vert[i].P()=mesh.vert[IndexOriginal].cP();
        }
        updateAllMeshAttributes(*PatchMesh);
    }

    void UpdatePatchPosition()
    {
        for (size_t i=0;i<Patches.size();i++)
        {
            if(!Patches[i].Active)continue;
            UpdatePatchPosition(i);
        }
    }

    bool SplitConcave(size_t IndexPatch)
    {

        assert(IndexPatch>=0);
        assert(IndexPatch<Patches.size());
        assert(Patches[IndexPatch].Active);
        MeshType *PatchMesh=Patches[IndexPatch].PatchMesh;
        VertexFieldTracer<MeshType> *PatchFieldTracer=Patches[IndexPatch].PatchFieldTracer;

        PatchSplitter<MeshType> PSplit(*PatchMesh,*PatchFieldTracer,Param.PrintDebug);

        std::vector<std::vector<size_t> > NewPatches,LocalPatches;
        std::vector<std::vector<size_t> > NewConvexPatchOriginalVert;
        std::vector<std::vector<size_t> > NewConcavePatchOriginalVert;

        bool Splitted=PSplit.SplitByConcave(Patches[IndexPatch].ConvexVertIndex,
                                            Patches[IndexPatch].ConcaveVertIndex,
                                            Patches[IndexPatch].TracingVert,
                                            Patches[IndexPatch].TracingDir,
                                            LocalPatches,NewPatches,
                                            Param.SplitAllConcave);//NewConvexPatchOriginalVert,NewConcavePatchOriginalVert);
        if (NewPatches.size()==1)return false;
        if (!Splitted)return false;


        if (Param.PrintDebug) {
            std::cout<<"Convave Tracing split into patches:"<<NewPatches.size()<<std::endl;
        }

        //update patches structures
        std::vector<size_t> NewIndexes;
        SplitPatchIntoSubPatches(IndexPatch,NewPatches,NewIndexes);

        UpdateFacePatches();
        UpdateBorderPatches();

        std::vector<std::vector<std::vector<size_t> > > NewPatchIndexes;
        NewPatchIndexes.push_back(NewPatches);
        UpdateSmoothFixedVert(NewPatchIndexes);

        //smooth to get which ones are concave or convex
        std::vector<CoordType> OldPos;
        for (size_t i=0;i<mesh.vert.size();i++)
            OldPos.push_back(mesh.vert[i].P());

        if (!Param.InterleaveSmooth)
            SmoothPatches(0.5,5,false);//do a quick smooth no reprojection
        else
            SmoothPatches(0.9,20,Param.Reproject);//do a quick smooth no reprojection

        PSplit.CopyPositionsFrom(mesh);
        PSplit.GetVertexTopology(LocalPatches,
                                 Patches[IndexPatch].ConvexVertIndex,
                                 Patches[IndexPatch].ConcaveVertIndex,
                                 NewConvexPatchOriginalVert,
                                 NewConcavePatchOriginalVert);

        //and recompute the field if the change of smooth is permanent
        if (Param.InterleaveSmooth)
            InitField();
        else
        {
            //set back old positions
            for (size_t i=0;i<mesh.vert.size();i++)
                mesh.vert[i].P()=OldPos[i];

            updateAllMeshAttributes(mesh);
        }
        InitIndexOnQ();

        UpdateSubPatchesVertexType(NewIndexes,NewConvexPatchOriginalVert,NewConcavePatchOriginalVert);

        //check consistency
        CheckPatches();

        //volor the patches
        ColorPatches();

        UpdatePatchPosition();

        return true;
    }

    void UpdateSmoothFixedVert(const std::vector<std::vector<std::vector<size_t> > > &NewPatches)
    {
        //for each old patch
        vcg::tri::UpdateSelection<MeshType>::FaceClear(mesh);
        for (size_t i=0;i<NewPatches.size();i++)
            //for each new patch
            for (size_t j=0;j<NewPatches[i].size();j++)
                //for each face
                for (size_t k=0;k<NewPatches[i][j].size();k++)
                {
                    size_t IndexF=NewPatches[i][j][k];
                    mesh.face[IndexF].SetS();
                }

        FixedVert.clear();
        vcg::tri::UpdateSelection<MeshType>::VertexFromFaceStrict(mesh);
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (mesh.vert[i].IsB()){FixedVert.insert(i);continue;}
            if (!mesh.vert[i].IsS()){FixedVert.insert(i);continue;}
        }
        vcg::tri::UpdateSelection<MeshType>::Clear(mesh);
    }

    bool SplitNonOk()
    {

        std::vector<size_t> To_Split;

        //size_t InitialPatchNum=Patches.size();
        for (size_t i=0;i<Patches.size();i++)
        {
            if (!Patches[i].Active)continue;
            //count the number of corners
            if (Patches[i].ConvexVertIndex.size()<3)
            {
                To_Split.push_back(i);
                continue;
            }
            if (Patches[i].ConvexVertIndex.size()>6)
            {
                To_Split.push_back(i);
                continue;
            }
            MeshType *CurrMesh=Patches[i].PatchMesh;
            size_t holes=vcg::tri::Clean<MeshType>::CountHoles(*CurrMesh);
            if (holes!=1)
            {
                To_Split.push_back(i);
                continue;
            }
            size_t genus=vcg::tri::Clean<MeshType>::MeshGenus(*CurrMesh);
            if (genus!=0)
            {
                To_Split.push_back(i);
                continue;
            }
            size_t nonManifV=vcg::tri::Clean<MeshType>::CountNonManifoldVertexFF(*CurrMesh);
            if (nonManifV!=0)
            {                
                if (Param.PrintDebug) {
                    std::cout<<"TO SPLIT"<<std::endl;
                }
                To_Split.push_back(i);
                continue;
            }
        }

        if (Param.PrintDebug) {
            std::cout<<"Patches Needed to be Resplitted "<<To_Split.size()<<std::endl;
        }

        std::vector<PatchSplitter<MeshType> > PSplit;
        for (size_t i=0;i<To_Split.size();i++)
        {
            size_t IndexPatch=To_Split[i];
            MeshType *PatchMesh=Patches[IndexPatch].PatchMesh;
            VertexFieldTracer<MeshType> *PatchFieldTracer=Patches[IndexPatch].PatchFieldTracer;
            PSplit.push_back(PatchSplitter<MeshType>(*PatchMesh,*PatchFieldTracer,Param.PrintDebug));
        }

        std::vector<std::vector<size_t> > NewIndexes(To_Split.size());
        std::vector<std::vector<std::vector<size_t> > > NewPatches(To_Split.size());
        std::vector<std::vector<std::vector<size_t> > >LocalPatches(To_Split.size());

        bool splitted=false;
        std::vector<bool> HadSplitted(To_Split.size(),false);
        for (size_t i=0;i<To_Split.size();i++)
        {
            if (Param.PrintDebug) {
                std::cout<<"Splitting Patch "<<To_Split[i]<<std::endl;
            }
            size_t IndexPatch=To_Split[i];

            HadSplitted[i]=PSplit[i].SplitByBorders(Patches[IndexPatch].ConvexVertIndex,
                                                    Patches[IndexPatch].TracingVert,
                                                    Patches[IndexPatch].TracingDir,
                                                    LocalPatches[i],NewPatches[i]);


            if (Param.PrintDebug) {
                std::cout<<"Border Tracing split into patches:"<<NewPatches[i].size()<<std::endl;
            }

            //update patches structures;
            splitted|=HadSplitted[i];

            if (HadSplitted[i])
            {
                SplitPatchIntoSubPatches(IndexPatch,NewPatches[i],NewIndexes[i]);
                assert(NewPatches[i].size()>0);
            }
        }
        if (!splitted)return false;

        //then smooth
        UpdateFacePatches();
        UpdateBorderPatches();

        if (Param.InterleaveSmooth)
        {
            UpdateSmoothFixedVert(NewPatches);
            SmoothPatches(0.9,20,Param.Reproject);
        }

        //get new vertex topology
        std::vector<std::vector<std::vector<size_t> > > NewConvexPatchOriginalVert(PSplit.size());
        std::vector<std::vector<std::vector<size_t> > > NewConcavePatchOriginalVert(PSplit.size());

        for (size_t i=0;i<PSplit.size();i++)
        {
            if (Param.PrintDebug)
                std::cout<<"Get Topology Step"<<std::endl;

            size_t IndexPatch=To_Split[i];
            PSplit[i].CopyPositionsFrom(mesh);
            if (HadSplitted[i])
                PSplit[i].GetVertexTopology(LocalPatches[i],
                                            Patches[IndexPatch].ConvexVertIndex,
                                            Patches[IndexPatch].ConcaveVertIndex,
                                            NewConvexPatchOriginalVert[i],
                                            NewConcavePatchOriginalVert[i]);
        }

        //and recompute the field if had changed
        if (Param.InterleaveSmooth)
            InitField();

        InitIndexOnQ();

        //create the new submeshes and update vertex type
        for (size_t i=0;i<PSplit.size();i++)
        {
            if (Param.PrintDebug)
                std::cout<<"Update Topology Step"<<std::endl;

            if (HadSplitted[i])
                UpdateSubPatchesVertexType(NewIndexes[i],
                                           NewConvexPatchOriginalVert[i],
                                           NewConcavePatchOriginalVert[i]);

            if (Param.PrintDebug)
                std::cout<<"Done"<<std::endl;

        }
        if (Param.PrintDebug)
            std::cout<<"check consistency "<<std::endl;

        CheckPatches();

        //volor the patches
        ColorPatches();

        UpdatePatchPosition();

        return true;
    }

    //    bool SplitSingleFacePatch(size_t IndexPatch)
    //    {

    //    }

    bool ForceSplit(size_t IndexPatch)
    {
        assert(IndexPatch>=0);
        assert(IndexPatch<Patches.size());
        assert(Patches[IndexPatch].Active);
        MeshType *PatchMesh=Patches[IndexPatch].PatchMesh;
        VertexFieldTracer<MeshType> *PatchFieldTracer=Patches[IndexPatch].PatchFieldTracer;

        PatchSplitter<MeshType> PSplit(*PatchMesh,*PatchFieldTracer,Param.PrintDebug);


        std::vector<std::vector<size_t> > LocalPatches,NewPatches;
        bool splitted=PSplit.SplitByBorders(Patches[IndexPatch].ConvexVertIndex,
                                            Patches[IndexPatch].TracingVert,
                                            Patches[IndexPatch].TracingDir,
                                            LocalPatches,NewPatches,4);

        std::vector<size_t> NewIndexes;
        SplitPatchIntoSubPatches(IndexPatch,NewPatches,NewIndexes);
        assert(NewPatches.size()>0);

        if (!splitted)return false;

        //then smooth
        UpdateFacePatches();
        UpdateBorderPatches();

        //        if (Param.InterleaveSmooth)
        //        {
        //            UpdateSmoothFixedVert(NewPatches);
        //            SmoothPatches(0.9,20,Param.Reproject);
        //        }

        //get new vertex topology
        std::vector<std::vector<size_t> > NewConvexPatchOriginalVert;
        std::vector<std::vector<size_t> > NewConcavePatchOriginalVert;

        if (Param.PrintDebug)
            std::cout<<"Get Topology Step"<<std::endl;

        PSplit.CopyPositionsFrom(mesh);
        PSplit.GetVertexTopology(LocalPatches,
                                 Patches[IndexPatch].ConvexVertIndex,
                                 Patches[IndexPatch].ConcaveVertIndex,
                                 NewConvexPatchOriginalVert,
                                 NewConcavePatchOriginalVert);
        //        }

        //        //and recompute the field if had changed
        //        if (Param.InterleaveSmooth)
        //            InitField();

        InitIndexOnQ();

        if (Param.PrintDebug)
            std::cout<<"Update Topology Step"<<std::endl;

        UpdateSubPatchesVertexType(NewIndexes,
                                   NewConvexPatchOriginalVert,
                                   NewConcavePatchOriginalVert);

        if (Param.PrintDebug)
            std::cout<<"Done"<<std::endl;


        if (Param.PrintDebug)
            std::cout<<"check consistency "<<std::endl;

        CheckPatches();

        //volor the patches
        ColorPatches();

        UpdatePatchPosition();

        return true;
    }

    void InitIndexOnQ()
    {
        //save vertex index on quality
        for (size_t i=0;i<mesh.vert.size();i++)
            mesh.vert[i].Q()=i;

        //save face index on quality
        for (size_t i=0;i<mesh.face.size();i++)
            mesh.face[i].Q()=i;
    }

    void RemoveNotTopologicallyOKPartitions()
    {
        std::vector<bool> IsOkPatch(Patches.size(),true);
        bool NeedErase=false;
        for (size_t i=0;i<Patches.size();i++)
        {
            if (!Patches[i].Active)continue;
            //check the holes
            MeshType *CurrMesh=Patches[i].PatchMesh;
            size_t holes=vcg::tri::Clean<MeshType>::CountHoles(*CurrMesh);
            if (holes!=1)
            {
                if (Param.PrintDebug) {
                    std::cout<<"WARNING - Mesh non Disk-Like - ERASED!"<<std::endl;
                }
                IsOkPatch[i]=false;
                NeedErase=true;
            }
            size_t nonManifV=vcg::tri::Clean<MeshType>::CountNonManifoldVertexFF(*CurrMesh);
            if (nonManifV!=0)
            {
                if (Param.PrintDebug) {
                    std::cout<<"WARNING - Non Manifold V partition - ERASED!"<<std::endl;
                }
                IsOkPatch[i]=false;
                NeedErase=true;
            }
        }

        if (!NeedErase)return;

        vcg::tri::UpdateSelection<MeshType>::Clear(mesh);
        for (size_t i=0;i<Patches.size();i++)
        {
            if (!Patches[i].Active)continue;
            if (!IsOkPatch[i])continue;
            for(size_t j=0;j<Patches[i].IndexF.size();j++)
                mesh.face[Patches[i].IndexF[j]].SetS();
        }
        vcg::tri::UpdateSelection<MeshType>::VertexFromFaceLoose(mesh);

        //create the new mesh
        MeshType NewMesh;
        vcg::tri::Append<MeshType,MeshType>::Mesh(NewMesh,mesh,true);
        mesh.Clear();
        vcg::tri::Append<MeshType,MeshType>::Mesh(mesh,NewMesh);
        InitIndexOnQ();
        updateAllMeshAttributes(mesh);

        //then create the vertex map
        std::map<size_t,size_t> VertexMap;
        for (size_t i=0;i<NewMesh.vert.size();i++)
        {
            size_t OldIndx=NewMesh.vert[i].Q();
            assert(VertexMap.count(OldIndx)==0);
            VertexMap[OldIndx]=i;
        }
        //then the face map
        std::map<size_t,size_t> FaceMap;
        for (size_t i=0;i<NewMesh.face.size();i++)
        {
            size_t OldIndx=NewMesh.face[i].Q();
            assert(FaceMap.count(OldIndx)==0);
            FaceMap[OldIndx]=i;
        }

        //then remove the patches
        std::vector<PatchInfo> NewPatches;
        for (size_t i=0;i<Patches.size();i++)
        {
            if (!Patches[i].Active)continue;
            if (!IsOkPatch[i])
            {
                delete(Patches[i].PatchMesh);
                delete(Patches[i].PatchFieldTracer);
                continue;
            }
            NewPatches.push_back(Patches[i]);
        }
        Patches=NewPatches;

        //remap indexes
        for (size_t i=0;i<Patches.size();i++)
        {
            if (!Patches[i].Active)continue;
            for(size_t j=0;j<Patches[i].IndexF.size();j++)
            {
                size_t IndexOld=Patches[i].IndexF[j];
                assert(FaceMap.count(IndexOld)>0);
                size_t IndexNew=FaceMap[IndexOld];
                assert(IndexNew>=0);
                assert(IndexNew<mesh.face.size());
                Patches[i].IndexF[j]=IndexNew;
            }

            std::vector<size_t> OldConvexIndex(Patches[i].ConvexVertIndex.begin(),
                                               Patches[i].ConvexVertIndex.end());
            std::vector<size_t> OldConcaveIndex(Patches[i].ConcaveVertIndex.begin(),
                                                Patches[i].ConcaveVertIndex.end());
            std::vector<size_t> OldFlatIndex(Patches[i].FlatVertIndex.begin(),
                                             Patches[i].FlatVertIndex.end());

            Patches[i].ConvexVertIndex.clear();
            Patches[i].ConcaveVertIndex.clear();
            Patches[i].FlatVertIndex.clear();

            for(size_t j=0;j<OldConvexIndex.size();j++)
            {
                size_t IndexOld=OldConvexIndex[j];
                assert(VertexMap.count(IndexOld)>0);
                size_t IndexNew=VertexMap[IndexOld];
                assert(IndexNew>=0);
                assert(IndexNew<mesh.vert.size());
                Patches[i].ConvexVertIndex.insert(IndexNew);
            }

            for(size_t j=0;j<OldConcaveIndex.size();j++)
            {
                size_t IndexOld=OldConcaveIndex[j];
                assert(VertexMap.count(IndexOld)>0);
                size_t IndexNew=VertexMap[IndexOld];
                assert(IndexNew>=0);
                assert(IndexNew<mesh.vert.size());
                Patches[i].ConcaveVertIndex.insert(IndexNew);
            }

            for(size_t j=0;j<OldFlatIndex.size();j++)
            {
                size_t IndexOld=OldFlatIndex[j];
                assert(VertexMap.count(IndexOld)>0);
                size_t IndexNew=VertexMap[IndexOld];
                assert(IndexNew>=0);
                assert(IndexNew<mesh.vert.size());
                Patches[i].FlatVertIndex.insert(IndexNew);
            }

            for(size_t j=0;j<Patches[i].PatchMesh->vert.size();j++)
            {
                size_t IndexOld=Patches[i].PatchMesh->vert[j].Q();
                assert(VertexMap.count(IndexOld)>0);
                size_t IndexNew=VertexMap[IndexOld];
                assert(IndexNew>=0);
                assert(IndexNew<mesh.vert.size());
                Patches[i].PatchMesh->vert[j].Q()=IndexNew;
            }

            for(size_t j=0;j<Patches[i].PatchMesh->face.size();j++)
            {
                size_t IndexOld=Patches[i].PatchMesh->face[j].Q();
                assert(FaceMap.count(IndexOld)>0);
                size_t IndexNew=FaceMap[IndexOld];
                assert(IndexNew>=0);
                assert(IndexNew<mesh.face.size());
                Patches[i].PatchMesh->face[j].Q()=IndexNew;
            }
        }

        UpdateFacePatches();
        UpdateBorderPatches();
    }

    void CheckVertQuality()
    {
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            size_t Q=mesh.vert[i].Q();
            assert(Q==i);
        }
    }

    void FixCorners()
    {
        if (Param.PrintDebug) {
            std::cout<<"fixing corners"<<std::endl;
        }
        CheckVertQuality();
        for (size_t i=0;i<Patches.size();i++)
        {
            if (!Patches[i].Active)continue;
            if ((Patches[i].ConvexVertIndex.size()>=3)&&
                    (Patches[i].ConvexVertIndex.size()<=6))continue;

            if (Param.PrintDebug) {
                std::cout<<"WARNING - Not good subdivision "<<Patches[i].ConvexVertIndex.size()<<" sides - FIXED!"<<std::endl;
            }

            for (size_t j=0;j<Patches[i].IndexF.size();j++)
                mesh.face[Patches[i].IndexF[j]].C()=vcg::Color4b(255,0,0,255);

            //get the currentmesh
            MeshType partition_mesh;
            GetPatchMesh(i,partition_mesh);

            //compute the angle of each vertex
            std::vector<ScalarType> Angles(partition_mesh.vert.size(),0);
            for (size_t j=0;j<partition_mesh.face.size();j++)
                for (size_t k=0;k<3;k++)
                {
                    size_t IndexV=vcg::tri::Index(partition_mesh,partition_mesh.face[j].V(k));
                    Angles[IndexV]+=vcg::face::WedgeAngleRad(partition_mesh.face[j],k);
                }

            //sort the angles
            std::vector<std::pair<ScalarType,size_t> > AngleVert;
            std::vector<std::pair<ScalarType,size_t> > AngleVertGlobal;
            for (size_t j=0;j<Angles.size();j++)
            {
                AngleVert.push_back(std::pair<ScalarType,size_t>(Angles[j],j));
                size_t IndexGlobal=partition_mesh.vert[j].Q();
                AngleVertGlobal.push_back(std::pair<ScalarType,size_t>(Angles[j],IndexGlobal));
            }

            std::sort(AngleVert.begin(),AngleVert.end());
            std::sort(AngleVertGlobal.begin(),AngleVertGlobal.end());
            std::reverse(AngleVertGlobal.begin(),AngleVertGlobal.end());

            if (Patches[i].ConvexVertIndex.size()<3)
            {

                size_t DifferenceV=3-Patches[i].ConvexVertIndex.size();
                if (Param.PrintDebug) {
                    std::cout<<"FIXED 2"<<std::endl;
                }
                //then add one by one until match
                for (size_t j=0;j<AngleVert.size();j++)
                {
                    //get the index
                    size_t currIndex=AngleVert[j].second;
                    //transform to global
                    currIndex=partition_mesh.vert[currIndex].Q();

                    if (Patches[i].ConvexVertIndex.count(currIndex)>0)continue;
                    Patches[i].ConvexVertIndex.insert(currIndex);
                    Patches[i].FlatVertIndex.erase(currIndex);
                    DifferenceV--;
                    if (DifferenceV==0)break;
                }
                //check
                assert(DifferenceV==0);//otherwise not possible fix (should be aloways possible)
            }
            else
                if (Patches[i].ConvexVertIndex.size()>6)
                {
                    size_t DifferenceV=Patches[i].ConvexVertIndex.size()-6;

                    //then add one by one until match
                    for (size_t j=0;j<AngleVertGlobal.size();j++)
                    {
                        //get the index
                        size_t currIndex=AngleVertGlobal[j].second;

                        if (Patches[i].ConvexVertIndex.count(currIndex)==0)continue;
                        Patches[i].ConvexVertIndex.erase(currIndex);
                        Patches[i].FlatVertIndex.insert(currIndex);
                        DifferenceV--;
                        if (DifferenceV==0)break;
                    }
                    assert(DifferenceV==0);//otherwise not possible fix (should be aloways possible)
                    assert(Patches[i].ConvexVertIndex.size()<7);
                }
            std::vector<size_t> Vertx(Patches[i].ConvexVertIndex.begin(),
                                      Patches[i].ConvexVertIndex.end());

            if (Param.PrintDebug) {
                std::cout<<"Vertices Fixed "<<std::endl;
                for (size_t i=0;i<Vertx.size();i++)
                    std::cout<<"Vert "<<Vertx[i]<<std::endl;
                std::cout<<"End "<<std::endl;
            }

            //        MeshType currMesh;
            //        GetPatchMesh(i,currMesh);

//#ifdef SAVE_MESHES_FOR_DEBUG
            //vcg::tri::io::ExporterOBJ<MeshType>::Save(currMesh,"results/tracer_patch2.obj", vcg::tri::io::Mask::IOM_NONE);
//#endif
        }
    }

//    void FindConnectedComponents(std::vector<std::vector<size_t> > &Components)
//    {
//        Components.clear();
//        std::set<size_t> explored;
//        //get connected components
//        for (size_t i=0;i<mesh.face.size();i++)
//        {
//            std::vector<size_t> stack;
//            size_t IndexF=i;
//            if (explored.count(IndexF)>0)continue;

//            stack.push_back(IndexF);
//            explored.insert(IndexF);
//            Components.resize(Components.size()+1);
//            do
//            {
//                size_t currF=stack.back();
//                stack.pop_back();

//                Components.back().push_back(currF);
//                for (size_t i=0;i<mesh.face[currF].VN();i++)
//                {
//                    if (vcg::face::IsBorder(mesh.face[currF],i))continue;

//                    int NextFIndex=vcg::tri::Index(mesh,mesh.face[currF].FFp(i));

//                    if (explored.count(NextFIndex)>0)continue;

//                    explored.insert(NextFIndex);
//                    stack.push_back(NextFIndex);
//                }
//            }while (!stack.empty());
//        }
//        //sort and make unique
//        for (size_t i=0;i<Components.size();i++)
//        {
//            std::sort(Components[i].begin(),Components[i].end());
//            auto last=std::unique(Components[i].begin(),Components[i].end());
//            Components[i].erase(last, Components[i].end());
//        }
//    }

    void RemoveSmallCC(ScalarType MinSize=0.001)
    {
        std::vector<std::vector<size_t> > ComponentsVect;
        ComponentsVect = findConnectedComponents(mesh);
        ScalarType BBsize=mesh.bbox.Diag();
        bool has_deleted=false;
        for (size_t i=0;i<ComponentsVect.size();i++)
        {
            vcg::Box3<ScalarType> BBox;
            for (size_t j=0;j<ComponentsVect[i].size();j++)
            {
                size_t IndexF=ComponentsVect[i][j];
                BBox.Add(mesh.face[IndexF].P(0));
                BBox.Add(mesh.face[IndexF].P(1));
                BBox.Add(mesh.face[IndexF].P(2));
            }
            if (BBox.Diag()>MinSize)continue;
            has_deleted=true;
            for (size_t j=0;j<ComponentsVect[i].size();j++)
                vcg::tri::Allocator<MeshType>::DeleteFace(mesh,mesh.face[ComponentsVect[i][j]]);
        }

        if (has_deleted)
        {
            if (Param.PrintDebug) {
                std::cout<<"WARNING: Removed Components"<<std::endl;
            }
            vcg::tri::Clean<MeshType>::RemoveUnreferencedVertex(mesh);
            vcg::tri::Allocator<MeshType>::CompactEveryVector(mesh);
            updateAllMeshAttributes(mesh);
        }

    }

    void FixSmallCC()
    {
        if (Param.PrintDebug) {
            std::cout<<"** fixing small connected components"<<std::endl;
        }

        std::vector<std::vector<size_t> > ComponentsVect;
        ComponentsVect = QuadRetopology::internal::findConnectedComponents(mesh);

        std::vector<std::set<size_t> > Components;
        for (size_t i=0;i<ComponentsVect.size();i++)
            Components.push_back(std::set<size_t>(ComponentsVect[i].begin(),
                                                  ComponentsVect[i].end()));

        //        //get connected components
        //        for (size_t i=0;i<mesh.face.size();i++)
        //        {
        //            std::vector<size_t> stack;
        //            size_t IndexF=i;
        //            if (explored.count(IndexF)>0)continue;

        //            stack.push_back(IndexF);
        //            explored.insert(IndexF);
        //            Components.resize(Components.size()+1);
        //            do
        //            {
        //                size_t currF=stack.back();
        //                stack.pop_back();

        //                Components.back().insert(currF);
        //                for (size_t i=0;i<mesh.face[currF].VN();i++)
        //                {
        //                    if (vcg::face::IsBorder(mesh.face[currF],i))continue;

        //                    int NextFIndex=vcg::tri::Index(mesh,mesh.face[currF].FFp(i));

        //                    if (explored.count(NextFIndex)>0)continue;

        //                    explored.insert(NextFIndex);
        //                    stack.push_back(NextFIndex);
        //                }
        //            }while (!stack.empty());
        //        }
        if (Param.PrintDebug)
            std::cout<<"There are "<<Components.size()<<" components"<<std::endl;

        std::vector<std::vector<size_t> > ComponentPatches;
        ComponentPatches.resize(Components.size());
        for (size_t i=0;i<Patches.size();i++)
        {
            if (!Patches[i].Active)continue;
            assert(Patches[i].IndexF.size()>0);
            size_t TestF=Patches[i].IndexF[0];
            for (size_t j=0;j<Components.size();j++)
            {
                if (Components[j].count(TestF)>0)
                    ComponentPatches[j].push_back(i);
            }
        }

        for (size_t i=0;i<ComponentPatches.size();i++)
        {
            if (ComponentPatches[i].size()>=3)continue;
            if (Param.PrintDebug)
                std::cout<<"One Patch must be splitted"<<std::endl;

            for (size_t j=0;j<ComponentPatches[i].size();j++)
            {
                size_t IndexPatch=ComponentPatches[i][j];
                if (Patches[IndexPatch].IndexF.size()==1)
                    assert(0);
                //SplitSingleFacePatch(IndexPatch);
                else
                    ForceSplit(IndexPatch);
            }
        }
        if (Param.PrintDebug) {
            std::cout<<"End"<<std::endl;
        }
        //        std::vector<std::set<size_t> > PatchComponent;
        //        for (size_t i=0;i<)
    }

    void RetrievePathesInfo(std::vector<std::vector<size_t> > &Partitions,
                            std::vector<std::vector<size_t> > &Corners,
                            std::vector<std::vector<std::vector<std::pair<size_t,size_t> > > > &BorderEdges)
    {
        UpdateFacePatches();
        UpdateBorderPatches();

        //set all faux
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<mesh.face[i].VN();j++)
                mesh.face[i].SetF(j);

        //set all border edges and  edges between patches as faux
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<mesh.face[i].VN();j++)
            {
                //if it is border clear faux
                if (vcg::face::IsBorder(mesh.face[i],j))mesh.face[i].ClearF(j);

                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
                std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
                if (BorderPatches.count(key)==0)continue;
                mesh.face[i].ClearF(j);
            }

        //select all vertices
        vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);
        for (size_t i=0;i<Patches.size();i++)
        {
            if (!Patches[i].Active)continue;
            std::vector<size_t> VertIndex(Patches[i].ConvexVertIndex.begin(),
                                          Patches[i].ConvexVertIndex.end());

            for (size_t j=0;j<VertIndex.size();j++)
                mesh.vert[VertIndex[j]].SetS();
        }

        Partitions.clear();
        Corners.clear();
        BorderEdges.clear();
        //then iterate over all patches
        for (size_t i=0;i<Patches.size();i++)
        {
            if (!Patches[i].Active)continue;
            vcg::face::Pos<FaceType> StartPos;
            bool found=false;
            Partitions.resize(Partitions.size()+1);
            Corners.resize(Corners.size()+1);
            BorderEdges.resize(BorderEdges.size()+1);

            for (size_t j=0;j<Patches[i].IndexF.size();j++)
            {
                size_t IndexF=Patches[i].IndexF[j];
                Partitions.back().push_back(IndexF);
                for (size_t k=0;k<3;k++)
                {
                    if (mesh.face[IndexF].IsF(k))continue;
                    size_t IndexV=vcg::tri::Index(mesh,mesh.face[IndexF].V(k));
                    if (Patches[i].ConvexVertIndex.count(IndexV)==0)continue;
                    StartPos=vcg::face::Pos<FaceType>(&mesh.face[IndexF],k);
                    found=true;
                    break;
                }
            }
            assert(found);
            assert(StartPos.V()->IsS());
            size_t IndexV=vcg::tri::Index(mesh,StartPos.V());
            assert(Patches[i].ConvexVertIndex.count(IndexV)>0);
            assert(!StartPos.IsFaux());

            //add the first corner
            Corners.back().push_back(IndexV);
            //create the first edge and subedge
            //add a new edge
            BorderEdges.back().resize(BorderEdges.back().size()+1);
            //add a subedge
            BorderEdges.back().back().resize(BorderEdges.back().back().size()+1);

            //then start iterating
            vcg::face::Pos<FaceType> CurrPos=StartPos;
            VertexType *StartV=StartPos.V();
            do
            {
                CurrPos.NextNotFaux();
                if (CurrPos.V()==StartV)break;

                size_t IndexF=vcg::tri::Index(mesh,CurrPos.F());
                size_t IndexE=CurrPos.E();
                BorderEdges.back().back().push_back(std::pair<size_t,size_t> (IndexF,IndexE));
                size_t IndexV=vcg::tri::Index(mesh,CurrPos.V());
                if (Patches[i].ConvexVertIndex.count(IndexV)>0)
                {
                    //add the corner
                    Corners.back().push_back(IndexV);
                    //add a new edge
                    BorderEdges.back().resize(BorderEdges.back().size()+1);
                    //add a new subedge
                    BorderEdges.back().back().resize(BorderEdges.back().back().size()+1);
                }
                else
                {
                    //in this case only add a sub edge
                    if (mesh.vert[IndexV].IsS())
                        BorderEdges.back().back().resize(BorderEdges.back().back().size()+1);
                }

            }while(true);
            if (Corners.back().size()!=Patches[i].ConvexVertIndex.size())
            {
                std::cout<<"size 0 "<<Corners.back().size()<<std::endl;
                std::cout<<"size 1 "<<Patches[i].ConvexVertIndex.size()<<std::endl;
                for (size_t j=0;j<Corners.back().size();j++)
                    std::cout<<"Corners "<<Corners.back()[j]<<std::endl;

#ifdef SAVE_MESHES_FOR_DEBUG
                vcg::tri::io::ExporterOBJ<MeshType>::Save(*Patches[i].PatchMesh,"results/tracer_patch.obj", vcg::tri::io::Mask::IOM_NONE);
#endif
                assert(0);
                //assert(Corners.back().size()==Patches[i].ConvexVertIndex.size());
                //assert(BorderEdges.back().size()==Patches[i].ConvexVertIndex.size());
            }
        }
        for (size_t i=0;i<Corners.size();i++)
        {
            assert(Corners[i].size()>=3);
            assert(Corners[i].size()<7);
        }
    }

    void FixHoles()
    {
        QuadMeshType SwapQuad;
        vcg::tri::Append<QuadMeshType,QuadMeshType>::Mesh(SwapQuad,quadMesh);

        //transform the quad mesh into tris
        //updateAllMeshAttributes()
        std::vector<int> ret=QuadRetopology::internal::splitFacesInTriangles(SwapQuad);

//#ifdef SAVE_MESHES_FOR_DEBUG
        //vcg::tri::io::ExporterOBJ<QuadMeshType>::Save(SwapQuad,"results/tracer_SwapQuad.obj", vcg::tri::io::Mask::IOM_NONE);
//#endif

        MeshType externalBound;
        vcg::tri::Append<MeshType,QuadMeshType>::Mesh(externalBound,SwapQuad);
        updateAllMeshAttributes(externalBound);
        MeshType TotalMesh;
        vcg::tri::Append<MeshType,MeshType>::Mesh(TotalMesh,externalBound);
        vcg::tri::Append<MeshType,MeshType>::Mesh(TotalMesh,mesh);
        //vcg::tri::Clean<MeshType>::RemoveDuplicateVertex(TotalMesh);
        vcg::tri::Clean<MeshType>::MergeCloseVertex(TotalMesh,0.000001);
        updateAllMeshAttributes(TotalMesh);

//#ifdef SAVE_MESHES_FOR_DEBUG
        //vcg::tri::io::ExporterOBJ<MeshType>::Save(TotalMesh,"results/tracer_total.obj", vcg::tri::io::Mask::IOM_NONE);
//#endif

        size_t num_holes=vcg::tri::Clean<MeshType>::CountHoles(TotalMesh);
        if (num_holes==0)return;

        std::cout<<"*** There are "<<num_holes<<" HOLES to be fixed!"<<std::endl;

        //first smooth around borders
        //vcg::tri::Smooth<MeshType>::VertexCoordLaplacian()
        std::cout<<"*** filling holes"<<std::endl;
        size_t FNum0=TotalMesh.face.size();
        vcg::tri::Hole<MeshType>::template EarCuttingFill<vcg::tri::TrivialEar<MeshType> >(TotalMesh,TotalMesh.face.size(),false);
        size_t FNum1=TotalMesh.face.size();
        std::cout<<"*** Done"<<std::endl;
        //create a new mesh only with the new
        vcg::tri::UpdateSelection<MeshType>::Clear(TotalMesh);
        for (size_t i=FNum0;i<FNum1;i++)
            TotalMesh.face[i].SetS();
        vcg::tri::UpdateSelection<MeshType>::VertexFromFaceLoose(TotalMesh);


        MeshType NewHoles;
        vcg::tri::Append<MeshType,MeshType>::Mesh(NewHoles,TotalMesh,true);
        updateAllMeshAttributes(NewHoles);

//#ifdef SAVE_MESHES_FOR_DEBUG
        //vcg::tri::io::ExporterOBJ<MeshType>::Save(NewHoles,"results/tracer_holes.obj", vcg::tri::io::Mask::IOM_NONE);
//#endif

        std::vector<std::vector<size_t> > Components;
        Components = QuadRetopology::internal::findConnectedComponents(NewHoles);
        std::cout<<"*** Components "<<Components.size()<<std::endl;

        assert(Components.size()==num_holes);

        for (size_t i=0;i<Components.size();i++)
        {
            Patches.resize(Patches.size()+1);
            for (size_t j=0;j<Components[i].size();j++)
            {
                size_t IndexF=Components[i][j];
                CoordType pos0=NewHoles.face[IndexF].P(0);
                CoordType pos1=NewHoles.face[IndexF].P(1);
                CoordType pos2=NewHoles.face[IndexF].P(2);
                size_t num0=mesh.vert.size();
                vcg::tri::Allocator<MeshType>::AddVertices(mesh,3);
                mesh.vert[num0].P()=pos0;
                //mesh.vert[num0].Q()=num0;
                mesh.vert[num0+1].P()=pos1;
                //mesh.vert[num0+1].Q()=num0+1;
                mesh.vert[num0+2].P()=pos2;
                //mesh.vert[num0+2].Q()=num0+2;
                vcg::tri::Allocator<MeshType>::AddFace(mesh,&mesh.vert[num0],
                                                       &mesh.vert[num0+1],
                                                       &mesh.vert[num0+2]);
                mesh.face.back().Q()=mesh.face.size()-1;
                Patches.back().IndexF.push_back(mesh.face.size()-1);
                Patches.back().Active=true;
            }
            InitPatchMesh(Patches.size()-1);
        }
        vcg::tri::Clean<MeshType>::RemoveDuplicateVertex(mesh);
        vcg::tri::Allocator<MeshType>::CompactEveryVector(mesh);

        updateAllMeshAttributes(mesh);

        size_t num_split=vcg::tri::Clean<MeshType>::SplitNonManifoldVertex(mesh,0);
        if (num_split>0)
            updateAllMeshAttributes(mesh);

        InitIndexOnQ();
    }

public:


#ifdef DRAWTRACE
    void GLDrawPatches()
    {

        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.99998);
        glLineWidth(10);
        for (size_t i=0;i<mesh.face.size();i++)
        {
            for (size_t j=0;j<mesh.face[i].VN();j++)
            {
                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
                std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
                bool IsBorder=vcg::face::IsBorder(mesh.face[i],j);
                if ((BorderPatches.count(key)==0)&&(!IsBorder))continue;

                vcg::glColor(vcg::Color4b(0,0,0,255));

                CoordType Pos0=mesh.face[i].P0(j);
                CoordType Pos1=mesh.face[i].P1(j);

                glBegin(GL_LINES);
                vcg::glVertex(Pos0);
                vcg::glVertex(Pos1);
                glEnd();
            }
        }
        glPopAttrib();
//        for (size_t i=0;i<Patches.size();i++)
//        {
//            if (!Patches[i].Active)continue;
//            GLDrawPatch(i);
//            //GLDrawPatchDir(i);
//        }
    }
#endif

    void BatchProcess(std::vector<std::vector<size_t> > &Partitions,
                      std::vector<std::vector<size_t> > &Corners,
                      std::vector<std::vector<std::vector<std::pair<size_t,size_t> > > > &BorderEdges,
                      Parameters &_Param)
    {
        Param=_Param;
        assert(Param.EdgeSizeFactor>=0.5);
        if (Param.PrintDebug) {
            std::cout<<"*****************************"<<std::endl;
            std::cout<<"***                       ***"<<std::endl;
            std::cout<<"***   GETTING GEO-CONSTR  ***"<<std::endl;
            std::cout<<"***                       ***"<<std::endl;
            std::cout<<"*****************************"<<std::endl;
        }
        //initial check
        updateAllMeshAttributes(mesh);
        int nonManifV=vcg::tri::Clean<MeshType>::CountNonManifoldVertexFF(mesh);
        int nonManifE=vcg::tri::Clean<MeshType>::CountNonManifoldEdgeFF(mesh);

        if (nonManifV > 0)
            std::cout<<"non manuf V:"<<nonManifV<<std::endl;
        if (nonManifE > 0)
            std::cout<<"non manuf E:"<<nonManifE<<std::endl;

        std::vector<CoordType> ConvexV,ConcaveV;
        GetGeometricCorners(ConvexV,ConcaveV);

        SelectFacesWithBorderEdge();
        vcg::PolygonalAlgorithm<MeshType>::Triangulate(mesh,true,true);
        updateAllMeshAttributes(mesh);

        OriginalMesh.Clear();
        vcg::tri::Append<MeshType,MeshType>::Mesh(OriginalMesh,mesh);
        updateAllMeshAttributes(OriginalMesh);
        Grid.Set(OriginalMesh.face.begin(),OriginalMesh.face.end());

        if (Param.InitialRemesh)
            RemeshMesh();

        //RemoveSmallCC();
        //    else
        //        SmoothMesh();

        size_t num_split=vcg::tri::Clean<MeshType>::SplitNonManifoldVertex(mesh,0);
        if (num_split>0)
            updateAllMeshAttributes(mesh);

        //initialize the field
        InitField();
        InitIndexOnQ();
        //return;
        SetInitialPatch(ConvexV,ConcaveV);

        InitPatchMesh(0);

        //split the concave patches first
        if (Param.PrintDebug) {
            std::cout<<"*****************************"<<std::endl;
            std::cout<<"***                       ***"<<std::endl;
            std::cout<<"*** SPLITTING CONCAVE     ***"<<std::endl;
            std::cout<<"***                       ***"<<std::endl;
            std::cout<<"*****************************"<<std::endl;
        }
        SplitConcave(0);

        bool splitted=false;
        do{
            splitted=false;
            for (size_t i=0;i<Patches.size();i++)
            {
                if (!Patches[i].Active)continue;
                if (Patches[i].ConcaveVertIndex.size()==0)continue;
                splitted|=SplitConcave(i);
            }
        }while (splitted);

        //safety check.. all concave non solved becomes flat
        for (size_t i=0;i<Patches.size();i++)
        {
            if (!Patches[i].Active)continue;
            if (Patches[i].ConcaveVertIndex.size()==0)continue;

            std::set<size_t>::iterator IteSetConcave;

            for (IteSetConcave=Patches[i].ConcaveVertIndex.begin();
                 IteSetConcave!=Patches[i].ConcaveVertIndex.end();
                 IteSetConcave++)
            {
                if (Param.PrintDebug) {
                    std::cout<<"WARNING - Non Traced Concave - MADE FLAT!"<<std::endl;
                }
                size_t IndConcave=(*IteSetConcave);
                //Patches[i].ConcaveVertIndex.erase(IndConcave);
                Patches[i].FlatVertIndex.insert(IndConcave);
            }
            Patches[i].ConcaveVertIndex.clear();

            //set to size 0
            for (size_t j=0;j<Patches[i].TracingDir.size();j++)
            {
                if (Patches[i].TracingDir[j].size()<2)continue;
                Patches[i].TracingDir[j].clear();
            }
        }

        if (Param.PrintDebug) {
            std::cout<<"*****************************"<<std::endl;
            std::cout<<"***                       ***"<<std::endl;
            std::cout<<"*** SPLITTING NON OK STEP ***"<<std::endl;
            std::cout<<"***                       ***"<<std::endl;
            std::cout<<"*****************************"<<std::endl;
        }

        //split the rest
        splitted=false;
        do{
            splitted=SplitNonOk();
        }while (splitted);


        if (Param.PrintDebug) {
            std::cout<<"*****************************"<<std::endl;
            std::cout<<"***                       ***"<<std::endl;
            std::cout<<"*** COLLECTING RESULTS    ***"<<std::endl;
            std::cout<<"***                       ***"<<std::endl;
            std::cout<<"*****************************"<<std::endl;
        }

        //final check
        nonManifV=vcg::tri::Clean<MeshType>::CountNonManifoldVertexFF(mesh);
        nonManifE=vcg::tri::Clean<MeshType>::CountNonManifoldEdgeFF(mesh);
        if (nonManifV > 0)
            std::cout<<"non manuf V:"<<nonManifV<<std::endl;
        if (nonManifE > 0)
            std::cout<<"non manuf E:"<<nonManifE<<std::endl;

        //vcg::tri::UpdateColor<MeshType>::PerFaceConstant(mesh);

        if (Param.PrintDebug) {
            std::cout<<"*****************************"<<std::endl;
            std::cout<<"***                       ***"<<std::endl;
            std::cout<<"***FIXING REMAINING ISSUES***"<<std::endl;
            std::cout<<"***                       ***"<<std::endl;
            std::cout<<"*****************************"<<std::endl;
        }

        RemoveNotTopologicallyOKPartitions();
        //FixSmallCC();
        //FixHoles();

        FixCorners();

        FixedVert.clear();
        if (Param.FinalSmooth)
            SmoothPatches(0.9,20,Param.Reproject);
        if (Param.CleanFolds && Param.InterleaveSmooth)
            CleanFoldsNeeded();

        if (Param.PrintDebug) {
            std::cout<<"*****************************"<<std::endl;
            std::cout<<"***                       ***"<<std::endl;
            std::cout<<"*** RETRIEVE PARTITIONS   ***"<<std::endl;
            std::cout<<"***                       ***"<<std::endl;
            std::cout<<"*****************************"<<std::endl;
        }

        RetrievePathesInfo(Partitions,Corners,BorderEdges);

        ColorPatches();

        nonManifV=vcg::tri::Clean<MeshType>::CountNonManifoldVertexFF(mesh);
        nonManifE=vcg::tri::Clean<MeshType>::CountNonManifoldEdgeFF(mesh);
        if (nonManifV > 0)
            std::cout<<"non manuf V:"<<nonManifV<<std::endl;
        if (nonManifE > 0)
            std::cout<<"non manuf E:"<<nonManifE<<std::endl;

    }

    PatchAssembler(MeshType &_mesh,QuadMeshType &_quadMesh):mesh(_mesh),quadMesh(_quadMesh)
    {}

    ~PatchAssembler()
    {
        for (size_t i=0;i<Patches.size();i++)
        {
            delete(Patches[i].PatchMesh);
            delete(Patches[i].PatchFieldTracer);
        }
    }
};


}
}
#endif //qr_PATCH_ASSEMBLER_H
