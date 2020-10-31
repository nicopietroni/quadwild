#include "qr_patch_tracer.h"

#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/complex/algorithms/clean.h>


namespace QuadRetopology {
namespace internal {

template <class PolyMeshType>
void QuadMeshTracer<PolyMeshType>::updatePolymeshAttributes()
{
    vcg::tri::Clean<PolyMeshType>::RemoveDuplicateVertex(PolyM);
    vcg::tri::Clean<PolyMeshType>::RemoveUnreferencedVertex(PolyM);
    vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(PolyM);
    vcg::tri::UpdateFlags<PolyMeshType>::VertexBorderFromFaceAdj(PolyM);
}

template <class PolyMeshType>
void QuadMeshTracer<PolyMeshType>::InitSingularities()
{
    TracedEdges.clear();
    TracingStack.clear();
    vcg::tri::UpdateQuality<PolyMeshType>::VertexValence(PolyM);
    for (size_t i=0;i<PolyM.face.size();i++) {
        if (PolyM.face[i].VN() == 4) {
            for (int j=0;j<PolyM.face[i].VN();j++)
            {
                if (PolyM.face[i].V(j)->IsB())continue;
                if (PolyM.face[i].V(j)->Q() == 4)continue;

                vcg::face::Pos<PFace> CurrP(&PolyM.face[i],j);

                size_t IndexV=vcg::tri::Index(PolyM,CurrP.V());
                VisitedVertices.insert(IndexV);

                CurrP.FlipV();
                TracingStack.push_back(CurrP);
            }
        }
        else {
            for (int j=0;j<PolyM.face[i].VN();j++)
            {
                vcg::face::Pos<PFace> CurrP(&PolyM.face[i],j);

                size_t IndexV=vcg::tri::Index(PolyM,CurrP.V());
                VisitedVertices.insert(IndexV);

                size_t IndexV1=vcg::tri::Index(PolyM,CurrP.VFlip());
                assert(IndexV1!=IndexV);
                // add to traced paths
                std::pair<size_t,size_t> key(std::min(IndexV,IndexV1),std::max(IndexV,IndexV1));
                TracedEdges.insert(key);

                typename PolyMeshType::FaceType* start = CurrP.F();
                do {
                    CurrP.FlipF();
                    CurrP.FlipE();
                    CurrP.FlipV();
                    if (start != CurrP.F()) {
                        TracingStack.push_back(CurrP);
                    }
                    CurrP.FlipV();
                }
                while (start != CurrP.F());
            }
        }
    }
    if (DebugMessages)
        std::cout<<"Inserted "<<TracingStack.size()<<" portals"<<std::endl;
}

template <class PolyMeshType>
int QuadMeshTracer<PolyMeshType>::SplitIntoPatches()
{
    FacePatch.resize(PolyM.face.size(),-1);
    std::vector<size_t> StackMoves;
    vcg::tri::UpdateFlags<PolyMeshType>::FaceClearS(PolyM);
    int currPartition=0;
    while (true)
    {
        //find next seed
        for (size_t i=0;i<PolyM.face.size();i++)
        {
            if (PolyM.face[i].IsS())
                continue;
            StackMoves.push_back(i);
            break;
        }
        if (StackMoves.empty())
            return currPartition;
        if (DebugMessages)
            std::cout<<"Partition "<<currPartition<<std::endl;
        while (!StackMoves.empty())
        {
            size_t currF=StackMoves.back();
            PolyM.face[currF].SetS();
            StackMoves.pop_back();
            FacePatch[currF]=currPartition;
            for (int j=0;j<PolyM.face[currF].VN();j++)
            {
                size_t IndexV0=vcg::tri::Index(PolyM,PolyM.face[currF].V0(j));
                size_t IndexV1=vcg::tri::Index(PolyM,PolyM.face[currF].V1(j));
                std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
                //found a border edge
                if (TracedEdges.count(key)>0)
                    continue;
                // go on the other size
                size_t nextF=vcg::tri::Index(PolyM,PolyM.face[currF].FFp(j));
                if (PolyM.face[nextF].IsS())
                    continue;
                if (nextF==currF)
                    continue;
                StackMoves.push_back(nextF);
            }
        }
        currPartition++;
    }
    return -1;
}

template <class PolyMeshType>
void QuadMeshTracer<PolyMeshType>::DoTrace()
{
    if (DebugMessages)
        std::cout<<TracingStack.size() << " need to be traced! "<<std::endl;

    assert(!TracingStack.empty());
    int step=0;
    do{
        vcg::face::Pos<PFace> CurrP=TracingStack.front();
        TracingStack.pop_front();
        size_t IndexV0=vcg::tri::Index(PolyM,CurrP.V());
        size_t IndexV1=vcg::tri::Index(PolyM,CurrP.VFlip());
        assert(IndexV1!=IndexV0);

        // add to traced paths
        std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
        TracedEdges.insert(key);

        //already visited
        if (VisitedVertices.count(IndexV0)>0)
            continue;

        //border reached
        if (PolyM.vert[IndexV0].IsB())
            continue;

        //must been visited already
        if (MotorCycle)
        {
            assert(VisitedVertices.count(IndexV1)>0);
        }

        if (MotorCycle)
            VisitedVertices.insert(IndexV0);

        //then go to next one
        CurrP.FlipE();
        CurrP.FlipF();
        CurrP.FlipE();
        CurrP.FlipV();

        TracingStack.push_front(CurrP);

        step++;
        if (((step%100)==0)&&(DebugMessages))
            std::cout<<"Traced "<<step<<" separatrices"<<std::endl;
    }while (!TracingStack.empty());
}

template <class PolyMeshType>
void QuadMeshTracer<PolyMeshType>::ColorMesh()
{
    int NumPartitions=0;
    for (size_t i=0;i<FacePatch.size();i++)
        NumPartitions=std::max(NumPartitions,FacePatch[i]);

    for (size_t i=0;i<PolyM.face.size();i++)
    {
        assert(FacePatch[i]>=0);
        vcg::Color4b currC=vcg::Color4b::Scatter(NumPartitions,FacePatch[i]);
        PolyM.face[i].C()=currC;
    }
}

template <class PolyMeshType>
void QuadMeshTracer<PolyMeshType>::SaveColoredMesh(std::string PatchPath)
{
    ColorMesh();

    vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(PolyM,PatchPath.c_str(),vcg::tri::io::Mask::IOM_FACECOLOR);
}

template <class PolyMeshType>
void QuadMeshTracer<PolyMeshType>::SavePatchDecomposition(std::string &PatchFile)
{
    FILE *f=fopen(PatchFile.c_str(),"wt");
    for (size_t i=0;i<FacePatch.size();i++)
        fprintf(f,"%d\n",FacePatch[i]);
    fclose(f);
}

template <class PolyMeshType>
void QuadMeshTracer<PolyMeshType>::SaveTracedEdgeMesh(std::string &EdgePath)
{
    PolyMeshType TestEdge;
    std::set<std::pair<size_t,size_t> >::iterator iteE;
    for (iteE=TracedEdges.begin();iteE!=TracedEdges.end();iteE++)
    {
        typename PolyMeshType::CoordType pos0=PolyM.vert[(*iteE).first].P();
        typename PolyMeshType::CoordType pos1=PolyM.vert[(*iteE).second].P();

        vcg::tri::Allocator<PolyMeshType>::AddEdge(TestEdge,pos0,pos1);
    }
    vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(TestEdge,EdgePath.c_str(),vcg::tri::io::Mask::IOM_EDGEINDEX);
}

template <class PolyMeshType>
void QuadMeshTracer<PolyMeshType>::TracePartitions()
{
    InitSingularities();

    DoTrace();

    SplitIntoPatches();
}

}
}
