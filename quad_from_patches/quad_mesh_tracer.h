#ifndef QUAD_MESH_TRACER_H
#define QUAD_MESH_TRACER_H

#include <set>
#include <vector>

#include <Eigen/Core>

#include <vcg/complex/complex.h>
#include <wrap/io_trimesh/export_obj.h>

template <class PolyMeshType>
class QuadMeshTracer
{
    typedef typename PolyMeshType::FaceType PFace;

    PolyMeshType &PolyM;

    std::deque<vcg::face::Pos<PFace> > TracingStack;
    std::set<size_t> VisitedVertices;
    std::set<std::pair<size_t,size_t> > TracedEdges;

public:
    QuadMeshTracer(PolyMeshType &_PolyM):PolyM(_PolyM){}

    bool MotorCycle=false;
    std::vector<int> FacePatch;
    bool DebugMessages=false;

    void updatePolymeshAttributes();
    void InitSingularities();
    int SplitIntoPatches();
    void DoTrace();

    void ColorMesh();

    void SaveColoredMesh(std::string PatchPath);
    void SavePatchDecomposition(std::string &PatchFile);
    void SaveTracedEdgeMesh(std::string &EdgePath);

    void TracePartitions();
};

#include "quad_mesh_tracer.cpp"

#endif // QUAD_MESH_TRACER_H
