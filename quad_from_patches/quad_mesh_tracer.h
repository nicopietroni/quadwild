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
