#pragma once
#include <vector>
#include <utility>
#include <functional>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>

namespace kt84 {

template <class MeshTrait>
inline bool append_quad_strip(                                          // Return: true if successful
    OpenMesh::PolyMesh_ArrayKernelT<MeshTrait>& mesh,                   // Target mesh
    OpenMesh::HalfedgeHandle                  & boundary_h_front,       // Halfedge where the first quad is attached to. Updated to the corresponding newly generated halfedge.
    OpenMesh::HalfedgeHandle                  & boundary_h_back)        // Halfedge where the last  quad is attached to. Updated to the corresponding newly generated halfedge.
{
    typedef OpenMesh::PolyMesh_ArrayKernelT<MeshTrait> Mesh;
    
    if (!mesh.is_boundary(boundary_h_front) || !mesh.is_boundary(boundary_h_back))
        return false;
    
    OpenMesh::VertexHandle v_prev, v_next;
    for (auto h = boundary_h_front; ; ) {
        auto h_next = mesh.next_halfedge_handle(h);
        auto v_from = mesh.from_vertex_handle(h);
        auto v_to   = mesh.to_vertex_handle  (h);
        if (h == boundary_h_front)
            v_prev = mesh.add_vertex(typename Mesh::Point());
        v_next = mesh.add_vertex(typename Mesh::Point());
        mesh.add_face(v_from, v_to, v_next, v_prev);
        if (h == boundary_h_back)
            break;
        v_prev = v_next;
        h      = h_next;
    }
    
    boundary_h_front = mesh.opposite_halfedge_handle(mesh.next_halfedge_handle(mesh.next_halfedge_handle(boundary_h_front)));
    boundary_h_back  = mesh.opposite_halfedge_handle(mesh.next_halfedge_handle(mesh.next_halfedge_handle(boundary_h_back )));
    
    return true;
}

}
