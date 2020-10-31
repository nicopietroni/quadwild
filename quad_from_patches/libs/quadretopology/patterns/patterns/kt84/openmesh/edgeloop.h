#pragma once
#include <vector>
#include <utility>
#include <functional>
#include <OpenMesh/Core/Mesh/Handles.hh>

namespace kt84 {

//template <class Mesh>     // Should be done with template alias in C++11
//struct EdgeLoop_VertexFunc : public std::function<void(Mesh&, OpenMesh::HalfedgeHandle, OpenMesh::VertexHandle)> {};


template <class Mesh>
inline void insert_edgeloop(                                     // Return: array of newly generated vertices
    Mesh&                       mesh,                                                           // Target mesh, assumed to be quad-only
    OpenMesh::HalfedgeHandle    h_start,                                                        // Boundary halfedge where the edgeloop starts
    double                      t = 0.5,                                                        // Position of splitting edges in the quad strip
    std::function<void(Mesh&, OpenMesh::HalfedgeHandle, OpenMesh::VertexHandle)> 
    func = std::function<void(Mesh&, OpenMesh::HalfedgeHandle, OpenMesh::VertexHandle)>())      // Function called when making new vertex
{
    // check if this edge loop forms a loop by walking toward h1 of e
    for (auto h = mesh.opposite_halfedge_handle(h_start); ; ) {
        auto h_opposite           = mesh.opposite_halfedge_handle(h);
        auto h_next_next          = mesh.next_halfedge_handle(mesh.next_halfedge_handle(h));
        auto h_next_next_opposite = mesh.opposite_halfedge_handle(h_next_next);
        
        if (mesh.is_boundary(h)) {
            // reached a boundary
            h_start = h_opposite;
            break;
        } else if (h_start == h_next_next) {
            // loop detected!
            break;
        }
        h = h_next_next_opposite;
    }
        
    // for each edge, insert a vertex at its middle
    for (auto h = h_start; ;) {
        auto h_next =
            mesh.opposite_halfedge_handle(
            mesh.next_halfedge_handle    (
            mesh.next_halfedge_handle    (h)));
            
        auto p0 = mesh.point(mesh.from_vertex_handle(h));
        auto p1 = mesh.point(mesh.  to_vertex_handle(h));
        auto v_new = mesh.add_vertex((1 - t) * p0 + t * p1);
        if (func)
            func(mesh, h, v_new);

        mesh.split_edge(mesh.edge_handle(h), v_new);
            
        if (mesh.is_boundary(h))
            break;
            
        if (h_next == h_start)
            break;
            
        h = h_next;
    }
        
    // connect inserted vertices
    for (auto h = h_start; ;) {
        auto h_prev =
            mesh.next_halfedge_handle(
            mesh.next_halfedge_handle(h));
            
        auto h_next = mesh.opposite_halfedge_handle(h_prev);
            
        mesh.insert_edge(h_prev, h);
            
        h = h_next;
        if (mesh.is_boundary(h) || h == h_start)
            break;
    }
}

template <class Mesh>
inline std::vector<OpenMesh::HalfedgeHandle> find_edgeloop(     // Return : array of halfedges belonging to the quad strip, starting with h_start
    Mesh&                       mesh,                           // Target mesh, assumed to be quad-only
    OpenMesh::HalfedgeHandle    h_start)                        // Boundary halfedge where the edgeloop starts
{
    // check if this edge loop forms a loop by walking toward h1 of e
    for (auto h = mesh.opposite_halfedge_handle(h_start); ; ) {
        auto h_opposite           = mesh.opposite_halfedge_handle(h);
        auto h_next_next          = mesh.next_halfedge_handle(mesh.next_halfedge_handle(h));
        auto h_next_next_opposite = mesh.opposite_halfedge_handle(h_next_next);
            
        if (mesh.is_boundary(h)) {
            // reached a boundary
            h_start = h_opposite;
            break;
        } else if (h_start == h_next_next) {
            // loop detected!
            break;
        }
        h = h_next_next_opposite;
    }
        
    // collect result
    std::vector<OpenMesh::HalfedgeHandle> result;
    for (auto h = h_start; ;) {
        result.push_back(h);
            
        if (mesh.is_boundary(h))
            break;
            
        auto h_next =
            mesh.opposite_halfedge_handle(
            mesh.next_halfedge_handle    (
            mesh.next_halfedge_handle    (h)));
            
        if (h_next == h_start) {
            result.push_back(h_start);          // indication of loop
            break;
        }
            
        h = h_next;
    }
        
    return result;
}

}
