#pragma once
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <map>
#include <vector>

namespace kt84 {

template <class Mesh>
inline void merge_nearby_vertices(Mesh& mesh, double axis_aligned_distance) {
    // WARNING! This doesn't take the manifoldness into account, could possibly result in some holes.
    
    typedef typename Mesh::Point Point;
    typedef OpenMesh::VertexHandle VHandle;
        
    // map from Point to VHandle for detecting nearby vertices
    auto pred = [&] (const Point& lhs, const Point& rhs) {
        for (int i = 0; i < Point::size_; ++i) {
            if (std::abs(lhs[i] - rhs[i]) <= axis_aligned_distance)
                continue;
            return lhs[i] < rhs[i];
        }
        return false;
    };
    std::map<Point, VHandle, decltype(pred)> point2vhandle(pred);
        
    // map from old vertex indices to new vertex indices
    std::vector<int> index_map(mesh.n_vertices(), -1);
        
    // detect nearby vertices
    for (auto v : mesh.vertices()) {
        auto p = mesh.point(v);
        auto found = point2vhandle.find(p);
        if (found == point2vhandle.end()) {
            point2vhandle.insert(std::make_pair(p, v));
            index_map[v.idx()] = v.idx();
        } else {
            index_map[v.idx()] = found->second.idx();
        }
    }
        
    // make sure that the mesh can store deletion flag
    mesh.request_vertex_status();
    mesh.request_face_status();
    mesh.request_edge_status();
    mesh.request_halfedge_status();
        
    // update all faces
    std::vector<std::vector<VHandle>> updated_faces;
    for (auto f : mesh.faces()) {
        std::vector<VHandle> face_vhandles;
        for (auto fv = mesh.fv_iter(f); fv.is_valid(); ++fv)
            face_vhandles.push_back(mesh.vertex_handle(index_map[fv->idx()]));
        updated_faces.push_back(face_vhandles);
        mesh.delete_face(f, false);
    }
    for (auto& face_vhandles : updated_faces) {
        mesh.add_face(face_vhandles);
    }
        
    // delete merged vertices
    for (auto i = 0u; i < mesh.n_vertices(); ++i) {
        if (index_map[i] != i)
            mesh.delete_vertex(mesh.vertex_handle(i));
    }
        
    mesh.garbage_collection();
}

}
