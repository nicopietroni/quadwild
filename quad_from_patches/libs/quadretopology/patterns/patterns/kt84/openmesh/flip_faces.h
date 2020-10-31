#pragma once
#include <vector>
#include <algorithm>

namespace kt84 {

template <class Mesh>
inline void flip_faces(Mesh& mesh) {
    auto copied = mesh;
    mesh.clear();
        
    // copy vertices
    for (auto v : copied.vertices())
        mesh.data(mesh.add_vertex(copied.point(v))) = copied.data(v);
        
    // copy faces
    for (auto f : copied.faces()) {
        std::vector<typename Mesh::VHandle> face_vhandles;
        for (auto fv = copied.fv_iter(f); fv.is_valid(); ++fv)
            face_vhandles.push_back(mesh.vertex_handle(fv->idx()));
        // flip face
        std::reverse(face_vhandles.begin(), face_vhandles.end());
        mesh.data(mesh.add_face(face_vhandles)) = copied.data(f);
    }
    // TODO: copy edge/halfedge data?
}

}
