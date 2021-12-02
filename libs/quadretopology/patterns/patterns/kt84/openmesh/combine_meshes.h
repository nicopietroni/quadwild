#pragma once
#include <vector>

namespace kt84 {
    template <class TMesh>
    inline TMesh combine_meshes(const std::vector<TMesh>& meshes) {
        TMesh result;
        int offset = 0;
        for (auto& mesh : meshes) {
            // add vertices
            for (auto v : mesh.vertices())
                result.add_vertex(mesh.point(v));
            
            // add faces
            for (auto f : mesh.faces()) {
                std::vector<typename TMesh::VHandle> face_vhandles;
                for (auto v = mesh.cfv_iter(f); v.is_valid(); ++v)
                    face_vhandles.push_back(result.vertex_handle(offset + v->idx()));
                result.add_face(face_vhandles);
            }
            
            offset += mesh.n_vertices();
        }
        return result;
    }
}
