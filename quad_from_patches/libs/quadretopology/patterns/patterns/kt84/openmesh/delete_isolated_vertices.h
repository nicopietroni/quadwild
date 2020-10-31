#pragma once
#include <vector>
#include <unordered_map>

namespace kt84 {
    template <class TMesh>
    inline TMesh delete_isolated_vertices(const TMesh& mesh) {
        TMesh result;
        std::unordered_map<int, typename TMesh::VHandle> vertex_map;
        for (auto f : mesh.faces()) {
            std::vector<typename TMesh::VHandle> face_vhandles; 
            for (auto v = mesh.cfv_iter(f); v.is_valid(); ++v) {
                auto found = vertex_map.find(v->idx());
                if (found == vertex_map.end()) {
                    auto v_new = result.add_vertex(mesh.point(*v));
                    vertex_map[v->idx()] = v_new;
                    face_vhandles.push_back(v_new);
                } else {
                    face_vhandles.push_back(vertex_map[v->idx()]);
                }
            }
            result.add_face(face_vhandles);
        }
        return result;
    }
}
