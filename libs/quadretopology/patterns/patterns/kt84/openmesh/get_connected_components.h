#pragma once
#include <vector>
#include <unordered_map>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

namespace kt84 {
    template <class TMesh>
    inline std::vector<TMesh> get_connected_components(const TMesh& mesh) {
        assert(mesh.has_vertex_status() && mesh.has_face_status());
        
        std::vector<TMesh> result;
        
        boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graph;
        for (auto e : mesh.edges()) {
            int f0 = mesh.face_handle(mesh.halfedge_handle(e, 0)).idx();
            int f1 = mesh.face_handle(mesh.halfedge_handle(e, 1)).idx();
            boost::add_edge(f0, f1, graph);
        }
        
        std::vector<int> f_component_id(mesh.n_faces());
        int n_components = boost::connected_components(graph, &f_component_id[0]);
        
        // generate array of meshes according to component id
        for (int i = 0; i < n_components; ++i) {
            TMesh component = mesh;
            for (auto f : component.faces()) {
                if (f_component_id[f.idx()] != i)
                    component.delete_face(f, true);
            }
            component.garbage_collection();
            result.push_back(component);
        }
        
        return result;
    }
}
