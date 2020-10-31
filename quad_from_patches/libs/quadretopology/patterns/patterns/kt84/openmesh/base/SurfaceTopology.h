#pragma once
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include "../DerivedPtrHolder.h"
#include "Utility.h"        // assumption: TMesh is derived from TMeshBase and Utility<TMeshBase, TMesh>

namespace kt84 {

template <class TMeshBase, class TMesh>
struct SurfaceTopology : public DerivedPtrHolder<TMesh, SurfaceTopology<TMeshBase, TMesh>> {
    struct Data {
        int num_connected_components;
        int num_boundary_loops;
        Data()
            : num_connected_components()
            , num_boundary_loops()
        {}
        bool is_closed() const { return num_boundary_loops == 0; }
    } surfaceTopology;
    
    void surfaceTopology_compute() {
        TMesh* mesh = DerivedPtrHolder<TMesh, SurfaceTopology<TMeshBase, TMesh>>::derived_ptr;
        
        // count connected components
        {
            boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graph;
            for (auto f : mesh->faces())
                boost::add_vertex(graph);
            for (auto e : mesh->edges()) {
                if (mesh->is_boundary(e)) continue;
                auto h0h1 = mesh->util_edge_to_halfedge_pair(e);
                auto f0 = mesh->face_handle(h0h1.first );
                auto f1 = mesh->face_handle(h0h1.second);
                boost::add_edge(f0.idx(), f1.idx(), graph);
            }
            std::vector<int> f_component_id(mesh->n_faces());
            surfaceTopology.num_connected_components = boost::connected_components(graph, &f_component_id[0]);
        }
        
        // count boundary loops
        {
            boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graph;
            for (auto v : mesh->vertices()) {
                if (!mesh->is_boundary(v)) continue;
                auto h0 = mesh->halfedge_handle(v);
                auto h1 = mesh->prev_halfedge_handle(h0);
                auto e0 = mesh->edge_handle(h0);
                auto e1 = mesh->edge_handle(h1);
                boost::add_edge(e0.idx(), e1.idx(), graph);
            }
            std::vector<int> e_component_id(mesh->n_edges());
            surfaceTopology.num_boundary_loops = boost::connected_components(graph, &e_component_id[0]);
        }
    }
};

}
