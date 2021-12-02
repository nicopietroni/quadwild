#pragma once
#include "../DerivedPtrHolder.h"

namespace kt84 {

struct EdgeLength_EdgeTraits {
    double edgeLength;
    EdgeLength_EdgeTraits() : edgeLength() {}
};

template <class TMeshBase, class TMesh>
struct EdgeLength : public DerivedPtrHolder<TMesh, EdgeLength<TMeshBase, TMesh>> {
    void edgeLength_compute() {
        TMesh* mesh = DerivedPtrHolder<TMesh, EdgeLength<TMeshBase, TMesh>>::derived_ptr;
        for (auto e : mesh->edges()) {
            auto p0 = mesh->point(mesh->to_vertex_handle(mesh->halfedge_handle(e, 0)));
            auto p1 = mesh->point(mesh->to_vertex_handle(mesh->halfedge_handle(e, 1)));
            mesh->data(e).edgeLength = (p1 - p0).norm();
        }
    }
};

}
