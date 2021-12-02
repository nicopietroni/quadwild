#pragma once
#include <limits>
#include <algorithm>
#include "../DerivedPtrHolder.h"

namespace kt84 {

struct NormalVariation_VertexTraits {
    double normalVariation;
    NormalVariation_VertexTraits()
        : normalVariation()
    {}
};

template <class TMeshBase, class TMesh>
struct NormalVariation : public DerivedPtrHolder<TMesh, NormalVariation<TMeshBase, TMesh>> {
    void normalVariation_compute() {
        TMesh* mesh = DerivedPtrHolder<TMesh, NormalVariation<TMeshBase, TMesh>>::derived_ptr;
        
        for (auto v : mesh->vertices()) {
            double& normalVariation = mesh->data(v).normalVariation;
            auto nv = mesh->normal(v);

            normalVariation = std::numeric_limits<double>::max();
            
            for (auto f = mesh->vf_iter(v); f.is_valid(); ++f) {
                auto nf = mesh->normal(*f);
                double d = nv | nf;
                normalVariation = std::min<double>(normalVariation, d);
            }
        }
    }
};

}
