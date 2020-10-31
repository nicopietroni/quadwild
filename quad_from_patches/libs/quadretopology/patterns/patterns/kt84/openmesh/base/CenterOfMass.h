#pragma once
#include "../vector_convert.h"
#include "../DerivedPtrHolder.h"

namespace kt84 {

template <class TMeshBase, class TMesh>
struct CenterOfMass : public DerivedPtrHolder<TMesh, CenterOfMass<TMeshBase, TMesh>> {
    static const int N = TMeshBase::Point::size_;
    Eigen::Matrix<double, N, 1> centerOfMass;
    
    void centerOfMass_compute() {
        TMesh* mesh = DerivedPtrHolder<TMesh, CenterOfMass<TMeshBase, TMesh>>::derived_ptr;
        centerOfMass.setZero();
        for (auto v : mesh->vertices())
            centerOfMass += o2e(mesh->point(v));
        centerOfMass /= mesh->n_vertices();
    }
};

}
