#pragma once
#include <Eigen/Geometry>
#include "../vector_convert.h"
#include "../DerivedPtrHolder.h"

namespace kt84 {

template <class TMeshBase, class TMesh>
struct BoundingBox : public DerivedPtrHolder<TMesh, BoundingBox<TMeshBase, TMesh>> {
    static const int N = TMeshBase::Point::size_;
    Eigen::AlignedBox<double, N> boundingBox;

    void boundingBox_compute() {
        TMesh* mesh = DerivedPtrHolder<TMesh, BoundingBox<TMeshBase, TMesh>>::derived_ptr;
        boundingBox.setEmpty();
        for (auto v : mesh->vertices())
            boundingBox.extend(o2e(mesh->point(v)));
    }
    double boundingBox_diagonal_norm() const {
        return boundingBox.diagonal().norm();
    }
};

}
