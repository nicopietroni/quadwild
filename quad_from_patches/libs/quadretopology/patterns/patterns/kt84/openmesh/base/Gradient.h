#pragma once
#include <functional>
#include <Eigen/Core>
#include "../vector_convert.h"
#include "../../eigen_util.h"
#include "../DerivedPtrHolder.h"

namespace kt84 {

struct Gradient_VertexTraits {
    double gradient_input;
};
template <int Dim>
struct Gradient_FaceTraits {
    Eigen::Matrix<double, Dim, 1> gradient_output;
};

template <class TMeshBase, class TMesh>
struct Gradient : public DerivedPtrHolder<TMesh, Gradient<TMeshBase, TMesh>> {
    void gradient_compute() {
        TMesh* mesh = DerivedPtrHolder<TMesh, Gradient<TMeshBase, TMesh>>::derived_ptr;
        
        for (auto f : mesh->faces()) {
            Eigen::Vector3d x[3];
            double          y[3];
            auto v = mesh->fv_iter(f);
            for (int i = 0; i < 3; ++i, ++v) {
                x[i] = o2e(mesh->point(*v));
                y[i] = mesh->data(*v).gradient_input;
            }
            mesh->data(f).gradient_output = eigen_util::compute_gradient(x[0], x[1], x[2], y[0], y[1], y[2]);
        }
    }
};

}
