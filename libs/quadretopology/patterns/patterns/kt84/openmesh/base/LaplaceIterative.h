#pragma once
#include <vector>
#include <Eigen/Core>
#include "../DerivedPtrHolder.h"

namespace kt84 {

template <int N>
struct LaplaceIterative_VertexTraits {
    struct Data {
        bool is_fixed;
        Eigen::Matrix<double, N, 1> value;
        Eigen::Matrix<double, N, 1> laplacian;
        
        Data()
            : is_fixed()
            , value    (Eigen::Matrix<double, N, 1>::Zero())
            , laplacian(Eigen::Matrix<double, N, 1>::Zero())
        {}
    } laplaceIterative;
};

struct LaplaceIterative_HalfedgeTraits {
    double laplaceIterative_weight;
    LaplaceIterative_HalfedgeTraits()
        : laplaceIterative_weight(1)
    {}
};

template <class TMeshBase, class TMesh, int N>
struct LaplaceIterative : public DerivedPtrHolder<TMesh, LaplaceIterative<TMeshBase, TMesh, N>> {
    typedef Eigen::Matrix<double, N, 1> Value;
    
    void laplaceIterative_compute(int num_iter = 1, double damping  = 0.25) {          // 0: no damping, 1: no change
        TMesh* mesh = get_mesh();
        
        int nv = mesh->n_vertices();
        std::vector<Value> values_new(nv, Value::Zero());
        
        for (int iter = 0; iter < num_iter; ++iter) {
            
            for (int i = 0; i < nv; ++i) {
                auto v = mesh->vertex_handle(i);
                auto& vdata = mesh->data(v).laplaceIterative;
                if (vdata.is_fixed) continue;
                
                // weighted sum among one-ring vertices
                Value value_sum = Value::Zero();
                double weight_sum = 0;
                for (auto h = mesh->voh_iter(v); h.is_valid(); ++h) {
                    double weight = mesh->data(*h).laplaceIterative_weight;
                    value_sum  += weight * mesh->data(mesh->to_vertex_handle(*h)).laplaceIterative.value;
                    weight_sum += weight;
                }
                
                // store new value to a temporary array
                value_sum /= weight_sum;
                values_new[i] = vdata.laplacian + value_sum;
            }
            
            // copy new value to the mesh
            for (int i = 0; i < nv; ++i) {
                auto v = mesh->vertex_handle(i);
                auto& vdata = mesh->data(v).laplaceIterative;
                if (vdata.is_fixed) continue;
                
                vdata.value = (1 - damping) * values_new[i] + damping * vdata.value;
            }
        }
    }
    void laplceIterative_set_laplacian_from_value() {
        TMesh* mesh = get_mesh();
        
        for (auto v : mesh->vertices()) {
            auto& vdata = mesh->data(v).laplaceIterative;
            
            vdata.laplacian = vdata.value;
            
            double weight_sum = 0;
            for (auto h = mesh->voh_iter(v); h; ++h)
                weight_sum += mesh->data(h).laplaceIterative_weight;
            
            for (auto h = mesh->voh_iter(v); h; ++h) {
                auto w = mesh->to_vertex_handle(h);
                auto& wdata = mesh->data(w).laplaceIterative;
                double weight = mesh->data(h).laplaceIterative_weight / weight_sum;
                
                vdata.laplacian -= weight * wdata.value;
            }
        }
    }
private:
    TMesh* get_mesh() const { return DerivedPtrHolder<TMesh, LaplaceIterative<TMeshBase, TMesh, N>>::derived_ptr; }
};

}
