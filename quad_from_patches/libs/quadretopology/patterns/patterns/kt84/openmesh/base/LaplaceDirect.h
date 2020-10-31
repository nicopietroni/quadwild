#pragma once
#include <vector>
#include <memory>
#include <Eigen/Sparse>
#include "../DerivedPtrHolder.h"

namespace kt84 {

template <int N>
struct LaplaceDirect_VertexTraits {
    struct Data {
        bool is_fixed;
        Eigen::Matrix<double, N, 1> value;
        Eigen::Matrix<double, N, 1> laplacian;
        
        Data()
            : is_fixed()
            , value    (Eigen::Matrix<double, N, 1>::Zero())
            , laplacian(Eigen::Matrix<double, N, 1>::Zero())
        {}
    } laplaceDirect;
};

struct LaplaceDirect_HalfedgeTraits {
    double laplaceDirect_weight;
    LaplaceDirect_HalfedgeTraits()
        : laplaceDirect_weight(1)
    {}
};

template <class TMeshBase, class TMesh, int N>
struct LaplaceDirect : public DerivedPtrHolder<TMesh, LaplaceDirect<TMeshBase, TMesh, N>> {
    typedef Eigen::Matrix<double,  N, 1> Value ;
    typedef Eigen::Matrix<double, -1, N> Vector;
    typedef Eigen::SparseMatrix <double> Matrix;
    
    double laplaceDirect_constraintWeight;
    Matrix L;
    std::shared_ptr<Eigen::SimplicialCholesky<Matrix>> solver;
    
    LaplaceDirect()
        : laplaceDirect_constraintWeight(1000.0)            // better to set large value as default? not really sure...
    {}
    
    void laplaceDirect_factorize() {
        TMesh* mesh = get_mesh();
        
        int nv = mesh->n_vertices();
        
        typedef Eigen::Triplet<double> Triplet;
        std::vector<Triplet> triplets;
        for (int i = 0; i < nv; ++i) {
            auto v = mesh->vertex_handle(i);
            
            triplets.push_back(Triplet(i, i, 1));
            
            double weight_sum = 0;
            for (auto h = mesh->voh_iter(v); h.is_valid(); ++h)
                weight_sum += mesh->data(*h).laplaceDirect_weight;
            
            for (auto h = mesh->voh_iter(v); h.is_valid(); ++h) {
                auto w = mesh->to_vertex_handle(*h);
                double weight = mesh->data(*h).laplaceDirect_weight / weight_sum;
                
                triplets.push_back(Triplet(i, w.idx(), -weight));
            }
            
            if (mesh->data(v).laplaceDirect.is_fixed)
                triplets.push_back(Triplet(i, i, laplaceDirect_constraintWeight));
        }
        
        L.resize(nv, nv);
        L.setFromTriplets(triplets.begin(), triplets.end());
        
        solver.reset(new Eigen::SimplicialCholesky<Matrix>(L.transpose() * L));
    }
    
    void laplaceDirect_solve() {
        TMesh* mesh = get_mesh();
        
        int nv = mesh->n_vertices();
        
        if (!solver || L.rows() != nv)
            laplaceDirect_factorize();
        
        // set right hand side
        Vector b(nv, N);
        for (int i = 0; i < nv; ++i) {
            auto& vdata = mesh->data(mesh->vertex_handle(i)).laplaceDirect;
            
            b.row(i) = vdata.laplacian.transpose();
            
            if (vdata.is_fixed)
                b.row(i) += laplaceDirect_constraintWeight * vdata.value.transpose();
        }
        
        // solve!
        Vector x = solver->solve(L.transpose() * b);
        
        // copy result to Data::value
        for (int i = 0; i < nv; ++i) {
            auto& vdata = mesh->data(mesh->vertex_handle(i)).laplaceDirect;
            
            if (!vdata.is_fixed)
                vdata.value = x.row(i).transpose();
        }
    }
    
    void laplaceDirect_set_laplacian_from_value() {
        TMesh* mesh = get_mesh();
        
        for (auto v : mesh->vertices()) {
            auto& vdata = mesh->data(v).laplaceDirect;
            
            vdata.laplacian = vdata.value;
            
            double weight_sum = 0;
            for (auto h = mesh->voh_iter(v); h.is_valid(); ++h)
                weight_sum += mesh->data(*h).laplaceDirect_weight;
            
            for (auto h = mesh->voh_iter(v); h.is_valid(); ++h) {
                auto w = mesh->to_vertex_handle(*h);
                auto& wdata = mesh->data(w).laplaceDirect;
                double weight = mesh->data(*h).laplaceDirect_weight / weight_sum;
                
                vdata.laplacian -= weight * wdata.value;
            }
        }
    }
private:
    TMesh* get_mesh() const { return DerivedPtrHolder<TMesh, LaplaceDirect<TMeshBase, TMesh, N>>::derived_ptr; }
};

}
