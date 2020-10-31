#pragma once
#include <memory>
#include <vector>
#include <flann/flann.hpp>
#include <kt84/container_cast.h>

namespace kt84 {

    template <class TMeshBase, class TMesh>
    struct KdTree : public DerivedPtrHolder<TMesh, KdTree<TMeshBase, TMesh>> {
        static const int N = TMeshBase::Point::size_;
        
        void kdtree_build() {
            TMesh* mesh = DerivedPtrHolder<TMesh, KdTree<TMeshBase, TMesh>>::derived_ptr;
            int n_vertices = mesh->n_vertices();
            kdtree_points.resize(N * n_vertices);
            for (int i = 0; i < n_vertices; ++i) {
                auto& p = mesh->point(mesh->vertex_handle(i));
                for (int j = 0; j < N; ++j)
                    kdtree_points[N * i + j] = p[j];
            }
            kdtree = std::make_shared<flann::Index<flann::L2<double>>>(
                flann::Matrix<double>(kdtree_points.data(), n_vertices, N),
                flann::KDTreeIndexParams());
            kdtree->buildIndex();
        }
        std::vector<typename TMeshBase::VHandle> kdtree_search(const typename TMeshBase::Point& point, int knn, std::vector<double>& knn_dist) const {
            std::vector<int> knn_idx(knn);
            knn_dist.resize(knn);
            kdtree->knnSearch(
                flann::Matrix<double>(const_cast<double*>(point.data()), 1, N),
                flann::Matrix<int   >(                 knn_idx .data() , 1, knn),
                flann::Matrix<double>(                 knn_dist.data() , 1, knn),
                knn, flann::SearchParams());
            return container_cast<typename TMeshBase::VHandle>(knn_idx);
        }
        typename TMeshBase::VHandle kdtree_search(const typename TMeshBase::Point& point, double& dist) const {
            std::vector<double> knn_dist;
            auto knn_idx = kdtree_search(point, 1, knn_dist);
            dist = knn_dist.front();
            return knn_idx .front();
        }
        std::vector<typename TMeshBase::VHandle> kdtree_search(const typename TMeshBase::Point& point, int knn) const {
            std::vector<double> knn_dist;
            return kdtree_search(point, knn, knn_dist);
        }
        typename TMeshBase::VHandle kdtree_search(const typename TMeshBase::Point& point) const {
            double dist;
            return kdtree_search(point, dist);
        }
        
        std::shared_ptr<flann::Index<flann::L2<double>>> kdtree;
        std::vector<double> kdtree_points;
    };
}
