#pragma once

/*
    Implementation of the following paper:
        Chen-shi Dong, Guo-zhao Wang
        Curvatures estimation on triangular mesh
        Journal of Zhejiang University Science
        August 2005, Volume 6, Issue 1 Supplement, pp 128-136
        DOI: 10.1007/BF02887228
        http://www.zju.edu.cn/jzus/2005/A05S1/A05S121.pdf
    
    which I found at the following forum posting:
        http://www.sidefx.com/index.php?option=com_forum&Itemid=172&page=viewtopic&t=13855&highlight=&sid=2c56b6ede4383079db2ce0b1d9f97be5
*/
#include "../DerivedPtrHolder.h"

namespace kt84 {

struct NormalCurvature_VertexTraits {
    double normalCurvature;
    NormalCurvature_VertexTraits()
        : normalCurvature()
    {}
};

template <class TMeshBase, class TMesh>
struct NormalCurvature : public DerivedPtrHolder<TMesh, NormalCurvature<TMeshBase, TMesh>> {
    void normalCurvature_compute() {
        TMesh* mesh = DerivedPtrHolder<TMesh, NormalCurvature<TMeshBase, TMesh>>::derived_ptr;
        
        for (auto v0 : mesh->vertices()) {
            double& k0 = mesh->data(v0).normalCurvature;
            k0 = 0;
            auto n0 = mesh->normal(v0);
            
            for (auto v1 = mesh->vv_iter(v0); v1.is_valid(); ++v1) {
                auto n1 = mesh->normal(*v1);
                auto dp = mesh->point(*v1) - mesh->point(v0);
                auto dn = n1 - n0;
                double k1 = (dp | dn) / (dp | dp);
                if (std::abs(k0) < std::abs(k1))
                    k0 = k1;
            }
        }
    }
};

}
