#pragma once
#include <cmath>
#include <iostream>
#include "../DerivedPtrHolder.h"
#include "EdgeLength.h"

namespace kt84 {

struct FaceArea_FaceTraits {
    double faceArea;
    FaceArea_FaceTraits() : faceArea() {}
};
typedef EdgeLength_EdgeTraits FaceArea_EdgeTraits;

template <class TMeshBase, class TMesh>
struct FaceArea
    : public DerivedPtrHolder<TMesh, FaceArea<TMeshBase, TMesh>>
    , public EdgeLength<TMeshBase, TMesh>
{
    void faceArea_compute() {
        TMesh* mesh = DerivedPtrHolder<TMesh, FaceArea<TMeshBase, TMesh>>::derived_ptr;
        
        mesh->edgeLength_compute();
        
        for (auto f : mesh->faces()) {
            if (mesh->valence(f) != 3) {
                std::cerr << "Error: faceArea_compute() is called on polygonal face!\n";
                assert(false);
            }
            auto fh = mesh->fh_iter(f);
            double a = mesh->data(mesh->edge_handle(*fh++)).edgeLength;
            double b = mesh->data(mesh->edge_handle(*fh++)).edgeLength;
            double c = mesh->data(mesh->edge_handle(*fh++)).edgeLength;
            assert(!fh.is_valid());
            // Heron's formula
            double s = (a + b + c) / 2;
            mesh->data(f).faceArea = std::sqrt(s * (s - a) * (s - b) * (s - c));
        }
    }
};

}
