#pragma once
#include "../DerivedPtrHolder.h"
#include "FaceArea.h"

namespace kt84 {

struct CotanWeight_HalfedgeTraits {
    double cotangent;             // Cotangent of corner opposing this halfedge (i.e., the corner corresponds to this->next_halfedge->to_vertex)
    CotanWeight_HalfedgeTraits()
        : cotangent(0)
    {}
};
typedef FaceArea_EdgeTraits CotanWeight_EdgeTraits;
typedef FaceArea_FaceTraits CotanWeight_FaceTraits;

template <class TMeshBase, class TMesh>
struct CotanWeight
    : public DerivedPtrHolder<TMesh, CotanWeight<TMeshBase, TMesh>>
    , public FaceArea<TMeshBase, TMesh>
{
    double cotanWeight(OpenMesh::PolyConnectivity::EHandle e) {
        TMesh* mesh = get_mesh();
        return 0.5 * (mesh->data(mesh->halfedge_handle(e, 0)).cotangent + 
                      mesh->data(mesh->halfedge_handle(e, 1)).cotangent);
    }
    void cotanWeight_compute() {
        TMesh* mesh = get_mesh();
        mesh->faceArea_compute();
        
        for (auto f : mesh->faces()) {
            // Intrinsig cotangent formula by Alec Jacobson: http://www.alecjacobson.com/weblog/?p=2405
            OpenMesh::PolyConnectivity::HHandle h[3];
            double l[3];
            int i = 0;
            for (auto fh = mesh->fh_iter(f); fh.is_valid(); ++fh, ++i) {
                h[i] = *fh;
                l[i] = mesh->data(mesh->edge_handle(h[i])).edgeLength;
            }
            double A = mesh->data(f).faceArea;
            for (i = 0; i < 3; ++i) {
                int i1 = (i + 1) % 3;
                int i2 = (i + 2) % 3;
                mesh->data(h[i]).cotangent = (l[i1] * l[i1] + l[i2] * l[i2] - l[i] * l[i]) / (4 * A);
            }
        }
    }
private:
    TMesh* get_mesh() const { return DerivedPtrHolder<TMesh, CotanWeight<TMeshBase, TMesh>>::derived_ptr; }
};

}
