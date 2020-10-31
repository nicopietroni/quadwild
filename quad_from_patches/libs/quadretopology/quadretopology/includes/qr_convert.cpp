#include "qr_convert.h"

#include <vcg/complex/complex.h>

namespace QuadRetopology {
namespace internal {

template<class PolyMeshType>
void VCGToEigen(
        PolyMeshType& vcgMesh,
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F,
        std::vector<int>& vMap,
        std::vector<int>& fMap,
        bool selectedOnly,
        int numVerticesPerFace,
        int dim)
{
    assert(dim >= 2);
    assert(numVerticesPerFace > 2);

    int nSelectedVertices = 0;
    for (size_t i = 0; i < vcgMesh.vert.size(); i++){
        if ((!selectedOnly || vcgMesh.vert[i].IsS()) && !vcgMesh.vert[i].IsD()) {
            nSelectedVertices++;
        }
    }
    int nSelectedFaces = 0;
    for (size_t i = 0; i < vcgMesh.face.size(); i++){
        if ((!selectedOnly || vcgMesh.face[i].IsS()) && !vcgMesh.face[i].IsD()) {
            nSelectedFaces++;
        }
    }

    V.resize(nSelectedVertices, dim);
    F.resize(nSelectedFaces, numVerticesPerFace);

    vMap.resize(vcgMesh.vert.size(), -1);
    int vId = 0;
    for (size_t i = 0; i < vcgMesh.vert.size(); i++){
        if ((!selectedOnly || vcgMesh.vert[i].IsS()) && !vcgMesh.vert[i].IsD()) {
            vMap[i] = vId;
            for (int j = 0; j < dim; j++) {
                V(vId, j) = vcgMesh.vert[i].P()[j];
            }
            vId++;
        }
    }

    fMap.resize(vcgMesh.face.size(), -1);
    int fId = 0;
    for (size_t i = 0; i < vcgMesh.face.size(); i++){
        if ((!selectedOnly || vcgMesh.face[i].IsS()) && !vcgMesh.face[i].IsD()) {
            fMap[i] = fId;
            for (int j = 0; j < vcgMesh.face[i].VN(); j++) {
                F(fId, j) = vMap[vcg::tri::Index(vcgMesh, vcgMesh.face[i].V(j))];
            }
            fId++;
        }
    }
}

template<class PolyMeshType>
void eigenToVCG(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        PolyMeshType& vcgMesh,
        int numVertices,
        int dim)
{
    assert(dim >= 2);
    assert(numVertices > 2);

    vcgMesh.Clear();

    vcg::tri::Allocator<PolyMeshType>::AddVertices(vcgMesh, V.rows());
    for (int i = 0; i < V.rows(); i++) {
        typename PolyMeshType::CoordType vv(V(i,0), V(i,1), V(i,2));
        vcgMesh.vert[static_cast<size_t>(i)].P() = vv;
    }

    vcg::tri::Allocator<PolyMeshType>::AddFaces(vcgMesh, static_cast<size_t>(F.rows()));
    for (int i = 0; i < F.rows(); i++) {
        vcgMesh.face[static_cast<size_t>(i)].Alloc(numVertices);
        for (int j = 0; j < numVertices; j++) {
            vcgMesh.face[static_cast<size_t>(i)].V(j) = &(vcgMesh.vert[static_cast<size_t>(F(i,j))]);
        }
    }
}

}
}
