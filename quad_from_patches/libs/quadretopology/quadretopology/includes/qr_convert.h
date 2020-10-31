#ifndef QR_CONVERT_H
#define QR_CONVERT_H

#include <Eigen/Core>
#include <vector>

namespace QuadRetopology {
namespace internal {

template<class PolyMeshType>
void VCGToEigen(
        PolyMeshType& vcgMesh,
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F,
        std::vector<int>& vMap,
        std::vector<int>& fMap,
        bool selectedOnly = false,
        int numVerticesPerFace = 3,
        int dim = 3);

template<class PolyMeshType>
void eigenToVCG(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        PolyMeshType& vcgMesh,
        int numVertices = 3,
        int dim = 3);

}
}

#include "qr_convert.cpp"

#endif // QR_CONVERT_H
