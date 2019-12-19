#ifndef QUADFROMPATCHES_CONVERT_H
#define QUADFROMPATCHES_CONVERT_H

#include <Eigen/Core>
#include <vector>

namespace qfp {

template<class PolyMesh>
void VCGToEigen(
        PolyMesh& vcgMesh,
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F,
        std::vector<int>& vMap,
        std::vector<int>& fMap,
        bool selectedOnly = false,
        int numVerticesPerFace = 3,
        int dim = 3);

template<class PolyMesh>
void eigenToVCG(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        PolyMesh& vcgMesh,
        int numVertices = 3,
        int dim = 3);

}

#include "convert.cpp"

#endif // QUADFROMPATCHES_CONVERT_H
