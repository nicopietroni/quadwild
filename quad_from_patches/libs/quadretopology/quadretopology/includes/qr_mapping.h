#ifndef QR_MAPPING_H
#define QR_MAPPING_H

#include <vector>

#include <Eigen/Core>

namespace QuadRetopology {
namespace internal {

void computeQuadrangulation(
        const Eigen::MatrixXd& chartV,
        const Eigen::MatrixXi& chartF,
        const Eigen::MatrixXd& patchV,
        const Eigen::MatrixXi& patchF,
        const std::vector<std::vector<std::vector<size_t>>>& chartSideVertices,
        const std::vector<std::vector<double>>& chartSideLength,
        const std::vector<std::vector<size_t>>& chartSideSubdivision,
        const std::vector<std::vector<size_t>>& patchSides,
        Eigen::MatrixXd& uvMapV,
        Eigen::MatrixXi& uvMapF,
        Eigen::MatrixXd& quadrangulationV,
        Eigen::MatrixXi& quadrangulationF);

}

}


#endif // QR_MAPPING_H
