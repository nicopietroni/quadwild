#ifndef QR_PATTERNS_H
#define QR_PATTERNS_H

#include <vector>

#include <Eigen/Core>

namespace QuadRetopology {
namespace internal {

template<class PolyMesh>
void computePattern(
        const Eigen::VectorXi &l,
        Eigen::MatrixXd& patchV,
        Eigen::MatrixXi& patchF,
        PolyMesh& patchMesh,
        std::vector<size_t>& borders,
        std::vector<size_t>& corners,
        std::vector<std::vector<size_t>>& sides);

}
}

#include "qr_patterns.cpp"

#endif // QR_PATTERNS_H
