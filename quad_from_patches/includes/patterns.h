#ifndef QUADFROMPATCHES_PATTERNS_H
#define QUADFROMPATCHES_PATTERNS_H

#include <vector>

#include <Eigen/Core>

namespace qfp {

void computePattern(
        const Eigen::VectorXi &l,
        Eigen::MatrixXd& patchV,
        Eigen::MatrixXi& patchF,
        std::vector<size_t>& borders,
        std::vector<size_t>& corners);

}

#endif // QUADFROMPATCHES_PATTERNS_H
