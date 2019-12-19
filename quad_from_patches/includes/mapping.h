#ifndef QUADFROMPATCHES_MAPPING_H
#define QUADFROMPATCHES_MAPPING_H

#include <vector>

#include <Eigen/Core>

namespace qfp {

std::vector<std::vector<size_t>> getPatchSides(
        Eigen::MatrixXd& patchV,
        Eigen::MatrixXi& patchF,
        std::vector<size_t>& borders,
        std::vector<size_t>& corners,
        const Eigen::VectorXi& l);

void computeQuadrangulation(
        const Eigen::MatrixXd& chartV,
        const Eigen::MatrixXi& chartF,
        const Eigen::MatrixXd& patchV,
        const Eigen::MatrixXi& patchF,        
        const std::vector<std::vector<size_t>>& chartSides,
        const std::vector<double>& chartSideLengths,
        const std::vector<std::vector<size_t>>& patchSides,
        Eigen::MatrixXd& uvMapV,
        Eigen::MatrixXi& uvMapF,
        Eigen::MatrixXd& quadrangulationV,
        Eigen::MatrixXi& quadrangulationF);

}


#endif // QUADFROMPATCHES_MAPPING_H
