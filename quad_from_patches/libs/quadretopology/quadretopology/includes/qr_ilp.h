#ifndef QR_ILP_H
#define QR_ILP_H

#include "qr_charts.h"
#include "qr_parameters.h"

namespace QuadRetopology {

enum ILPStatus { SOLUTIONFOUND, SOLUTIONWRONG, INFEASIBLE };

namespace internal {

inline std::vector<int> solveILP(
        const ChartData& chartData,
        const std::vector<size_t>& fixedSubsides,
        const std::vector<double>& chartEdgeLength,
        const ILPMethod& method,
        const double alpha,
        const bool isometry,
        const bool regularityForQuadrilaterals,
        const bool regularityForNonQuadrilaterals,
        const double regularityNonQuadrilateralWeight,
        const bool feasibilityFix,
        const bool hardParityConstraint,
        const double timeLimit,
        const double minimumGap,
        double& gap,
        ILPStatus& status);

}
}

#include "qr_ilp.cpp"

#endif // QR_ILP_H
