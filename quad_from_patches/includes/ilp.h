#ifndef QUADFROMPATCHES_ILP_H
#define QUADFROMPATCHES_ILP_H

#include "charts.h"
#include "common.h"

namespace qfp {

enum ILPStatus { SOLUTIONFOUND, SOLUTIONWRONG, INFEASIBLE };

std::vector<int> solveILP(
        const ChartData& chartData,
        const std::vector<double>& edgeFactor,
        const ILPMethod& method,
        const double alpha,
        const bool isometry,
        const bool regularityForQuadrilaterals,
        const bool regularityForNonQuadrilaterals,
        const double nonQuadrilateralSimilarityFactor,
        const bool hardParityConstraint,
        const double timeLimit,
        const double minimumGap,
        double& gap,
        ILPStatus& status);

}

#include "ilp.cpp"

#endif // QUADFROMPATCHES_ILP_H
