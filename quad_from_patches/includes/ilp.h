#ifndef QUADFROMPATCHES_ILP_H
#define QUADFROMPATCHES_ILP_H

#include "charts.h"
#include "common.h"

namespace qfp {

enum ILPStatus { SOLUTIONFOUND, SOLUTIONWRONG, INFEASIBLE };

std::vector<int> solveILP(
        const ChartData& chartData,
        const std::vector<double>& edgeFactor,
        const double alpha,
        const ILPMethod& method,
        const bool regularity,
        const bool regularityForNonQuadrilaterals,
        const double nonQuadrilateralSimilarityFactor,
        const double timeLimit,
        double& gap,
        ILPStatus& status);

}

#include "ilp.cpp"

#endif // QUADFROMPATCHES_ILP_H
