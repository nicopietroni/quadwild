#ifndef QR_ILP_H
#define QR_ILP_H

#include "qr_charts.h"
#include "qr_parameters.h"

#include <gurobi_c++.h>

#define ILP_FIND_SUBDIVISION -1
#define ILP_IGNORE -2

namespace QuadRetopology {

enum ILPStatus { SOLUTIONFOUND, SOLUTIONWRONG, INFEASIBLE };

namespace internal {

void solveILP(
        const ChartData& chartData,     
        const std::vector<double>& chartEdgeLength,
        const ILPMethod& method,
        const double alpha,
        const bool isometry,
        const bool regularityQuadrilaterals,
        const bool regularityNonQuadrilaterals,
        const double regularityNonQuadrilateralsWeight,
        const bool alignSingularities,
        const double alignSingularitiesWeight,
        const int repeatLosingConstraintsIterations,
        const bool repeatLosingConstraintsQuads,
        const bool repeatLosingConstraintsNonQuads,
        const bool repeatLosingConstraintsAlign,
        const bool feasibilityFix,
        const bool hardParityConstraint,
        const double timeLimit,
        const double gapLimit,
        const std::vector<float>& callbackTimeLimit,
        const std::vector<float>& callbackGapLimit,
        double& gap,
        ILPStatus& status,
        std::vector<int>& ilpResults);

}
}

#include "qr_ilp.cpp"

#endif // QR_ILP_H
