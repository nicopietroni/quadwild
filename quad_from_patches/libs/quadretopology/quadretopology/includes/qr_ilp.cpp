#include "qr_ilp.h"

#include <gurobi_c++.h>

#define MIN_SUBDIVISION_VALUE 1
#define FEASIBILITY_FIX_COST 1000000.0

namespace QuadRetopology {
namespace internal {

class GurobiCallBack : public GRBCallback
{
    const std::vector<float>& times;
    const std::vector<float>& gaps;
    size_t index;
public:
    GurobiCallBack(const std::vector<float>& times, const std::vector<float>& gaps);
protected:
    void callback() override;
};

void getChartSubsideSum(
        const ChartData& chartData,
        const size_t& cId,
        const std::vector<GRBVar>& vars,
        const std::vector<bool>& isFixed,
        const std::vector<int>& ilpResults,
        const bool hardParityConstraint,
        std::vector<GRBLinExpr>& chartSubsideSum);

void getChartSubsideSumResults(
        const ChartData& chartData,
        const size_t& cId,
        const std::vector<int>& results,
        std::vector<int>& chartSubsideSum);

GRBQuadExpr getGurobiCostTermInteger(GRBModel& model, const ILPMethod& method, const GRBLinExpr& value);
GRBQuadExpr getGurobiCostTermContinuous(GRBModel& model, const ILPMethod& method, const GRBLinExpr& value);

GRBLinExpr getGurobiAbsInteger(GRBModel& model, const GRBLinExpr& value);
GRBLinExpr getGurobiAbsContinuous(GRBModel& model, const GRBLinExpr& value);

inline void solveILP(
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
        std::vector<int>& ilpResults)
{
    using namespace std;

    vector<int> result;

    std::vector<bool> constraintRespected;

    bool allRespectConstraints = false;
    int it = 0;
    do {
        std::cout << std::endl << std::endl << " --- OPTIMIZATION: iteration " << it << " ----" << std::endl;
        try {
            vector<bool> isFixed(chartData.subsides.size(), false);
            vector<bool> isComputable(chartData.subsides.size(), true);
            for (size_t subsideId = 0; subsideId < chartData.subsides.size(); ++subsideId) {
                if (ilpResults[subsideId] == ILP_IGNORE) {
                    isComputable[subsideId] = false;
                }
                else if (ilpResults[subsideId] >= 0) {
                    isFixed[subsideId] = true;
                }
            }

            result = ilpResults;

            GRBEnv env = GRBEnv();

            GRBModel model = GRBModel(env);

            GurobiCallBack cb(callbackTimeLimit, callbackGapLimit);
            if (!callbackTimeLimit.empty() && !callbackGapLimit.empty() && callbackTimeLimit.size() == callbackGapLimit.size()) {
                model.setCallback(&cb);
            }

#ifdef GUROBI_NON_VERBOSE
            model.set(GRB_IntParam_OutputFlag, 0);
#endif

            model.set(GRB_DoubleParam_TimeLimit, timeLimit);
            model.set(GRB_DoubleParam_MIPGap, gapLimit);

            // Create variables
            GRBQuadExpr obj = 0;
            GRBQuadExpr supportObj = 0;

            vector<GRBVar> vars(chartData.subsides.size());

            size_t constraintRespectedId = 0;

            for (size_t subsideId = 0; subsideId < chartData.subsides.size(); subsideId++) {
                if (!isComputable[subsideId]) {
                    continue;
                }

                const ChartSubside& subside = chartData.subsides[subsideId];

                //If it is not a border (free)
                if (!isFixed[subsideId]) {
                    vars[subsideId] = model.addVar(MIN_SUBDIVISION_VALUE, GRB_INFINITY, 0.0, GRB_INTEGER, "s" + to_string(subsideId));
                }
                else if (feasibilityFix && subside.isOnBorder) {
                    const int fixedSize = hardParityConstraint ? ilpResults[subsideId] : std::round(ilpResults[subsideId] / 2.0);

                    vars[subsideId] = model.addVar(std::max(fixedSize - 1, MIN_SUBDIVISION_VALUE), fixedSize + 1, 0.0, GRB_INTEGER, "s" + to_string(subsideId));
                    GRBLinExpr value = vars[subsideId] - fixedSize;

                    supportObj += getGurobiAbsInteger(model, value) * FEASIBILITY_FIX_COST * (1.0 / fixedSize);

                    isFixed[subsideId] = false;
                }
            }

            std::cout << chartData.subsides.size() << " subsides!" << std::endl;

            for (size_t cId = 0; cId < chartData.charts.size(); cId++) {
                const Chart& chart = chartData.charts[cId];

                bool computable = true;

                for (size_t i = 0; i < chart.chartSubsides.size(); i++) {
                    const size_t subsideId = chart.chartSubsides[i];
                    if (!isComputable[subsideId])
                        computable = false;
                }

                if (computable && chart.faces.size() > 0) {
                    size_t nSides = chart.chartSides.size();


                    size_t numRegularityTerms = 0;
                    size_t numIsometryTerms = 0;

                    GRBQuadExpr regExpr = 0;
                    GRBQuadExpr isoExpr = 0;

                    /* ------------------------ ISOMETRY ------------------------ */

                    if (isometry) {
                        for (size_t i = 0; i < chart.chartSubsides.size(); i++) {
                            const size_t subsideId = chart.chartSubsides[i];
                            const ChartSubside& subside = chartData.subsides[subsideId];

                            //If it is not fixed (free)
                            if (!isFixed[subsideId]) {
                                double edgeLength = chartEdgeLength[cId];

                                double sideSubdivision = subside.length / edgeLength;
                                if (!hardParityConstraint) {
                                    sideSubdivision /= 2.0;
                                }

                                sideSubdivision = std::max(static_cast<double>(MIN_SUBDIVISION_VALUE), sideSubdivision);

                                GRBLinExpr value = vars[subsideId] - sideSubdivision;

                                isoExpr += getGurobiCostTermContinuous(model, method, value);
                                numIsometryTerms++;
                            }
                        }
                    }


                    /* ------------------------ REGULARITY ------------------------ */

                    std::vector<GRBLinExpr> chartSubsideSum;
                    getChartSubsideSum(chartData, cId, vars, isFixed, ilpResults, hardParityConstraint, chartSubsideSum);

                    for (size_t j = 0; j < nSides; j++) {
                        GRBLinExpr value = 0.0;
                        bool valueComputed = false;

                        //Regularity for quad case
                        if (nSides == 4 && regularityQuadrilaterals) {
                            const GRBLinExpr& subside0Sum = chartSubsideSum[j];
                            const GRBLinExpr& subside2Sum = chartSubsideSum[(j+2)%nSides];

                            value = subside0Sum - subside2Sum;
                            valueComputed = true;
                        }
                        //Regularity for triangular case
                        else if (nSides == 3 && regularityNonQuadrilaterals) {
                            const GRBLinExpr& subside0Sum = chartSubsideSum[j];
                            const GRBLinExpr& subside1Sum = chartSubsideSum[(j+1)%nSides];
                            const GRBLinExpr& subside2Sum = chartSubsideSum[(j+2)%nSides];

                            GRBVar c = model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER);
                            model.addConstr(subside0Sum + 1 <= subside1Sum + subside2Sum + c);

                            value = c;
                            valueComputed = true;
                        }
                        //Regularity for pentagonal case
                        else if (nSides == 5 && regularityNonQuadrilaterals) {
                            const GRBLinExpr& subside0Sum = chartSubsideSum[j];
                            const GRBLinExpr& subside1Sum = chartSubsideSum[(j+1)%nSides];
                            const GRBLinExpr& subside2Sum = chartSubsideSum[(j+2)%nSides];
                            const GRBLinExpr& subside3Sum = chartSubsideSum[(j+3)%nSides];
                            const GRBLinExpr& subside4Sum = chartSubsideSum[(j+4)%nSides];

                            GRBVar c = model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER);
                            model.addConstr(subside0Sum + subside1Sum + 1 <= subside2Sum + subside3Sum + subside4Sum + c);

                            value = c;
                            valueComputed = true;
                        }
                        //Regularity for hexagonal case
                        else if (nSides == 6 && regularityNonQuadrilaterals) {
                            const GRBLinExpr& subside0Sum = chartSubsideSum[j];
                            const GRBLinExpr& subside2Sum = chartSubsideSum[(j+2)%nSides];
                            const GRBLinExpr& subside4Sum = chartSubsideSum[(j+4)%nSides];

                            GRBVar c = model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER);
                            model.addConstr(subside0Sum + 1 <= subside2Sum + subside4Sum + c);

                            GRBLinExpr parityEquation = subside0Sum + subside2Sum + subside4Sum;
                            GRBVar hexParity = model.addVar(3, GRB_INFINITY, 0.0, GRB_INTEGER);
                            GRBVar hexFree = model.addVar(0, 1, 0.0, GRB_INTEGER);
                            model.addConstr(hexParity * 2 == parityEquation + hexFree);

                            value = (c + hexFree) / 2.0;
                            valueComputed = true;
                        }

                        if (valueComputed) {
                            if (constraintRespected.empty() || constraintRespected[constraintRespectedId]) {
                                if (nSides == 4) {
                                    regExpr += getGurobiCostTermInteger(model, ILPMethod::ABS, value) / nSides;
                                }
                                else {
                                    regExpr += getGurobiCostTermInteger(model, ILPMethod::ABS, value) / nSides * regularityNonQuadrilateralsWeight;
                                }

                                numRegularityTerms++;
                            }

                            constraintRespectedId++;
                        }

                    }



                    /* ------------------------ ALIGN SINGULARITIES ------------------------ */

                    if (alignSingularities && (nSides == 3 || nSides == 5 || nSides == 6)) {
                        for (size_t j = 0; j < nSides; j++) {
                            GRBLinExpr value1 = 0.0;
                            GRBLinExpr value2 = 0.0;

                            int currentChartId = cId;
                            int currentChartSideId = j;
                            size_t currentChartNSides = chartData.charts[currentChartId].chartSides.size();

                            bool valueComputed = false;
                            do {
                                const Chart& currentChart = chartData.charts[currentChartId];
                                const ChartSide& currentChartSide = currentChart.chartSides[currentChartSideId];

                                if (currentChartSide.subsides.size() != 1) {
                                    currentChartId = -1;
                                    currentChartSideId = -1;
                                    currentChartNSides = 0;
                                }
                                else {
                                    int adjOppositeSideId = currentChartSideId;
                                    if (currentChartId != static_cast<int>(cId)) {
                                        adjOppositeSideId = (currentChartSideId + 2) % chartData.charts[currentChartId].chartSides.size();
                                    }
                                    const ChartSide& adjOppositeSide = currentChart.chartSides[adjOppositeSideId];

                                    const std::array<int, 2>& incidentCharts = chartData.subsides[adjOppositeSide.subsides[0]].incidentCharts;
                                    const std::array<int, 2>& incidentChartSides = chartData.subsides[adjOppositeSide.subsides[0]].incidentChartSideId;

                                    if (incidentCharts[0] == static_cast<int>(currentChartId)) {
                                        currentChartId = incidentCharts[1];
                                        currentChartSideId = incidentChartSides[1];
                                    }
                                    else if (incidentCharts[1] == static_cast<int>(currentChartId)) {
                                        currentChartId = incidentCharts[0];
                                        currentChartSideId = incidentChartSides[0];
                                    }

                                    if (currentChartId > -1) {
                                        currentChartNSides = chartData.charts[currentChartId].chartSides.size();
                                    }
                                }

                            } while (currentChartId != static_cast<int>(cId) && currentChartId > -1 && currentChartNSides == 4);

                            if (currentChartId != static_cast<int>(cId) && currentChartId > -1 && (currentChartNSides == 3 || currentChartNSides == 5 || currentChartNSides == 6)) {
                                bool currentComputable = true;
                                for (size_t i = 0; i < chartData.charts[currentChartId].chartSubsides.size(); i++) {
                                    const size_t subsideId = chartData.charts[currentChartId].chartSubsides[i];
                                    if (!isComputable[subsideId])
                                        currentComputable = false;
                                }

                                if (currentComputable) {
                                    GRBLinExpr valueDown1 = 0.0;
                                    GRBLinExpr valueUp1 = 0.0;
                                    GRBLinExpr valueDown2 = 0.0;
                                    GRBLinExpr valueUp2 = 0.0;

                                    //Singularity alignment for triangular case
                                    if (nSides == 3) {
                                        const GRBLinExpr& subside0Sum = chartSubsideSum[j];
                                        const GRBLinExpr& subside1Sum = chartSubsideSum[(j+1)%nSides];
                                        const GRBLinExpr& subside2Sum = chartSubsideSum[(j+2)%nSides];

                                        GRBLinExpr down = subside0Sum + subside2Sum - subside1Sum;
                                        GRBLinExpr up = subside1Sum + subside0Sum - subside2Sum;

                                        valueDown1 = down;
                                        valueUp1 = up;
                                    }
                                    //Singularity alignment for pentagonal case
                                    else if (nSides == 5) {
                                        const GRBLinExpr& subside0Sum = chartSubsideSum[j];
                                        const GRBLinExpr& subside1Sum = chartSubsideSum[(j+1)%nSides];
                                        const GRBLinExpr& subside2Sum = chartSubsideSum[(j+2)%nSides];
                                        const GRBLinExpr& subside3Sum = chartSubsideSum[(j+3)%nSides];
                                        const GRBLinExpr& subside4Sum = chartSubsideSum[(j+4)%nSides];

                                        GRBLinExpr down = (subside0Sum + subside1Sum + subside2Sum) - (subside3Sum + subside4Sum);
                                        GRBLinExpr up = (subside3Sum + subside4Sum + subside0Sum) - (subside1Sum + subside2Sum);

                                        valueDown1 = down;
                                        valueUp1 = up;
                                    }
                                    //Singularity alignment for hexagonal case
                                    else if (nSides == 6) {
                                        const GRBLinExpr& subside0Sum = chartSubsideSum[j];
                                        const GRBLinExpr& subside2Sum = chartSubsideSum[(j+2)%nSides];
                                        const GRBLinExpr& subside4Sum = chartSubsideSum[(j+4)%nSides];

                                        GRBLinExpr down = subside0Sum + subside2Sum - subside4Sum;
                                        GRBLinExpr up = subside4Sum + subside0Sum - subside2Sum;

                                        valueDown1 = down;
                                        valueUp1 = up;
                                    }




                                    std::vector<GRBLinExpr> adjChartSubsideSum;
                                    getChartSubsideSum(chartData, currentChartId, vars, isFixed, ilpResults, hardParityConstraint, adjChartSubsideSum);

                                    //Singularity alignment for triangular case
                                    if (currentChartNSides == 3) {
                                        const GRBLinExpr& subside0Sum = adjChartSubsideSum[currentChartSideId];
                                        const GRBLinExpr& subside1Sum = adjChartSubsideSum[(currentChartSideId+1)%currentChartNSides];
                                        const GRBLinExpr& subside2Sum = adjChartSubsideSum[(currentChartSideId+2)%currentChartNSides];

                                        GRBLinExpr down = subside0Sum + subside2Sum - subside1Sum;
                                        GRBLinExpr up = subside1Sum + subside0Sum - subside2Sum;

                                        valueDown2 = down;
                                        valueUp2 = up;
                                    }
                                    //Singularity alignment for pentagonal case
                                    else if (currentChartNSides == 5) {
                                        const GRBLinExpr& subside0Sum = adjChartSubsideSum[currentChartSideId];
                                        const GRBLinExpr& subside1Sum = adjChartSubsideSum[(currentChartSideId+1)%currentChartNSides];
                                        const GRBLinExpr& subside2Sum = adjChartSubsideSum[(currentChartSideId+2)%currentChartNSides];
                                        const GRBLinExpr& subside3Sum = adjChartSubsideSum[(currentChartSideId+3)%currentChartNSides];
                                        const GRBLinExpr& subside4Sum = adjChartSubsideSum[(currentChartSideId+4)%currentChartNSides];

                                        GRBLinExpr down = (subside0Sum + subside1Sum + subside2Sum) - (subside3Sum + subside4Sum);
                                        GRBLinExpr up = (subside3Sum + subside4Sum + subside0Sum) - (subside1Sum + subside2Sum);

                                        valueDown2 = down;
                                        valueUp2 = up;
                                    }
                                    //Singularity alignment for hexagonal case
                                    else if (currentChartNSides == 6) {
                                        const GRBLinExpr& subside0Sum = adjChartSubsideSum[currentChartSideId];
                                        const GRBLinExpr& subside2Sum = adjChartSubsideSum[(currentChartSideId+2)%currentChartNSides];
                                        const GRBLinExpr& subside4Sum = adjChartSubsideSum[(currentChartSideId+4)%currentChartNSides];

                                        GRBLinExpr down = subside0Sum + subside2Sum - subside4Sum;
                                        GRBLinExpr up = subside4Sum + subside0Sum - subside2Sum;

                                        valueDown2 = down;
                                        valueUp2 = up;
                                    }

                                    value1 = valueUp1 - valueDown2;
                                    value2 = valueDown1 - valueUp2;
                                    valueComputed = true;
                                }
                            }

                            if (valueComputed) {
                                if (constraintRespected.empty() || constraintRespected[constraintRespectedId]) {
                                    regExpr += (
                                        0.5 * getGurobiCostTermInteger(model, ILPMethod::ABS, value1) +
                                        0.5 * getGurobiCostTermInteger(model, ILPMethod::ABS, value2)
                                        )  / nSides * alignSingularitiesWeight;
                                }

                                constraintRespectedId++;
                            }
                        }
                    }

                    const double isoWeight = alpha;
                    const double regWeight = (1 - alpha);

                    if (numIsometryTerms > 0)
                        obj += isoWeight * isoExpr / numIsometryTerms;
                    if (numRegularityTerms > 0)
                        obj += regWeight * regExpr / numRegularityTerms;


                    //Even side size sum constraint in a chart
                    if (hardParityConstraint) {
                        GRBLinExpr sumExp = 0;
                        for (const size_t& subsideId : chart.chartSubsides) {
                            const ChartSubside& subside = chartData.subsides[subsideId];

                            if (!isFixed[subsideId]) {
                                assert(isComputable[subsideId] && !isFixed[subsideId]);
                                sumExp += vars[subsideId];
                            }
                            else {                                
                                const int fixedSize = hardParityConstraint ? ilpResults[subsideId] : std::round(ilpResults[subsideId] / 2.0);
                                sumExp += fixedSize;
                            }
                        }

                        GRBVar free = model.addVar(2, GRB_INFINITY, 0.0, GRB_INTEGER);
                        model.addConstr(free * 2 == sumExp);
                    }

                    if (chart.chartSides.size() < 3 || chart.chartSides.size() > 6) {
                        std::cout << "Chart " << cId << " has " << chart.chartSides.size() << " sides." << std::endl;
#ifdef ASSERT_FOR_NUMBER_SIDES
                        assert(chart.chartSides.size() >= 3 && chart.chartSides.size() <= 6);
#endif
                    }

                }
            }

            //Set objective function
            model.setObjective(obj + supportObj, GRB_MINIMIZE);

            //model.write("out.lp");

            //Optimize model
            model.optimize();

            for (size_t subsideId = 0; subsideId < chartData.subsides.size(); subsideId++) {
                const ChartSubside& subside = chartData.subsides[subsideId];
                if (isComputable[subsideId]) {
                    if (!isFixed[subsideId]) {
                        assert(isComputable[subsideId] && !isFixed[subsideId]);
                        result[subsideId] = static_cast<int>(std::round(vars[subsideId].get(GRB_DoubleAttr_X)));
                    }
                    else {
                        const int fixedSize = hardParityConstraint ? ilpResults[subsideId] : std::round(ilpResults[subsideId] / 2.0);
                        result[subsideId] = fixedSize;
                    }
    
                    if (!hardParityConstraint) {
                        result[subsideId] *= 2;
                    }
                }
            }


            if (it < repeatLosingConstraintsIterations) {
                int numLostConstraintsQuad = 0;
                int numLostConstraintsNonQuad = 0;
                int numLostConstraintsAlign = 0;

                if (constraintRespected.empty()) {
                    constraintRespected.resize(constraintRespectedId, true);
                }

                allRespectConstraints = true;

                constraintRespectedId = 0;

                for (size_t cId = 0; cId < chartData.charts.size(); cId++) {
                    const Chart& chart = chartData.charts[cId];

                    bool computable = true;

                    for (size_t i = 0; i < chart.chartSubsides.size(); i++) {
                        const size_t subsideId = chart.chartSubsides[i];
                        if (!isComputable[subsideId])
                            computable = false;
                    }

                    if (computable && chart.faces.size() > 0) {
                        size_t nSides = chart.chartSides.size();

                        std::vector<int> chartSubsideSumResults;
                        getChartSubsideSumResults(chartData, cId, result, chartSubsideSumResults);

                        for (size_t j = 0; j < nSides; j++) {
                            int value = 0.0;
                            bool valueComputed = false;
                            bool respected = false;

                            //Regularity for quad case
                            if (nSides == 4 && regularityQuadrilaterals) {
                                const int& subside0Sum = chartSubsideSumResults[j];
                                const int& subside2Sum = chartSubsideSumResults[(j+2)%nSides];

                                value = subside0Sum - subside2Sum;
                                respected = value == 0;
                                valueComputed = true;
                            }
                            //Regularity for triangular case
                            else if (nSides == 3 && regularityNonQuadrilaterals) {
                                const int& subside0Sum = chartSubsideSumResults[j];
                                const int& subside1Sum = chartSubsideSumResults[(j+1)%nSides];
                                const int& subside2Sum = chartSubsideSumResults[(j+2)%nSides];

                                respected = subside0Sum < subside1Sum + subside2Sum;
                                valueComputed = true;
                            }
                            //Regularity for pentagonal case
                            else if (nSides == 5 && regularityNonQuadrilaterals) {
                                const int& subside0Sum = chartSubsideSumResults[j];
                                const int& subside1Sum = chartSubsideSumResults[(j+1)%nSides];
                                const int& subside2Sum = chartSubsideSumResults[(j+2)%nSides];
                                const int& subside3Sum = chartSubsideSumResults[(j+3)%nSides];
                                const int& subside4Sum = chartSubsideSumResults[(j+4)%nSides];

                                respected = subside0Sum + subside1Sum < subside2Sum + subside3Sum + subside4Sum;
                                valueComputed = true;
                            }
                            //Regularity for hexagonal case
                            else if (nSides == 6 && regularityNonQuadrilaterals) {
                                const int& subside0Sum = chartSubsideSumResults[j];
                                const int& subside2Sum = chartSubsideSumResults[(j+2)%nSides];
                                const int& subside4Sum = chartSubsideSumResults[(j+4)%nSides];

                                respected = subside0Sum < subside2Sum + subside4Sum && ((subside0Sum + subside2Sum + subside4Sum) % 2 == 0);
                                valueComputed = true;
                            }

                            if (valueComputed) {
                                if (!respected) {
                                    if (nSides == 4 && repeatLosingConstraintsQuads) {
                                        numLostConstraintsQuad++;
                                        allRespectConstraints = false;
                                        constraintRespected[constraintRespectedId] = false;
                                    }
                                    else if (((nSides == 3 || nSides == 5 || nSides == 6) && repeatLosingConstraintsNonQuads)) {
                                        numLostConstraintsNonQuad++;
                                        allRespectConstraints = false;
                                        constraintRespected[constraintRespectedId] = false;
                                    }
                                }

                                constraintRespectedId++;
                            }
                        }

                        if (alignSingularities && (nSides == 3 || nSides == 5 || nSides == 6)) {
                            for (size_t j = 0; j < nSides; j++) {
                                int value1 = 0.0;
                                int value2 = 0.0;

                                int currentChartId = cId;
                                int currentChartSideId = j;
                                size_t currentChartNSides = chartData.charts[currentChartId].chartSides.size();

                                bool valueComputed = false;
                                bool respected = false;
                                do {
                                    const Chart& currentChart = chartData.charts[currentChartId];
                                    const ChartSide& currentChartSide = currentChart.chartSides[currentChartSideId];

                                    if (currentChartSide.subsides.size() != 1) {
                                        currentChartId = -1;
                                        currentChartSideId = -1;
                                        currentChartNSides = 0;
                                    }
                                    else {
                                        int adjOppositeSideId = currentChartSideId;
                                        if (currentChartId != static_cast<int>(cId)) {
                                            adjOppositeSideId = (currentChartSideId + 2) % chartData.charts[currentChartId].chartSides.size();
                                        }
                                        const ChartSide& adjOppositeSide = currentChart.chartSides[adjOppositeSideId];

                                        const std::array<int, 2>& incidentCharts = chartData.subsides[adjOppositeSide.subsides[0]].incidentCharts;
                                        const std::array<int, 2>& incidentChartSides = chartData.subsides[adjOppositeSide.subsides[0]].incidentChartSideId;

                                        if (incidentCharts[0] == static_cast<int>(currentChartId)) {
                                            currentChartId = incidentCharts[1];
                                            currentChartSideId = incidentChartSides[1];
                                        }
                                        else if (incidentCharts[1] == static_cast<int>(currentChartId)) {
                                            currentChartId = incidentCharts[0];
                                            currentChartSideId = incidentChartSides[0];
                                        }

                                        if (currentChartId > -1) {
                                            currentChartNSides = chartData.charts[currentChartId].chartSides.size();
                                        }
                                    }

                                } while (currentChartId != static_cast<int>(cId) && currentChartId > -1 && currentChartNSides == 4);

                                if (currentChartId != static_cast<int>(cId) && currentChartId > -1 && (currentChartNSides == 3 || currentChartNSides == 5 || currentChartNSides == 6)) {
                                    bool currentComputable = true;
                                    for (size_t i = 0; i < chartData.charts[currentChartId].chartSubsides.size(); i++) {
                                        const size_t subsideId = chartData.charts[currentChartId].chartSubsides[i];
                                        if (!isComputable[subsideId])
                                            currentComputable = false;
                                    }

                                    if (currentComputable) {

                                        int valueDown1 = 0.0;
                                        int valueUp1 = 0.0;
                                        int valueDown2 = 0.0;
                                        int valueUp2 = 0.0;

                                        //Singularity alignment for triangular case
                                        if (nSides == 3) {
                                            const int& subside0Sum = chartSubsideSumResults[j];
                                            const int& subside1Sum = chartSubsideSumResults[(j+1)%nSides];
                                            const int& subside2Sum = chartSubsideSumResults[(j+2)%nSides];

                                            int down = subside0Sum + subside2Sum - subside1Sum;
                                            int up = subside1Sum + subside0Sum - subside2Sum;

                                            valueDown1 = down;
                                            valueUp1 = up;
                                        }
                                        //Singularity alignment for pentagonal case
                                        else if (nSides == 5) {
                                            const int& subside0Sum = chartSubsideSumResults[j];
                                            const int& subside1Sum = chartSubsideSumResults[(j+1)%nSides];
                                            const int& subside2Sum = chartSubsideSumResults[(j+2)%nSides];
                                            const int& subside3Sum = chartSubsideSumResults[(j+3)%nSides];
                                            const int& subside4Sum = chartSubsideSumResults[(j+4)%nSides];

                                            int down = (subside0Sum + subside1Sum + subside2Sum) - (subside3Sum + subside4Sum);
                                            int up = (subside3Sum + subside4Sum + subside0Sum) - (subside1Sum + subside2Sum);

                                            valueDown1 = down;
                                            valueUp1 = up;
                                        }
                                        //Singularity alignment for hexagonal case
                                        else if (nSides == 6) {
                                            const int& subside0Sum = chartSubsideSumResults[j];
                                            const int& subside2Sum = chartSubsideSumResults[(j+2)%nSides];
                                            const int& subside4Sum = chartSubsideSumResults[(j+4)%nSides];

                                            int down = subside0Sum + subside2Sum - subside4Sum;
                                            int up = subside4Sum + subside0Sum - subside2Sum;

                                            valueDown1 = down;
                                            valueUp1 = up;
                                        }




                                        std::vector<int> adjChartSubsideSumResults;
                                        getChartSubsideSumResults(chartData, currentChartId, result, adjChartSubsideSumResults);

                                        //Singularity alignment for triangular case
                                        if (currentChartNSides == 3) {
                                            const int& subside0Sum = adjChartSubsideSumResults[currentChartSideId];
                                            const int& subside1Sum = adjChartSubsideSumResults[(currentChartSideId+1)%currentChartNSides];
                                            const int& subside2Sum = adjChartSubsideSumResults[(currentChartSideId+2)%currentChartNSides];

                                            int down = subside0Sum + subside2Sum - subside1Sum;
                                            int up = subside1Sum + subside0Sum - subside2Sum;

                                            valueDown2 = down;
                                            valueUp2 = up;
                                        }
                                        //Singularity alignment for pentagonal case
                                        else if (currentChartNSides == 5) {
                                            const int& subside0Sum = adjChartSubsideSumResults[currentChartSideId];
                                            const int& subside1Sum = adjChartSubsideSumResults[(currentChartSideId+1)%currentChartNSides];
                                            const int& subside2Sum = adjChartSubsideSumResults[(currentChartSideId+2)%currentChartNSides];
                                            const int& subside3Sum = adjChartSubsideSumResults[(currentChartSideId+3)%currentChartNSides];
                                            const int& subside4Sum = adjChartSubsideSumResults[(currentChartSideId+4)%currentChartNSides];

                                            int down = (subside0Sum + subside1Sum + subside2Sum) - (subside3Sum + subside4Sum);
                                            int up = (subside3Sum + subside4Sum + subside0Sum) - (subside1Sum + subside2Sum);

                                            valueDown2 = down;
                                            valueUp2 = up;
                                        }
                                        //Singularity alignment for hexagonal case
                                        else if (currentChartNSides == 6) {
                                            const int& subside0Sum = adjChartSubsideSumResults[currentChartSideId];
                                            const int& subside2Sum = adjChartSubsideSumResults[(currentChartSideId+2)%currentChartNSides];
                                            const int& subside4Sum = adjChartSubsideSumResults[(currentChartSideId+4)%currentChartNSides];

                                            int down = subside0Sum + subside2Sum - subside4Sum;
                                            int up = subside4Sum + subside0Sum - subside2Sum;

                                            valueDown2 = down;
                                            valueUp2 = up;
                                        }

                                        value1 = valueUp1 - valueDown2;
                                        value2 = valueDown1 - valueUp2;

                                        respected = value1 == 0 && value2 == 0;

                                        valueComputed = true;
                                    }
                                }

                                if (valueComputed) {
                                    if (!respected) {
                                        if (((nSides == 3 || nSides == 5 || nSides == 6) && repeatLosingConstraintsAlign)) {
                                            allRespectConstraints = false;
                                            constraintRespected[constraintRespectedId] = false;
                                            numLostConstraintsAlign++;
                                        }
                                    }

                                    constraintRespectedId++;
                                }
                            }
                        }
                    }
                }

                std::cout << "Lost contraints quad: " << numLostConstraintsQuad << std::endl;
                std::cout << "Lost contraints non quad: " << numLostConstraintsNonQuad << std::endl;
                std::cout << "Lost contraints alignment: " << numLostConstraintsAlign << std::endl << std::endl << std::endl;
            }

#ifndef GUROBI_NON_VERBOSE
            cout << "Support obj: " << supportObj.getValue() << endl;
            cout << "Obj: " << obj.getValue() << endl;
#endif
            gap = model.get(GRB_DoubleAttr_MIPGap);

            status = ILPStatus::SOLUTIONFOUND;

            model.reset();


            for (size_t cId = 0; cId < chartData.charts.size(); cId++) {
                const Chart& chart = chartData.charts[cId];

                bool computable = true;

                for (size_t i = 0; i < chart.chartSubsides.size(); i++) {
                    const size_t subsideId = chart.chartSubsides[i];
                    if (!isComputable[subsideId])
                        computable = false;
                }

                if (computable && chart.faces.size() > 0) {
                    int sizeSum = 0;
                    for (const size_t& subsideId : chart.chartSubsides) {
                        sizeSum += result[subsideId];
                    }

                    if (sizeSum % 2 == 1) {
                        std::cout << "Error not even, chart: " << cId << " -> ";
                        for (const size_t& subsideId : chart.chartSubsides) {
                            std::cout << result[subsideId] << " ";
                        }
                        std::cout << " = " << sizeSum << std::endl;

                        status = ILPStatus::SOLUTIONWRONG;
                    }
                }
            }

        } catch(GRBException e) {
            cout << "Error code = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;

            status = ILPStatus::INFEASIBLE;
        }

        it++;
    }
    while (status == ILPStatus::SOLUTIONFOUND && it < repeatLosingConstraintsIterations + 1 && !allRespectConstraints);

    ilpResults = result;
}

inline void getChartSubsideSum(
        const ChartData& chartData,
        const size_t& cId,
        const std::vector<GRBVar>& vars,
        const std::vector<bool>& isFixed,
        const std::vector<int>& ilpResults,
        const bool hardParityConstraint,
        std::vector<GRBLinExpr>& chartSubsideSum)
{
    const Chart& chart = chartData.charts[cId];

    size_t nSides = chart.chartSides.size();

    chartSubsideSum.resize(nSides, 0.0);

    for (size_t i = 0; i < nSides; i++) {
        const ChartSide& side = chart.chartSides[i];
        GRBLinExpr subsideSum = 0;
        for (const size_t& subsideId : side.subsides) {
            const ChartSubside& subside = chartData.subsides[subsideId];

            if (!isFixed[subsideId]) {
                assert(ilpResults[subsideId] != ILP_IGNORE && !isFixed[subsideId]);
                subsideSum += vars[subsideId];
            }
            else {                
                const int fixedSize = hardParityConstraint ? ilpResults[subsideId] : std::round(ilpResults[subsideId] / 2.0);
                subsideSum += fixedSize;
            }
        }

        chartSubsideSum[i] = subsideSum;
    }
}

inline void getChartSubsideSumResults(
        const ChartData& chartData,
        const size_t& cId,
        const std::vector<int>& results,
        std::vector<int>& chartSubsideSum)
{
    const Chart& chart = chartData.charts[cId];

    size_t nSides = chart.chartSides.size();

    chartSubsideSum.resize(nSides, 0.0);

    for (size_t i = 0; i < nSides; i++) {        
        const ChartSide& side = chart.chartSides[i];
        int subsideSum = 0;
        for (const size_t& subsideId : side.subsides) {
            subsideSum += results[subsideId];
        }

        chartSubsideSum[i] = subsideSum;
    }
}

inline GRBQuadExpr getGurobiCostTermInteger(GRBModel& model, const ILPMethod& method, const GRBLinExpr& value)
{
    GRBQuadExpr expr;

    if (method == LEASTSQUARES) {
        expr = value * value;
    }
    else {
        expr = getGurobiAbsInteger(model, value);
    }

    return expr;
}

inline GRBQuadExpr getGurobiCostTermContinuous(GRBModel& model, const ILPMethod& method, const GRBLinExpr& value)
{
    GRBQuadExpr expr;

    if (method == LEASTSQUARES) {
        expr = value * value;
    }
    else {
        expr = getGurobiAbsContinuous(model, value);
    }

    return expr;
}

inline GRBLinExpr getGurobiAbsInteger(GRBModel& model, const GRBLinExpr& value) {
    GRBVar diff = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER);
    GRBVar abs = model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER);
    model.addConstr(diff == value);
    model.addGenConstrAbs(abs, diff);

    return abs;
}

inline GRBLinExpr getGurobiAbsContinuous(GRBModel& model, const GRBLinExpr& value) {
    GRBVar diff = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
    GRBVar abs = model.addVar(0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
    model.addConstr(diff == value);
    model.addGenConstrAbs(abs, diff);

    return abs;
}

inline GurobiCallBack::GurobiCallBack(const std::vector<float>& times, const std::vector<float>& gaps) :
    times(times), gaps(gaps), index(0)
{

}

inline void GurobiCallBack::callback() {
    try {
        GRBCallback::callback();

        if (where == GRB_CB_MIP) {
            const double runtime = getDoubleInfo(GRB_CB_RUNTIME);

            while (index < times.size() && times[index] <= runtime) {
                index++;
            }

            if (index > 0) {
                const double objbst = getDoubleInfo(GRB_CB_MIP_OBJBST);
                const double objbnd = getDoubleInfo(GRB_CB_MIP_OBJBND);

                if (objbst == GRB_INFINITY)
                    return;

                const double currentGap = std::fabs(objbnd - objbst) / std::fabs(objbst);

                if (currentGap < gaps[index - 1]) {
                    std::cout << "Stop early - Gap: " << currentGap << " < " << gaps[index - 1] << " gap achieved after " << times[index - 1] << " seconds" << std::endl;
                    abort();
                }
            }
        }
    } catch (GRBException e) {
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch (...) {
        std::cout << "Error during callback" << std::endl;
    }
}

}
}


