#include <gurobi_c++.h>

#include "qr_ilp.h"

#define MIN_SUBDIVISION_VALUE 1
#define FEASIBILITY_FIX_COST 1000000.0

namespace QuadRetopology {
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
        ILPStatus& status)
{
    using namespace std;

    vector<int> result(chartData.subsides.size(), -1);
    vector<bool> isFixed(chartData.subsides.size(), false);
    for (size_t subsideId : fixedSubsides) {
        isFixed[subsideId] = true;
    }

    try {
        GRBEnv env = GRBEnv();

        GRBModel model = GRBModel(env);

#ifdef GUROBI_NON_VERBOSE
        model.set(GRB_IntParam_OutputFlag, 0);
#endif

        model.set(GRB_DoubleParam_TimeLimit, timeLimit);
        model.set(GRB_DoubleParam_MIPGap, minimumGap);

        // Create variables
        GRBQuadExpr obj = 0;

        vector<GRBVar> vars(chartData.subsides.size());
        vector<GRBVar> diff;
        vector<GRBVar> abs;

        vector<GRBVar> free(chartData.charts.size());

        std::vector<bool> hasVariable(chartData.subsides.size(), false);
        for (size_t subsideId = 0; subsideId < chartData.subsides.size(); subsideId++) {
            const ChartSubside& subside = chartData.subsides[subsideId];

            //If it is not a border (free)
            if (!isFixed[subsideId]) {
                vars[subsideId] = model.addVar(MIN_SUBDIVISION_VALUE, GRB_INFINITY, 0.0, GRB_INTEGER, "s" + to_string(subsideId));

                hasVariable[subsideId] = true;
            }
            else if (feasibilityFix && subside.isOnBorder) {
                vars[subsideId] = model.addVar(std::max(subside.size - 1, MIN_SUBDIVISION_VALUE), subside.size + 1, 0.0, GRB_INTEGER, "s" + to_string(subsideId));

                size_t dId = diff.size();
                size_t aId = abs.size();
                diff.push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER, "do" + to_string(dId)));
                abs.push_back(model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER, "ao" + to_string(aId)));
                model.addConstr(diff[dId] == vars[subsideId] - subside.size, "doc" + to_string(dId));
                model.addGenConstrAbs(abs[aId], diff[dId], "aoc" + to_string(aId));

                obj += abs[aId] * FEASIBILITY_FIX_COST;

                hasVariable[subsideId] = true;
            }
        }

        std::cout << chartData.subsides.size() << " subsides!" << std::endl;

        const double isoWeight = alpha;
        const double regWeight = (1 - alpha);

        for (size_t cId = 0; cId < chartData.charts.size(); cId++) {
            const Chart& chart = chartData.charts[cId];
            if (chart.faces.size() > 0) {
                size_t nSides = chart.chartSides.size();

                size_t numRegularityTerms = 0;
                size_t numIsometryTerms = 0;

                GRBQuadExpr regExpr = 0;
                GRBQuadExpr isoExpr = 0;

                //Isometry
                if (isometry) {
                    for (size_t i = 0; i < chart.chartSubsides.size(); i++) {
                        const size_t subsideId = chart.chartSubsides[i];
                        const ChartSubside& subside = chartData.subsides[subsideId];

                        //If it is not fixed (free)
                        if (hasVariable[subsideId]) {
                            numIsometryTerms++;

                            double edgeLength = chartEdgeLength[cId];

                            double sideSubdivision = subside.length / edgeLength;
                            if (!hardParityConstraint)
                                sideSubdivision /= 2.0;

                            sideSubdivision = std::max(static_cast<double>(MIN_SUBDIVISION_VALUE), sideSubdivision);

                            size_t dId = diff.size();
                            size_t aId = abs.size();

                            GRBLinExpr value = vars[subsideId] - sideSubdivision;

                            if (method == LEASTSQUARES) {
                                isoExpr += value * value;
                            }
                            else {
                                diff.push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "d" + to_string(dId)));
                                abs.push_back(model.addVar(0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "a" + to_string(aId)));

                                model.addConstr(diff[dId] == value, "dc" + to_string(dId));
                                model.addGenConstrAbs(abs[aId], diff[dId], "ac" + to_string(aId));

                                isoExpr += abs[aId];
                            }
                        }
                    }
                }

                //Regularity for quad case
                if (nSides == 4 && regularityForQuadrilaterals) {
                    for (size_t j = 0; j < 2; j++) {
                        const ChartSide& side1 = chart.chartSides[j];
                        const ChartSide& side2 = chart.chartSides[(j+2)%nSides];

                        GRBLinExpr subside1Sum = 0;
                        for (const size_t& subsideId : side1.subsides) {
                            const ChartSubside& subside = chartData.subsides[subsideId];
                            if (hasVariable[subsideId]) {
                                subside1Sum += vars[subsideId];
                            }
                            else {
                                subside1Sum += subside.size;
                            }
                        }
                        GRBLinExpr subside2Sum = 0;
                        for (const size_t& subsideId : side2.subsides) {
                            const ChartSubside& subside = chartData.subsides[subsideId];
                            if (hasVariable[subsideId]) {
                                subside2Sum += vars[subsideId];
                            }
                            else {
                                subside2Sum += subside.size;
                            }
                        }

                        numRegularityTerms++;

                        GRBLinExpr value = subside1Sum - subside2Sum;

                        if (method == LEASTSQUARES) {
                            regExpr += (value * value) / 2.0;
                        }
                        else {
                            size_t dId = diff.size();
                            size_t aId = abs.size();

                            diff.push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER, "d" + to_string(dId)));
                            abs.push_back(model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER, "a" + to_string(aId)));

                            model.addConstr(diff[dId] == value, "dc" + to_string(dId));
                            model.addGenConstrAbs(abs[aId], diff[dId], "ac" + to_string(aId));

                            regExpr += abs[aId] / 2.0;
                        }
                    }
                }
                else if (nSides == 3 && regularityForNonQuadrilaterals) {
                    //Regularity for triangular case
                    for (size_t j = 0; j < nSides; j++) {
                        const ChartSide& side1 = chart.chartSides[j];
                        const ChartSide& side2 = chart.chartSides[(j+1) % nSides];
                        const ChartSide& side3 = chart.chartSides[(j+2) % nSides];

                        GRBLinExpr subside1Sum = 0;
                        for (const size_t& subsideId : side1.subsides) {
                            const ChartSubside& subside = chartData.subsides[subsideId];
                            if (hasVariable[subsideId]) {
                                subside1Sum += vars[subsideId];
                            }
                            else {
                                subside1Sum += subside.size;
                            }
                        }
                        GRBLinExpr subside2Sum = 0;
                        for (const size_t& subsideId : side2.subsides) {
                            const ChartSubside& subside = chartData.subsides[subsideId];
                            if (hasVariable[subsideId]) {
                                subside2Sum += vars[subsideId];
                            }
                            else {
                                subside2Sum += subside.size;
                            }
                        }
                        GRBLinExpr subside3Sum = 0;
                        for (const size_t& subsideId : side3.subsides) {
                            const ChartSubside& subside = chartData.subsides[subsideId];
                            if (hasVariable[subsideId]) {
                                subside3Sum += vars[subsideId];
                            }
                            else {
                                subside3Sum += subside.size;
                            }
                        }

                        numRegularityTerms++;

                        GRBVar c = model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER, "c" + to_string(numRegularityTerms));
                        model.addConstr(subside1Sum + 1 <= subside2Sum + subside3Sum + c, "cc" + to_string(numRegularityTerms));

                        GRBLinExpr value = c;

                        if (method == LEASTSQUARES) {
                            regExpr += (value * value) * regularityNonQuadrilateralWeight / nSides;
                        }
                        else {
                            size_t dId = diff.size();
                            size_t aId = abs.size();

                            diff.push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER, "d" + to_string(dId)));
                            abs.push_back(model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER, "a" + to_string(aId)));

                            model.addConstr(diff[dId] == value, "dc" + to_string(dId));
                            model.addGenConstrAbs(abs[aId], diff[dId], "ac" + to_string(aId));

                            regExpr += abs[aId] * regularityNonQuadrilateralWeight / nSides;
                        }
                    }
                }
                //Regularity for pentagonal case
                else if (nSides == 5 && regularityForNonQuadrilaterals) {
                    for (size_t j = 0; j < nSides; j++) {
                        const ChartSide& side1 = chart.chartSides[j];
                        const ChartSide& side2 = chart.chartSides[(j+1) % nSides];
                        const ChartSide& side3 = chart.chartSides[(j+2) % nSides];
                        const ChartSide& side4 = chart.chartSides[(j+3) % nSides];
                        const ChartSide& side5 = chart.chartSides[(j+4) % nSides];

                        GRBLinExpr subside1Sum = 0;
                        for (const size_t& subsideId : side1.subsides) {
                            const ChartSubside& subside = chartData.subsides[subsideId];
                            if (hasVariable[subsideId]) {
                                subside1Sum += vars[subsideId];
                            }
                            else {
                                subside1Sum += subside.size;
                            }
                        }
                        GRBLinExpr subside2Sum = 0;
                        for (const size_t& subsideId : side2.subsides) {
                            const ChartSubside& subside = chartData.subsides[subsideId];
                            if (hasVariable[subsideId]) {
                                subside2Sum += vars[subsideId];
                            }
                            else {
                                subside2Sum += subside.size;
                            }
                        }
                        GRBLinExpr subside3Sum = 0;
                        for (const size_t& subsideId : side3.subsides) {
                            const ChartSubside& subside = chartData.subsides[subsideId];
                            if (hasVariable[subsideId]) {
                                subside3Sum += vars[subsideId];
                            }
                            else {
                                subside3Sum += subside.size;
                            }
                        }

                        GRBLinExpr subside4Sum = 0;
                        for (const size_t& subsideId : side4.subsides) {
                            const ChartSubside& subside = chartData.subsides[subsideId];
                            if (hasVariable[subsideId]) {
                                subside4Sum += vars[subsideId];
                            }
                            else {
                                subside4Sum += subside.size;
                            }
                        }

                        GRBLinExpr subside5Sum = 0;
                        for (const size_t& subsideId : side5.subsides) {
                            const ChartSubside& subside = chartData.subsides[subsideId];
                            if (hasVariable[subsideId]) {
                                subside5Sum += vars[subsideId];
                            }
                            else {
                                subside5Sum += subside.size;
                            }
                        }

                        numRegularityTerms++;

                        GRBVar c = model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER, "c" + to_string(numRegularityTerms));
                        model.addConstr(subside1Sum + subside2Sum + 1 <= subside3Sum + subside4Sum + subside5Sum + c, "cc" + to_string(numRegularityTerms));

                        GRBLinExpr value = c;

                        if (method == LEASTSQUARES) {
                            regExpr += (value * value) * regularityNonQuadrilateralWeight / nSides;
                        }
                        else {
                            size_t dId = diff.size();
                            size_t aId = abs.size();

                            diff.push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER, "d" + to_string(dId)));
                            abs.push_back(model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER, "a" + to_string(aId)));

                            model.addConstr(diff[dId] == value, "dc" + to_string(dId));
                            model.addGenConstrAbs(abs[aId], diff[dId], "ac" + to_string(aId));

                            regExpr += abs[aId] * regularityNonQuadrilateralWeight / nSides;
                        }
                    }
                }
                //Regularity for hexagonal case
                else if (nSides == 6 && regularityForNonQuadrilaterals) {
                    for (size_t j = 0; j < nSides; j++) {
                        const ChartSide& side1 = chart.chartSides[j];
                        const ChartSide& side2 = chart.chartSides[(j+1) % nSides];
                        const ChartSide& side3 = chart.chartSides[(j+2) % nSides];
                        const ChartSide& side4 = chart.chartSides[(j+3) % nSides];
                        const ChartSide& side5 = chart.chartSides[(j+4) % nSides];
                        const ChartSide& side6 = chart.chartSides[(j+5) % nSides];

                        GRBLinExpr subside1Sum = 0;
                        for (const size_t& subsideId : side1.subsides) {
                            const ChartSubside& subside = chartData.subsides[subsideId];
                            if (hasVariable[subsideId]) {
                                subside1Sum += vars[subsideId];
                            }
                            else {
                                subside1Sum += subside.size;
                            }
                        }
                        GRBLinExpr subside2Sum = 0;
                        for (const size_t& subsideId : side2.subsides) {
                            const ChartSubside& subside = chartData.subsides[subsideId];
                            if (hasVariable[subsideId]) {
                                subside2Sum += vars[subsideId];
                            }
                            else {
                                subside2Sum += subside.size;
                            }
                        }
                        GRBLinExpr subside3Sum = 0;
                        for (const size_t& subsideId : side3.subsides) {
                            const ChartSubside& subside = chartData.subsides[subsideId];
                            if (hasVariable[subsideId]) {
                                subside3Sum += vars[subsideId];
                            }
                            else {
                                subside3Sum += subside.size;
                            }
                        }

                        GRBLinExpr subside4Sum = 0;
                        for (const size_t& subsideId : side4.subsides) {
                            const ChartSubside& subside = chartData.subsides[subsideId];
                            if (hasVariable[subsideId]) {
                                subside4Sum += vars[subsideId];
                            }
                            else {
                                subside4Sum += subside.size;
                            }
                        }

                        GRBLinExpr subside5Sum = 0;
                        for (const size_t& subsideId : side5.subsides) {
                            const ChartSubside& subside = chartData.subsides[subsideId];
                            if (hasVariable[subsideId]) {
                                subside5Sum += vars[subsideId];
                            }
                            else {
                                subside5Sum += subside.size;
                            }
                        }

                        GRBLinExpr subside6Sum = 0;
                        for (const size_t& subsideId : side6.subsides) {
                            const ChartSubside& subside = chartData.subsides[subsideId];
                            if (hasVariable[subsideId]) {
                                subside6Sum += vars[subsideId];
                            }
                            else {
                                subside6Sum += subside.size;
                            }
                        }

                        numRegularityTerms++;

                        GRBVar c = model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER, "c" + to_string(numRegularityTerms));
                        model.addConstr(subside1Sum + 1 <= subside3Sum + subside5Sum + c, "cc" + to_string(numRegularityTerms));

                        GRBLinExpr parityEquation = subside1Sum + subside3Sum + subside5Sum;
                        GRBVar hexParity = model.addVar(3, GRB_INFINITY, 0.0, GRB_INTEGER, "pc" + to_string(numRegularityTerms));
                        GRBVar hexFree = model.addVar(0, 1, 0.0, GRB_INTEGER, "pf" + to_string(numRegularityTerms));
                        model.addConstr(hexParity * 2 == parityEquation + hexFree);

                        GRBLinExpr value = (c + hexFree) / 2.0;

                        if (method == LEASTSQUARES) {
                            regExpr += (value * value) * regularityNonQuadrilateralWeight / nSides;
                        }
                        else {
                            size_t dId = diff.size();
                            size_t aId = abs.size();

                            diff.push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER, "d" + to_string(dId)));
                            abs.push_back(model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER, "a" + to_string(aId)));

                            model.addConstr(diff[dId] == value, "dc" + to_string(dId));
                            model.addGenConstrAbs(abs[aId], diff[dId], "ac" + to_string(aId));

                            regExpr += abs[aId] * regularityNonQuadrilateralWeight / nSides;
                        }
                    }
                }

                if (numRegularityTerms > 0)
                    obj += regWeight * regExpr / numRegularityTerms;
                if (numIsometryTerms > 0)
                    obj += isoWeight * isoExpr / numIsometryTerms;

                //Even side size sum constraint in a chart
                if (hardParityConstraint) {
                    GRBLinExpr sumExp = 0;
                    for (const size_t& subsideId : chart.chartSubsides) {
                        const ChartSubside& subside = chartData.subsides[subsideId];
                        if (hasVariable[subsideId]) {
                            sumExp += vars[subsideId];
                        }
                        else {
                            sumExp += subside.size;
                        }
                    }

                    free[cId] = model.addVar(2, GRB_INFINITY, 0.0, GRB_INTEGER, "f" + to_string(cId));
                    model.addConstr(free[cId] * 2 == sumExp);
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
        model.setObjective(obj, GRB_MINIMIZE);

//        model.write("out.lp");

        //Optimize model
        model.optimize();

        for (size_t subsideId = 0; subsideId < chartData.subsides.size(); subsideId++) {
            const ChartSubside& subside = chartData.subsides[subsideId];
            if (hasVariable[subsideId]) {
                result[subsideId] = static_cast<int>(std::round(vars[subsideId].get(GRB_DoubleAttr_X)));

                if (!hardParityConstraint) {
                    result[subsideId] *= 2;
                }
            }
            else {
                result[subsideId] = subside.size;
            }
        }

        cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

        gap = model.get(GRB_DoubleAttr_MIPGap);

        status = ILPStatus::SOLUTIONFOUND;


        for (size_t cId = 0; cId < chartData.charts.size(); cId++) {
            const Chart& chart = chartData.charts[cId];

            if (chart.faces.size() > 0) {
                int sizeSum = 0;
                for (const size_t& subsideId : chart.chartSubsides) {
                    sizeSum += result[subsideId];
                }

                if (sizeSum % 2 == 1) {
                    std::cout << "Error not even, chart: " << cId << " -> ";
                    for (const size_t& subsideId : chart.chartSubsides) {
                        std::cout << result[subsideId] << " ";
                    }
                    std::cout << " = " << sizeSum << " - FREE: " << free[cId].get(GRB_DoubleAttr_X) << std::endl;

                    status = ILPStatus::SOLUTIONWRONG;
                }
            }
        }

    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;

        status = ILPStatus::INFEASIBLE;
    }

    return result;
}

}
}


