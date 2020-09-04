#include "ilp.h"

#include <gurobi_c++.h>

#include "utils.h"

#define MINSIDEVALUE 1
#define AVERAGELENGTHSMOOTHITERATIONS 1
#define MAXCOST 1000000.0

namespace qfp {

inline std::vector<int> solveILP(
        const ChartData& chartData,
        const std::vector<double>& edgeFactor,
        const ILPMethod& method,
        const double alpha,
        const bool isometry,
        const bool regularityForQuadrilaterals,
        const bool regularityForNonQuadrilaterals,
        const double regularityNonQuadrilateralWeight,
        const bool hardParityConstraint,
        const double timeLimit,
        const double minimumGap,
        double& gap,
        ILPStatus& status)
{
    using namespace std;

    vector<int> result(chartData.subSides.size(), -1);

    try {
        GRBEnv env = GRBEnv();

        GRBModel model = GRBModel(env);

//        model.set(GRB_IntParam_OutputFlag, 0);

        model.set(GRB_DoubleParam_TimeLimit, timeLimit);
        model.set(GRB_DoubleParam_MIPGap, minimumGap);

        // Create variables
        GRBQuadExpr obj = 0;

        vector<GRBVar> vars(chartData.subSides.size());
        vector<GRBVar> diff;
        vector<GRBVar> abs;

        vector<GRBVar> free(chartData.charts.size());

//        //Calculate ideal size
//        std::vector<double> idealMinSize(chartData.subSides.size(), std::numeric_limits<double>::max());
//        std::vector<double> idealMaxSize(chartData.subSides.size(), 0.0);
//        for (size_t cId = 0; cId < chartData.charts.size(); cId++) {
//            const Chart& chart = chartData.charts[cId];
//            if (chart.faces.size() > 0) {
//                for (size_t i = 0; i < chart.chartSubSides.size(); i++) {
//                    const size_t subsideId = chart.chartSubSides[i];
//                    const ChartSubSide& subside = chartData.subSides[subsideId];

//                    //If it is not fixed (free)
//                    if (!subside.isFixed) {
//                        double edgeLength = edgeFactor[cId];

//                        double sideSubdivision = subside.length / edgeLength;

//                        idealMinSize[subsideId] = std::min(sideSubdivision, idealMinSize[subsideId]);
//                        idealMaxSize[subsideId] = std::max(sideSubdivision, idealMaxSize[subsideId]);
//                    }
//                }
//            }
//        }

        for (size_t i = 0; i < chartData.subSides.size(); i++) {
            const ChartSubSide& subside = chartData.subSides[i];

            double minValue = MINSIDEVALUE;
            double maxValue = GRB_INFINITY;

//            size_t minValue = std::max(static_cast<size_t>(MINSIDEVALUE), static_cast<size_t>(std::round(idealMinSize[i] / 2)));
//            size_t maxValue = std::min(std::max(static_cast<size_t>(4), static_cast<size_t>(std::round(idealMaxSize[i] * 2))), minValue + 2);

            //If it is not a border (free)
            if (!subside.isFixed) {
                vars[i] = model.addVar(minValue, maxValue, 0.0, GRB_INTEGER, "s" + to_string(i));
            }
        }

        std::cout << chartData.subSides.size() << " subsides!" << std::endl;


        const double isoWeight = alpha;
        const double regWeight = (1 - alpha);

        for (size_t cId = 0; cId < chartData.charts.size(); cId++) {
            const Chart& chart = chartData.charts[cId];
            if (chart.faces.size() > 0) {
                size_t nSides = chart.chartSides.size();

                size_t numRegularityCosts = 0;
                size_t numIsometryCosts = 0;

                GRBQuadExpr regExpr = 0;
                GRBQuadExpr isoExpr = 0;

                //Isometry
                if (isometry) {
                    for (size_t i = 0; i < chart.chartSubSides.size(); i++) {
                        const size_t subsideId = chart.chartSubSides[i];
                        const ChartSubSide& subside = chartData.subSides[subsideId];

                        //If it is not fixed (free)
                        if (!subside.isFixed) {
                            numIsometryCosts++;

                            double edgeLength = edgeFactor[cId];

                            double sideSubdivision = subside.length / edgeLength;
                            if (!hardParityConstraint)
                                sideSubdivision /= 2.0;

                            sideSubdivision = std::max(static_cast<double>(MINSIDEVALUE), sideSubdivision);

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
                if (nSides == 4) {
                    if (regularityForQuadrilaterals) {
                        for (size_t j = 0; j < 2; j++) {
                            const ChartSide& side1 = chart.chartSides[j];
                            const ChartSide& side2 = chart.chartSides[(j+2)%nSides];

                            GRBLinExpr subSide1Sum = 0;
                            for (const size_t& subSideId : side1.subsides) {
                                const ChartSubSide& subSide = chartData.subSides[subSideId];
                                if (subSide.isFixed) {
                                    subSide1Sum += subSide.size;
                                }
                                else {
                                    subSide1Sum += vars[subSideId];
                                }
                            }
                            GRBLinExpr subSide2Sum = 0;
                            for (const size_t& subSideId : side2.subsides) {
                                const ChartSubSide& subSide = chartData.subSides[subSideId];
                                if (subSide.isFixed) {
                                    subSide2Sum += subSide.size;
                                }
                                else {
                                    subSide2Sum += vars[subSideId];
                                }
                            }

                            numRegularityCosts++;

                            GRBLinExpr value = subSide1Sum - subSide2Sum;

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
                }
                else if (regularityForNonQuadrilaterals) {
                    //Regularity for triangular case
                    if (nSides == 3) {
                        for (size_t j = 0; j < nSides; j++) {
                            const ChartSide& side1 = chart.chartSides[j];
                            const ChartSide& side2 = chart.chartSides[(j+1) % nSides];
                            const ChartSide& side3 = chart.chartSides[(j+2) % nSides];

                            GRBLinExpr subSide1Sum = 0;
                            for (const size_t& subSideId : side1.subsides) {
                                const ChartSubSide& subSide = chartData.subSides[subSideId];
                                if (subSide.isFixed) {
                                    subSide1Sum += subSide.size;
                                }
                                else {
                                    subSide1Sum += vars[subSideId];
                                }
                            }
                            GRBLinExpr subSide2Sum = 0;
                            for (const size_t& subSideId : side2.subsides) {
                                const ChartSubSide& subSide = chartData.subSides[subSideId];
                                if (subSide.isFixed) {
                                    subSide2Sum += subSide.size;
                                }
                                else {
                                    subSide2Sum += vars[subSideId];
                                }
                            }
                            GRBLinExpr subSide3Sum = 0;
                            for (const size_t& subSideId : side3.subsides) {
                                const ChartSubSide& subSide = chartData.subSides[subSideId];
                                if (subSide.isFixed) {
                                    subSide3Sum += subSide.size;
                                }
                                else {
                                    subSide3Sum += vars[subSideId];
                                }
                            }

                            numRegularityCosts++;

                            GRBVar c = model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER, "c" + to_string(numRegularityCosts));
                            model.addConstr(subSide1Sum + 1 <= subSide2Sum + subSide3Sum + c, "cc" + to_string(numRegularityCosts));

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
                    else if (nSides == 5) {
                        for (size_t j = 0; j < nSides; j++) {
                            const ChartSide& side1 = chart.chartSides[j];
                            const ChartSide& side2 = chart.chartSides[(j+1) % nSides];
                            const ChartSide& side3 = chart.chartSides[(j+2) % nSides];
                            const ChartSide& side4 = chart.chartSides[(j+3) % nSides];
                            const ChartSide& side5 = chart.chartSides[(j+4) % nSides];

                            GRBLinExpr subSide1Sum = 0;
                            for (const size_t& subSideId : side1.subsides) {
                                const ChartSubSide& subSide = chartData.subSides[subSideId];
                                if (subSide.isFixed) {
                                    subSide1Sum += subSide.size;
                                }
                                else {
                                    subSide1Sum += vars[subSideId];
                                }
                            }
                            GRBLinExpr subSide2Sum = 0;
                            for (const size_t& subSideId : side2.subsides) {
                                const ChartSubSide& subSide = chartData.subSides[subSideId];
                                if (subSide.isFixed) {
                                    subSide2Sum += subSide.size;
                                }
                                else {
                                    subSide2Sum += vars[subSideId];
                                }
                            }
                            GRBLinExpr subSide3Sum = 0;
                            for (const size_t& subSideId : side3.subsides) {
                                const ChartSubSide& subSide = chartData.subSides[subSideId];
                                if (subSide.isFixed) {
                                    subSide3Sum += subSide.size;
                                }
                                else {
                                    subSide3Sum += vars[subSideId];
                                }
                            }

                            GRBLinExpr subSide4Sum = 0;
                            for (const size_t& subSideId : side4.subsides) {
                                const ChartSubSide& subSide = chartData.subSides[subSideId];
                                if (subSide.isFixed) {
                                    subSide4Sum += subSide.size;
                                }
                                else {
                                    subSide4Sum += vars[subSideId];
                                }
                            }

                            GRBLinExpr subSide5Sum = 0;
                            for (const size_t& subSideId : side5.subsides) {
                                const ChartSubSide& subSide = chartData.subSides[subSideId];
                                if (subSide.isFixed) {
                                    subSide5Sum += subSide.size;
                                }
                                else {
                                    subSide5Sum += vars[subSideId];
                                }
                            }

                            numRegularityCosts++;

                            GRBVar c = model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER, "c" + to_string(numRegularityCosts));
                            model.addConstr(subSide1Sum + subSide2Sum + 1 <= subSide3Sum + subSide4Sum + subSide5Sum + c, "cc" + to_string(numRegularityCosts));

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
                    else if (nSides == 6) {
                        for (size_t j = 0; j < nSides; j++) {
                            const ChartSide& side1 = chart.chartSides[j];
                            const ChartSide& side2 = chart.chartSides[(j+1) % nSides];
                            const ChartSide& side3 = chart.chartSides[(j+2) % nSides];
                            const ChartSide& side4 = chart.chartSides[(j+3) % nSides];
                            const ChartSide& side5 = chart.chartSides[(j+4) % nSides];
                            const ChartSide& side6 = chart.chartSides[(j+5) % nSides];

                            GRBLinExpr subSide1Sum = 0;
                            for (const size_t& subSideId : side1.subsides) {
                                const ChartSubSide& subSide = chartData.subSides[subSideId];
                                if (subSide.isFixed) {
                                    subSide1Sum += subSide.size;
                                }
                                else {
                                    subSide1Sum += vars[subSideId];
                                }
                            }
                            GRBLinExpr subSide2Sum = 0;
                            for (const size_t& subSideId : side2.subsides) {
                                const ChartSubSide& subSide = chartData.subSides[subSideId];
                                if (subSide.isFixed) {
                                    subSide2Sum += subSide.size;
                                }
                                else {
                                    subSide2Sum += vars[subSideId];
                                }
                            }
                            GRBLinExpr subSide3Sum = 0;
                            for (const size_t& subSideId : side3.subsides) {
                                const ChartSubSide& subSide = chartData.subSides[subSideId];
                                if (subSide.isFixed) {
                                    subSide3Sum += subSide.size;
                                }
                                else {
                                    subSide3Sum += vars[subSideId];
                                }
                            }

                            GRBLinExpr subSide4Sum = 0;
                            for (const size_t& subSideId : side4.subsides) {
                                const ChartSubSide& subSide = chartData.subSides[subSideId];
                                if (subSide.isFixed) {
                                    subSide4Sum += subSide.size;
                                }
                                else {
                                    subSide4Sum += vars[subSideId];
                                }
                            }

                            GRBLinExpr subSide5Sum = 0;
                            for (const size_t& subSideId : side5.subsides) {
                                const ChartSubSide& subSide = chartData.subSides[subSideId];
                                if (subSide.isFixed) {
                                    subSide5Sum += subSide.size;
                                }
                                else {
                                    subSide5Sum += vars[subSideId];
                                }
                            }

                            GRBLinExpr subSide6Sum = 0;
                            for (const size_t& subSideId : side6.subsides) {
                                const ChartSubSide& subSide = chartData.subSides[subSideId];
                                if (subSide.isFixed) {
                                    subSide6Sum += subSide.size;
                                }
                                else {
                                    subSide6Sum += vars[subSideId];
                                }
                            }

                            numRegularityCosts++;

                            GRBVar c = model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER, "c" + to_string(numRegularityCosts));
                            model.addConstr(subSide1Sum + 1 <= subSide3Sum + subSide5Sum + c, "cc" + to_string(numRegularityCosts));

                            GRBLinExpr parityEquation = subSide1Sum + subSide3Sum + subSide5Sum;
                            GRBVar hexParity = model.addVar(3, GRB_INFINITY, 0.0, GRB_INTEGER, "pc" + to_string(numRegularityCosts));
                            GRBVar hexFree = model.addVar(0, 1, 0.0, GRB_INTEGER, "pf" + to_string(numRegularityCosts));
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
                }

                if (numRegularityCosts > 0)
                    obj += regWeight * regExpr / numRegularityCosts;
                if (numIsometryCosts > 0)
                    obj += isoWeight * isoExpr / numIsometryCosts;


                //Even side size sum constraint in a chart
                if (chart.chartSides.size() < 3 || chart.chartSides.size() > 6) {
                    std::cout << "Chart " << cId << " has " << chart.chartSides.size() << " sides." << std::endl;
                    continue;
                }

                if (hardParityConstraint) {
                    GRBLinExpr sumExp = 0;
                    for (const size_t& subSideId : chart.chartSubSides) {
                        const ChartSubSide& subSide = chartData.subSides[subSideId];
                        if (subSide.isFixed) {
                            sumExp += subSide.size;
                        }
                        else {
                            sumExp += vars[subSideId];
                        }
                    }

                    free[cId] = model.addVar(2, GRB_INFINITY, 0.0, GRB_INTEGER, "f" + to_string(cId));
                    model.addConstr(free[cId] * 2 == sumExp);
                }

            }
        }

//        model.update();

        //Set objective function
        model.setObjective(obj, GRB_MINIMIZE);

//        model.write("out.lp");

        //Optimize model
        model.optimize();

        for (size_t i = 0; i < chartData.subSides.size(); i++) {
            const ChartSubSide& subSide = chartData.subSides[i];
            if (subSide.isFixed) {
                result[i] = subSide.size;
            }
            else {
                result[i] = static_cast<int>(std::round(vars[i].get(GRB_DoubleAttr_X)));
                if (!hardParityConstraint) {
                    result[i] *= 2;
                }
            }
        }

        cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

        gap = model.get(GRB_DoubleAttr_MIPGap);

        status = ILPStatus::SOLUTIONFOUND;


        for (size_t i = 0; i < chartData.charts.size(); i++) {
            const Chart& chart = chartData.charts[i];

            if (chart.faces.size() > 0) {
                int sizeSum = 0;
                for (const size_t& subSideId : chart.chartSubSides) {
                    sizeSum += result[subSideId];
                }

                if (sizeSum % 2 == 1) {
                    std::cout << "Error not even, chart: " << i << " -> ";
                    for (const size_t& subSideId : chart.chartSubSides) {
                        std::cout << result[subSideId] << " ";
                    }
                    std::cout << " = " << sizeSum << " - FREE: " << free[i].get(GRB_DoubleAttr_X) << std::endl;

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


