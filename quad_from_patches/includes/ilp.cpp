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
        const double nonQuadrilateralSimilarityFactor,
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

        //Calculate ideal size
        std::vector<double> idealSize(chartData.subSides.size(), 0.0);
        std::vector<size_t> numberTerms(chartData.subSides.size(), 0);
        for (size_t cId = 0; cId < chartData.charts.size(); cId++) {
            const Chart& chart = chartData.charts[cId];
            if (chart.faces.size() > 0) {
                for (size_t i = 0; i < chart.chartSubSides.size(); i++) {
                    const size_t subsideId = chart.chartSubSides[i];
                    const ChartSubSide& subside = chartData.subSides[subsideId];

                    //If it is not fixed (free)
                    if (!subside.isFixed) {
                        double edgeLength = edgeFactor[cId];

                        double sideSubdivision = subside.length / edgeLength;

                        idealSize[subsideId] += sideSubdivision;
                        numberTerms[subsideId]++;
                    }
                }
            }
        }
        for (size_t i = 0; i < chartData.subSides.size(); i++) {
            idealSize[i] = idealSize[i] / numberTerms[i];
        }

        for (size_t i = 0; i < chartData.subSides.size(); i++) {
            const ChartSubSide& subside = chartData.subSides[i];

            size_t minValue = MINSIDEVALUE;
            size_t maxValue = GRB_INFINITY;

//            size_t minValue = std::max(static_cast<size_t>(MINSIDEVALUE), static_cast<size_t>(std::round(idealSize[i] / 2)));
//            size_t maxValue = std::min(std::max(static_cast<size_t>(4), static_cast<size_t>(std::round(idealSize[i] * 2))), minValue + 4);

            //If it is not a border (free)
            if (!subside.isFixed) {
                vars[i] = model.addVar(minValue, maxValue, 0.0, GRB_INTEGER, "s" + to_string(i));
            }
        }

        std::cout << chartData.subSides.size() << " subsides!" << std::endl;


        const double isoCost = alpha;
        const double regCost = (1 - alpha);

        for (size_t cId = 0; cId < chartData.charts.size(); cId++) {
            const Chart& chart = chartData.charts[cId];
            if (chart.faces.size() > 0) {
                size_t nSides = chart.chartSides.size();

                size_t numRegularityConstraints = 0;
                size_t numIsometryConstraints = 0;

                GRBQuadExpr regExpr = 0;
                GRBQuadExpr isoExpr = 0;

                //Isometry
                if (isometry) {
                    for (size_t i = 0; i < chart.chartSubSides.size(); i++) {
                        const size_t subsideId = chart.chartSubSides[i];
                        const ChartSubSide& subside = chartData.subSides[subsideId];

                        //If it is not fixed (free)
                        if (!subside.isFixed) {
                            numIsometryConstraints++;

                            double edgeLength = edgeFactor[cId];

                            int sideSubdivision = static_cast<int>(std::round(subside.length / edgeLength));

                            size_t dId = diff.size();
                            size_t aId = abs.size();

                            if (method == LEASTSQUARES) {
                                isoExpr += ((vars[subsideId] - sideSubdivision) * (vars[subsideId] - sideSubdivision));
                            }
                            else {
                                diff.push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER, "d" + to_string(dId)));
                                abs.push_back(model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER, "a" + to_string(aId)));

                                model.addConstr(diff[dId] == vars[subsideId] - sideSubdivision, "dc" + to_string(dId));
                                model.addGenConstrAbs(abs[aId], diff[dId], "ac" + to_string(aId));

                                isoExpr += abs[aId];
                            }
                        }
                    }
                }

                //Regularity for quad case
                if (nSides == 4 && regularityForQuadrilaterals) {
                    for (size_t j = 0; j <= 1; j++) {
                        bool areFixed = true;

                        const ChartSide& side1 = chart.chartSides[j];
                        const ChartSide& side2 = chart.chartSides[(j+2)%4];

                        GRBLinExpr subSide1Sum = 0;
                        for (const size_t& subSideId : side1.subsides) {
                            const ChartSubSide& subSide = chartData.subSides[subSideId];
                            if (subSide.isFixed) {
                                subSide1Sum += subSide.size;
                            }
                            else {
                                subSide1Sum += vars[subSideId];
                                areFixed = false;
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
                                areFixed = false;
                            }
                        }

                        if (!areFixed) {
                            numRegularityConstraints++;

                            if (method == LEASTSQUARES) {
                                regExpr += (subSide1Sum - subSide2Sum) * (subSide1Sum - subSide2Sum);
                            }
                            else {
                                size_t dId = diff.size();
                                size_t aId = abs.size();

                                diff.push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER, "d" + to_string(dId)));
                                abs.push_back(model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER, "a" + to_string(aId)));

                                model.addConstr(diff[dId] == subSide1Sum - subSide2Sum, "dc" + to_string(dId));
                                model.addGenConstrAbs(abs[aId], diff[dId], "ac" + to_string(aId));

                                regExpr += abs[aId];
                            }
                        }

                    }
                }
                //Regularity for non-quad case
                else if (nSides != 4 && regularityForNonQuadrilaterals) {
                    for (size_t j = 0; j < nSides; j++) {
                        bool areFixed = true;

                        const ChartSide& side1 = chart.chartSides[j];
                        double length1 = 0;

                        GRBLinExpr subSide1Sum = 0;
                        for (const size_t& subSideId : side1.subsides) {
                            const ChartSubSide& subSide = chartData.subSides[subSideId];

                            length1 += subSide.length;
                            if (subSide.isFixed) {
                                subSide1Sum += subSide.size;
                            }
                            else {
                                subSide1Sum += vars[subSideId];
                                areFixed = false;
                            }
                        }

                        for (size_t k = 0; k < nSides; k++) {
                            const ChartSide& side2 = chart.chartSides[k];
                            double length2 = 0;

                            GRBLinExpr subSide2Sum = 0;
                            for (const size_t& subSideId : side2.subsides) {
                                const ChartSubSide& subSide = chartData.subSides[subSideId];

                                length2 += subSide.length;
                                if (subSide.isFixed) {
                                    subSide2Sum += subSide.size;
                                }
                                else {
                                    subSide2Sum += vars[subSideId];
                                    areFixed = false;
                                }
                            }


                            if (!areFixed) {
                                double factor = std::max(length1, length2) / std::min(length1, length2);
                                assert(factor >= 1);

                                if (factor <= nonQuadrilateralSimilarityFactor) {
                                    numRegularityConstraints++;
                                    if (method == LEASTSQUARES) {
                                        regExpr += ((subSide1Sum - subSide2Sum) * (subSide1Sum - subSide2Sum));
                                    }
                                    else {
                                        size_t dId = diff.size();
                                        size_t aId = abs.size();

                                        diff.push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER, "d" + to_string(dId)));
                                        abs.push_back(model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER, "a" + to_string(aId)));

                                        model.addConstr(diff[dId] == subSide1Sum - subSide2Sum, "dc" + to_string(dId));
                                        model.addGenConstrAbs(abs[aId], diff[dId], "ac" + to_string(aId));

                                        regExpr += abs[aId];
                                    }
                                }
                            }
                        }

                    }
                }

                if (numRegularityConstraints > 0)
                    obj += regCost * regExpr / numRegularityConstraints;
                if (numIsometryConstraints > 0)
                    obj += isoCost * isoExpr / numIsometryConstraints;


                //Even side size sum constraint in a chart
                if (chart.chartSides.size() < 3 || chart.chartSides.size() > 6) {
                    std::cout << "Chart " << cId << " has " << chart.chartSides.size() << " sides." << std::endl;
                    continue;
                }

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

                if (hardParityConstraint) {
                    free[cId] = model.addVar(2, GRB_INFINITY, 0.0, GRB_INTEGER, "f" + to_string(cId));
                    model.addConstr(free[cId]*2 == sumExp);
                }
                else {
                    size_t dId = diff.size();
                    size_t aId = abs.size();

                    free[cId] = model.addVar(2, GRB_INFINITY, 0.0, GRB_INTEGER, "f" + to_string(cId));

                    diff.push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER, "d" + to_string(dId)));
                    abs.push_back(model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER, "a" + to_string(aId)));

                    model.addConstr(diff[dId] == free[cId]*2 - sumExp, "dc" + to_string(dId));
                    model.addGenConstrAbs(abs[aId], diff[dId], "ac" + to_string(aId));

                    obj += abs[aId] * MAXCOST;
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


