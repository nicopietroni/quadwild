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
        const double alpha,
        const ILPMethod& method,
        const bool regularity,
        const double timeLimit,
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

        // Create variables
        GRBQuadExpr obj = 0;

        vector<GRBVar> vars(chartData.subSides.size());
        vector<GRBVar> diff;
        vector<GRBVar> abs;

        vector<GRBVar> free(chartData.charts.size());

        for (size_t i = 0; i < chartData.subSides.size(); i++) {
            const ChartSubSide& subside = chartData.subSides[i];

            //If it is not a border (free)
            if (!subside.isFixed) {
                assert(subside.incidentCharts[0] >= 0 && subside.incidentCharts[1] >= 0);

                vars[i] = model.addVar(MINSIDEVALUE, GRB_INFINITY, 0.0, GRB_INTEGER, "s" + to_string(i));
            }
        }

        std::cout << chartData.subSides.size() << " subsides!" << std::endl;

        if (regularity) {
            GRBQuadExpr regExpr = 0;
            const double regCost = alpha;
            //Regularity
            for (size_t i = 0; i < chartData.subSides.size(); i++) {
                const ChartSubSide& subside = chartData.subSides[i];

                //If it is not a border (free)
                if (!subside.isFixed) {
                    size_t nIncidents = 0;
                    double incidentChartAverageLength;

                    if (subside.incidentCharts[0] >= 0) {
                        assert(edgeFactor[subside.incidentCharts[0]] > 0);
                        incidentChartAverageLength += edgeFactor[subside.incidentCharts[0]];
                        nIncidents++;
                    }
                    if (subside.incidentCharts[1] >= 0) {
                        assert(edgeFactor[subside.incidentCharts[1]] > 0);
                        incidentChartAverageLength += edgeFactor[subside.incidentCharts[1]];
                        nIncidents++;
                    }

                    assert(nIncidents > 0);
                    incidentChartAverageLength /= nIncidents;

                    int sideSubdivision = static_cast<int>(std::round(subside.length / incidentChartAverageLength));

                    size_t dId = diff.size();
                    size_t aId = abs.size();

                    if (method == LEASTSQUARES) {
                        regExpr += (vars[i] - sideSubdivision) * (vars[i] - sideSubdivision);
                    }
                    else {
                        diff.push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER, "d" + to_string(dId)));
                        abs.push_back(model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER, "a" + to_string(aId)));

                        model.addConstr(diff[dId] == vars[i] - sideSubdivision, "dc" + to_string(dId));
                        model.addGenConstrAbs(abs[aId], diff[dId], "ac" + to_string(aId));

                        regExpr += abs[aId];
                    }
                }
            }
            obj += regCost * regExpr;
        }

        GRBQuadExpr isoExpr = 0;
        const double isoCost = (1-alpha);
        for (size_t i = 0; i < chartData.charts.size(); i++) {
            const Chart& chart = chartData.charts[i];
            if (chart.faces.size() > 0) {
                size_t nSides = chart.chartSides.size();

                //Quad case
                if (nSides == 4) {

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
                            if (method == LEASTSQUARES) {
                                isoExpr += (subSide1Sum - subSide2Sum) * (subSide1Sum - subSide2Sum);
                            }
                            else {
                                size_t dId = diff.size();
                                size_t aId = abs.size();

                                diff.push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER, "d" + to_string(dId)));
                                abs.push_back(model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER, "a" + to_string(aId)));

                                model.addConstr(diff[dId] == subSide1Sum - subSide2Sum, "dc" + to_string(dId));
                                model.addGenConstrAbs(abs[aId], diff[dId], "ac" + to_string(aId));

                                isoExpr += abs[aId];
                            }
                        }

                    }
                }
            }

            //Even side size sum constraint in a chart
            for (size_t i = 0; i < chartData.charts.size(); i++) {
                const Chart& chart = chartData.charts[i];
                if (chart.faces.size() > 0) {
                    if (chart.chartSides.size() < 3 || chart.chartSides.size() > 6) {
                        std::cout << "Chart " << i << " has " << chart.chartSides.size() << " sides." << std::endl;
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
                    free[i] = model.addVar(2, GRB_INFINITY, 0.0, GRB_INTEGER, "f" + to_string(i));
                    model.addConstr(free[i]*2 == sumExp);
                }
            }
        }
        obj += isoCost * isoExpr;


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


