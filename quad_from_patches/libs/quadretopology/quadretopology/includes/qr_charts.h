#ifndef QR_CHARTS_H
#define QR_CHARTS_H

#include <vector>
#include <array>
#include <set>
#include <cmath>

namespace QuadRetopology {

struct ChartSubside {
    std::array<int, 2> incidentCharts;
    std::array<int, 2> incidentChartSideId;
    std::vector<size_t> vertices;

    double length;
    int size;

    bool isOnBorder;
};

struct ChartSide {
    std::vector<size_t> vertices;
    std::vector<size_t> subsides;
    std::vector<bool> reversedSubside;

    double length;
    int size;
};

struct Chart {
    std::vector<size_t> faces;
    std::vector<size_t> borderFaces;

    std::vector<size_t> adjacentCharts;


    std::vector<ChartSide> chartSides;
    std::vector<size_t> chartSubsides;

    int label;
};

struct ChartData {
    std::set<int> labels;
    std::vector<Chart> charts;
    std::vector<ChartSubside> subsides;
};

template<class TriangleMeshType>
void findChartFacesAndBorderFaces(
        TriangleMeshType& mesh,
        const std::vector<int>& faceLabel,
        ChartData& chartData);

}

#include "qr_charts.cpp"

#endif // QR_CHARTS_H
