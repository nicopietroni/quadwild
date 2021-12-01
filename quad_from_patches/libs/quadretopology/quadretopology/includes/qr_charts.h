/***************************************************************************/
/* Copyright(C) 2021


The authors of

Reliable Feature-Line Driven Quad-Remeshing
Siggraph 2021


 All rights reserved.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
****************************************************************************/

#ifndef QR_CHARTS_H
#define QR_CHARTS_H

#include <vector>
#include <array>
#include <set>
#include <cmath>

namespace QuadRetopology {

struct ChartSubside {
    std::array<int, 2> incidentCharts;
    std::array<int, 2> incidentChartSubsideId;
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
