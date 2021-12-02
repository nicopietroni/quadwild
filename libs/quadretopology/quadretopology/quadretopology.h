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

#ifndef QUADRETOPOLOGY_H
#define QUADRETOPOLOGY_H

#include <vector>

#include <Eigen/Core>

#include <unordered_map>

#include "includes/qr_charts.h"
#include "includes/qr_ilp.h"
#include "includes/qr_parameters.h"

namespace QuadRetopology {

template<class TriangleMesh>
ChartData computeChartData(
        TriangleMesh& mesh,
        const std::vector<std::vector<size_t>>& meshPartitions,
        const std::vector<std::vector<size_t>>& meshCorners);
template<class TriangleMeshType>
ChartData computeChartData(
        TriangleMeshType& mesh,
        const std::vector<int>& faceLabel,
        const std::vector<std::vector<size_t>>& corners);

std::vector<double> computeChartEdgeLength(
        const ChartData& chartData,
        const size_t& iterations,
        const std::vector<int>& ilpResults,
        const double& weight);

void findSubdivisions(
        const ChartData& chartData,
        const std::vector<double>& chartEdgeLength,
        const Parameters& parameters,
        double& gap,
        std::vector<int>& ilpResults);

void findSubdivisions(
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
        const double minimumGap,
        double& gap,
        std::vector<int>& ilpResults);

template<class TriangleMeshType, class PolyMeshType>
void quadrangulate(
        TriangleMeshType& newSurface,
        const ChartData& chartData,
        const std::vector<size_t> fixedPositionSubsides,
        const std::vector<int>& ilpResult,
        const Parameters& parameters,
        PolyMeshType& quadrangulation,
        std::vector<int>& quadrangulationFaceLabel,
        std::vector<std::vector<size_t>>& quadrangulationPartitions,
        std::vector<std::vector<size_t>>& quadrangulationCorners);
template<class TriangleMeshType, class PolyMeshType>
void quadrangulate(
        TriangleMeshType& newSurface,
        const ChartData& chartData,
        const std::vector<size_t> fixedPositionSubsides,
        const std::vector<int>& ilpResult,
        const int chartSmoothingIterations,
        const int quadrangulationFixedSmoothingIterations,
        const int quadrangulationNonFixedSmoothingIterations,
        const bool doubletsRemoval,
        PolyMeshType& quadrangulation,
        std::vector<int>& quadrangulationFaceLabel,
        std::vector<std::vector<size_t>>& quadrangulationPartitions,
        std::vector<std::vector<size_t>>& quadrangulationCorners);

}

#include "quadretopology.cpp"

#endif // QUADRETOPOLOGY_H
