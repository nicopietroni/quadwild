#ifndef QUADRETOPOLOGY_H
#define QUADRETOPOLOGY_H

#include <vector>

#include <Eigen/Core>

#include <unordered_map>

#include "includes/qr_charts.h"
#include "includes/qr_ilp.h"
#include "includes/qr_parameters.h"

namespace QuadRetopology {

template<class TriangleMeshType, class PolyMeshType>
std::vector<int> patchDecomposition(
        TriangleMeshType& newSurface,
        PolyMeshType& preservedSurface,
        std::vector<std::vector<size_t>>& partitions,
        std::vector<std::vector<size_t>>& corners,
        const Parameters& parameters);
template<class TriangleMeshType, class PolyMeshType>
std::vector<int> patchDecomposition(
        TriangleMeshType& newSurface,
        PolyMeshType& preservedSurface,
        std::vector<std::vector<size_t>>& partitions,
        std::vector<std::vector<size_t>>& corners,
        const bool initialRemeshing,
        const double edgeFactor,
        const bool reproject,
        const bool splitConcaves,
        const bool finalSmoothing);

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
        const std::vector<size_t> fixedSubsides,
        const size_t& iterations,
        const double& weight);

std::vector<int> findSubdivisions(
        const ChartData& chartData,
        const std::vector<size_t> fixedSubsides,
        const std::vector<double>& chartEdgeLength,
        const Parameters& parameters,
        double& gap);
std::vector<int> findSubdivisions(
        const ChartData& chartData,
        const std::vector<size_t> fixedSubsides,
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
        const double gapLimit,
        const double minimumGap,
        double& gap);

template<class TriangleMeshType, class PolyMeshType>
void quadrangulate(
        TriangleMeshType& newSurface,
        const ChartData& chartData,
        const std::vector<size_t> fixedSubsides,
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
        const std::vector<size_t> fixedSubsides,
        const std::vector<int>& ilpResult,
        const int chartSmoothingIterations,
        const int quadrangulationFixedSmoothingIterations,
        const int quadrangulationNonFixedSmoothingIterations,
        const bool doubletsRemoval,
        PolyMeshType& quadrangulation,
        std::vector<int>& quadrangulationFaceLabel,
        std::vector<std::vector<size_t>>& quadrangulationPartitions,
        std::vector<std::vector<size_t>>& quadrangulationCorners);

template<class PolyMeshType, class TriangleMeshType>
void computeResult(
        PolyMeshType& preservedSurface,
        PolyMeshType& quadrangulatedNewSurface,
        PolyMeshType& result,
        TriangleMeshType& targetBoolean,
        const Parameters& parameters,
        std::vector<int>& resultPreservedVertexMap,
        std::vector<int>& resultPreservedFaceMap);
template<class PolyMeshType, class TriangleMeshType>
void computeResult(
        PolyMeshType& preservedSurface,
        PolyMeshType& quadrangulatedNewSurface,
        PolyMeshType& result,
        TriangleMeshType& targetBoolean,
        const bool doubletRemoval,
        const int resultSmoothingIterations,
        const double resultSmoothingNRing,
        const int resultSmoothingLaplacianIterations,
        const double resultSmoothingLaplacianNRing,
        std::vector<int>& resultPreservedVertexMap,
        std::vector<int>& resultPreservedFaceMap);
}

#include "quadretopology.cpp"

#endif // QUADRETOPOLOGY_H
