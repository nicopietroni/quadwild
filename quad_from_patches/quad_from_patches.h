#ifndef QUADFROMPATCHES_H
#define QUADFROMPATCHES_H

#include "includes/common.h"
#include "includes/charts.h"

namespace qfp {

template<class PolyMesh, class TriangleMesh>
void quadrangulationFromPatches(
    TriangleMesh& trimesh,
    const std::vector<std::vector<size_t>>& trimeshPartitions,
    const std::vector<std::vector<size_t>>& trimeshCorners,
    const std::vector<double>& avgLength,
    const Parameters& parameters,
    PolyMesh& quadmesh,
    std::vector<std::vector<size_t>>& quadmeshPartitions,
    std::vector<std::vector<size_t>>& quadmeshCorners,
    std::vector<int>& ilpResult);

template<class TriangleMesh>
ChartData getChartData(
        TriangleMesh& trimesh,
        const std::vector<std::vector<size_t>>& trimeshPartitions,
        const std::vector<std::vector<size_t>>& trimeshCorners);

std::vector<int> findSubdivisions(
        const ChartData& chartData,
        const std::vector<double>& edgeFactor,
        const double alpha,
        const double timeLimit,
        const double gapLimit,
        const ILPMethod& method);

template<class TriangleMesh, class PolyMesh>
void quadrangulate(
        TriangleMesh& trimesh,
        const ChartData& chartData,
        const std::vector<int>& ilpResult,
        const int chartSmoothingIterations,
        const int quadrangulationSmoothingIterations,
        PolyMesh& quadmesh,
        std::vector<std::vector<size_t>>& quadmeshPartitions,
        std::vector<std::vector<size_t>>& quadmeshCorners);


}

#include "quad_from_patches.cpp"

#endif // QUADFROMPATCHES_H
