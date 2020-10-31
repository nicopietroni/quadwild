#include "quad_from_patches.h"

#include <quadretopology/quad_retopology.h>

#ifdef SAVE_MESHES_FOR_DEBUG
#include <igl/writeOBJ.h>
#include <wrap/io_trimesh/export_obj.h>
#endif


namespace qfp {

template<class PolyMesh, class TriangleMesh>
void quadrangulationFromPatches(
    TriangleMesh& trimesh,
    const std::vector<std::vector<size_t>>& trimeshPartitions,
    const std::vector<std::vector<size_t>>& trimeshCorners,
    const std::vector<double>& chartEdgeLength,
    const QuadRetopology::Parameters& parameters,
    PolyMesh& quadmesh,
    std::vector<std::vector<size_t>>& quadmeshPartitions,
    std::vector<std::vector<size_t>>& quadmeshCorners,
    std::vector<int>& ilpResult)
{
    assert(trimeshPartitions.size() == trimeshCorners.size() && chartEdgeLength.size() == trimeshPartitions.size());

    //Get chart data
    QuadRetopology::ChartData chartData = QuadRetopology::computeChartData(
            trimesh,
            trimeshPartitions,
            trimeshCorners);

    //Solve ILP to find best side size
    std::vector<size_t> fixedSubsides;
    double gap;
    ilpResult = QuadRetopology::findSubdivisions(
            chartData,
            fixedSubsides,
            chartEdgeLength,
            parameters,
            gap);

    //Quadrangulate
    std::vector<int> quadmeshLabel;
    QuadRetopology::quadrangulate(
            trimesh,
            chartData,
            fixedSubsides,
            ilpResult,
            parameters,
            quadmesh,
            quadmeshLabel,
            quadmeshPartitions,
            quadmeshCorners);
}

}
