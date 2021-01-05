#include "quad_from_patches.h"

#include <quadretopology/quadretopology.h>

#include <random>

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
    const float fixToIsometryRatio,
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

    //Fix subsides
    std::vector<size_t> fixedSubdivisionSubsides;

    int numToFix = std::round(chartData.subsides.size() * fixToIsometryRatio);
    if (numToFix > 0) {
        std::vector<bool> subsideAlreadyFixed(chartData.subsides.size(), false);
        std::vector<bool> chartAlreadyFixed(chartData.charts.size(), false);

        int numFixed = 0;
        bool foundCandidate;
        do {
            foundCandidate = false;

            std::vector<size_t> candidateSubsides;
            for (size_t subsideId = 0; subsideId < chartData.subsides.size(); subsideId++) {
                std::array<int, 2> incidentCharts = chartData.subsides[subsideId].incidentCharts;

                if (!subsideAlreadyFixed[subsideId] &&
                    (incidentCharts[0] == -1 || !chartAlreadyFixed[incidentCharts[0]]) &&
                    (incidentCharts[1] == -1 || !chartAlreadyFixed[incidentCharts[1]]))
                {
                    candidateSubsides.push_back(subsideId);
                }
            }

            if (numToFix > numFixed && !candidateSubsides.empty()) {
                std::mt19937 generator(0);
                std::uniform_int_distribution<int> distribution(0, candidateSubsides.size());

                size_t targetSubsideId = candidateSubsides[distribution(generator)];

                std::array<int, 2> incidentCharts = chartData.subsides[targetSubsideId].incidentCharts;
                const QuadRetopology::ChartSubside& subside = chartData.subsides[targetSubsideId];

                double edgeLength = 0.0;
                int numIncident = 0;
                for (int cId : incidentCharts) {
                    if (cId >= 0) {
                        edgeLength += chartEdgeLength[cId];
                        numIncident++;
                    }
                }
                edgeLength /= numIncident;

                double sideSubdivision = subside.length / edgeLength;
                if (!parameters.hardParityConstraint) {
                    sideSubdivision /= 2.0;
                }

                sideSubdivision = std::max(static_cast<double>(MIN_SUBDIVISION_VALUE), sideSubdivision);

                chartData.subsides[targetSubsideId].size = sideSubdivision;
                fixedSubdivisionSubsides.push_back(targetSubsideId);

                subsideAlreadyFixed[targetSubsideId] = true;
                if (incidentCharts[0] != -1)
                    chartAlreadyFixed[incidentCharts[0]] = true;
                if (incidentCharts[1] != -1)
                    chartAlreadyFixed[incidentCharts[1]] = true;

                foundCandidate = true;

                numFixed++;
            }
        }
        while (foundCandidate);

        std::cout << "Fixed " << numFixed << " / " << chartData.subsides.size() << " subsides. Target was " << numToFix << std::endl;
    }

    //Solve ILP to find best side size
    double gap;
    ilpResult = QuadRetopology::findSubdivisions(
            chartData,
            fixedSubdivisionSubsides,
            chartEdgeLength,
            parameters,
            gap);

    //Quadrangulate
    std::vector<size_t> fixedPositionSubsides;
    std::vector<int> quadmeshLabel;
    QuadRetopology::quadrangulate(
            trimesh,
            chartData,
            fixedPositionSubsides,
            ilpResult,
            parameters,
            quadmesh,
            quadmeshLabel,
            quadmeshPartitions,
            quadmeshCorners);
}

}
