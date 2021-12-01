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

#include "qr_charts.h"

#include "assert.h"

#include <map>
#include <unordered_set>

#include "qr_utils.h"

#include <vcg/complex/complex.h>
#include <wrap/io_trimesh/export.h>

#define MAXITERATIONS 100000

namespace QuadRetopology {
namespace internal {

template<class TriangleMeshType>
void findChartFacesAndBorderFaces(
        TriangleMeshType& mesh,
        const std::vector<int>& faceLabel,
        ChartData& chartData)
{
    if (mesh.face.size() == 0)
        return;

    std::set<int>& labels = chartData.labels;
    std::vector<Chart>& charts = chartData.charts;

    for (size_t fId = 0; fId < mesh.face.size(); fId++) {
        if (!mesh.face[fId].IsD() && faceLabel.at(fId) >= 0)
            labels.insert(faceLabel.at(fId));
    }

    int maxChartLabel = *labels.rbegin();
    charts.resize(maxChartLabel + 1);

    std::vector<bool> visited(mesh.face.size(), false);

    for (size_t i = 0; i < mesh.face.size(); i++) {
        if (!mesh.face[i].IsD() && !visited[i]) {
            //Region growing to get chart
            std::stack<size_t> stack;
            stack.push(i);

            Chart chart;

            chart.label = faceLabel[i];

            if (chart.label >= 0) {
                std::set<size_t> borderFacesSet;
                do {
                    size_t fId = stack.top();
                    stack.pop();
                    assert(faceLabel[fId] == chart.label);

                    if (!visited[fId]) {
                        chart.faces.push_back(fId);

                        typename TriangleMeshType::FaceType* currentFacePointer = &mesh.face[fId];
                        vcg::face::Pos<typename TriangleMeshType::FaceType> pos(currentFacePointer, 0);
                        for (int k = 0; k < 3; k++) {
                            pos.FlipF();
                            size_t adjFace = vcg::tri::Index(mesh, pos.F());

                            //Saving border faces
                            if (currentFacePointer == pos.F() || faceLabel[adjFace] != chart.label) {
                                borderFacesSet.insert(fId);
                            }
                            else {
                                if (!visited[adjFace]) {
                                    stack.push(adjFace);
                                }
                            }
                            pos.FlipF();

                            //Next edge
                            pos.FlipV();
                            pos.FlipE();
                        }

                        visited[fId] = true;
                    }
                }
                while (!stack.empty());

                assert(borderFacesSet.size() >= 1);

                std::copy(borderFacesSet.begin(), borderFacesSet.end(), std::back_inserter(chart.borderFaces));

                charts[chart.label] = chart;
            }

        }
    }
}

}
}
