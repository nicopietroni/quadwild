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

#ifndef QUADFROMPATCHES_H
#define QUADFROMPATCHES_H

#include <quadretopology/quadretopology.h>

namespace qfp {

template<class PolyMesh, class TriangleMesh>
void quadrangulationFromPatches(
    TriangleMesh& trimesh,
    const std::vector<std::vector<size_t>>& trimeshPartitions,
    const std::vector<std::vector<size_t>>& trimeshCorners,
    const std::vector<double>& chartEdgeLength,
    const QuadRetopology::Parameters& parameters,
    const int fixedChartClusters,
    PolyMesh& quadmesh,
    std::vector<std::vector<size_t>>& quadmeshPartitions,
    std::vector<std::vector<size_t>>& quadmeshCorners,
    std::vector<int>& ilpResult);


}

#include "quad_from_patches.cpp"

#endif // QUADFROMPATCHES_H
