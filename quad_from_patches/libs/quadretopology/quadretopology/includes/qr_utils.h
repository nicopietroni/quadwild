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

#ifndef QR_UTILS_H
#define QR_UTILS_H

#include <vector>
#include "qr_charts.h"

namespace QuadRetopology {
namespace internal {

template <class PolyMeshType>
int removeUnreferencedVertices(PolyMeshType& m, bool onlySelected);

template <class PolyMeshType, class ScalarType>
int clusterVertices(PolyMeshType &m, const bool onlySelected, const ScalarType radius);

template <class PolyMeshType>
int removeDuplicateVertices(PolyMeshType& mesh, const bool onlySelected);

template <class PolyMeshType>
int removeDegenerateFaces(PolyMeshType& mesh, bool onlySelected, bool alwaysDelete);

template <class PolyMeshType>
int removeDoubletFaces(PolyMeshType& mesh, bool onlySelected, bool recursive);

static std::vector<size_t> dummySizetVector;

template <class MeshType>
void updateAllMeshAttributes(MeshType &mesh);

template <class MeshType>
std::vector<std::vector<size_t>> findConnectedComponents(
        const MeshType &mesh);

template<class MeshType>
void LaplacianPos(MeshType &poly_m,std::vector<typename MeshType::CoordType> &AvVert);

template <class MeshType>
void LaplacianGeodesicSmoothing(
        MeshType &poly_m,
        int nstep,
        const double maxDistance,
        const double minDumpS = 0.5,
        std::vector<size_t>& smoothedVertices = dummySizetVector);

std::vector<size_t> findVertexChainPath(
        const size_t& vCurrentId,
        const std::vector<std::vector<size_t>>& vertexNextMap);


template <class MeshType>
bool isTriangleMesh(MeshType& mesh);

template <class MeshType>
bool isQuadMesh(MeshType& mesh);

template <class MeshType>
std::vector<int> splitFacesInTriangles(MeshType& mesh);

template<class PolyMeshType>
typename PolyMeshType::ScalarType averageEdgeLength(PolyMeshType& mesh, const std::vector<size_t>& faces);

template<class PolyMeshType>
typename PolyMeshType::ScalarType averageEdgeLength(PolyMeshType& mesh);


}
}

#include "qr_utils.cpp"

#endif // QR_UTILS_H
