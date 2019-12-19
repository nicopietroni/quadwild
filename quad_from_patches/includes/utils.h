#ifndef QUADFROMPATCHES_UTILS_H
#define QUADFROMPATCHES_UTILS_H

#include <vector>

namespace qfp {
static std::vector<size_t> dummySizetVector;

template <class Mesh>
void updateAllMeshAttributes(Mesh &mesh);

template <class Mesh>
std::vector<std::vector<size_t>> findConnectedComponents(
        const Mesh &mesh);

template<class Mesh>
void LaplacianPos(Mesh &poly_m,std::vector<typename Mesh::CoordType> &AvVert);

template <class Mesh>
void LaplacianGeodesic(
        Mesh &poly_m,
        int nstep,
        const double maxDistance,
        const double minDumpS = 0.5,
        std::vector<size_t>& smoothedVertices = dummySizetVector);

std::vector<size_t> findVertexChainPath(
        const size_t& vCurrentId,
        const std::vector<std::vector<size_t>>& vertexNextMap);

template <class Mesh>
bool isTriangleMesh(Mesh& mesh);

template <class Mesh>
bool isQuadMesh(Mesh& mesh);

}

#include "utils.cpp"

#endif // QUADFROMPATCHES_UTILS_H
