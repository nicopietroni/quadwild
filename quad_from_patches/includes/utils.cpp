#include "utils.h"

#include <vector>

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/complex/algorithms/polygonal_algorithms.h>
#include <vcg/complex/algorithms/geodesic.h>
#include <vcg/space/distance3.h>

namespace qfp {

bool findVertexChainPathRecursive(
        const size_t& vCurrentId,
        const size_t& vStartId,
        const std::vector<std::vector<size_t>>& vertexNextMap,
        std::vector<size_t>& nextConfiguration);

template <class Mesh>
void updateAllMeshAttributes(Mesh &mesh)
{
    vcg::tri::UpdateNormal<Mesh>::PerFaceNormalized(mesh);
    vcg::tri::UpdateNormal<Mesh>::PerVertexNormalized(mesh);
    vcg::tri::UpdateBounding<Mesh>::Box(mesh);
    vcg::tri::UpdateTopology<Mesh>::FaceFace(mesh);
    vcg::tri::UpdateTopology<Mesh>::VertexFace(mesh);
    vcg::tri::UpdateFlags<Mesh>::FaceBorderFromFF(mesh);
    vcg::tri::UpdateFlags<Mesh>::VertexBorderFromNone(mesh);
}

template <class Mesh>
std::vector<std::vector<size_t>> findConnectedComponents(
        const Mesh &mesh)
{
    std::vector<std::vector<size_t>> components;

    std::set<size_t> visited;

    for (size_t i = 0; i < mesh.face.size(); i++) {
        std::stack<size_t> stack;

        //If face has not been visited
        if (visited.find(i) == visited.end()) {
            stack.push(i);

            std::vector<size_t> facesInComponents;
            do {
                size_t fId = stack.top();
                stack.pop();

                if (visited.find(fId) == visited.end()) {
                    facesInComponents.push_back(fId);
                    visited.insert(fId);

                    for (size_t j = 0; j < mesh.face[fId].VN(); j++) {
                        if (!vcg::face::IsBorder(mesh.face[fId], j)) {
                            size_t adjFId = vcg::tri::Index(mesh, mesh.face[fId].cFFp(j));

                            if (visited.find(adjFId) == visited.end()) {
                                stack.push(adjFId);
                            }
                        }
                    }
                }
            } while (!stack.empty());
            components.push_back(facesInComponents);
        }
    }

#ifdef SAVEMESHESFORDEBUG
    //Check unique faces in components
    for (size_t i = 0; i < components.size(); i++)
    {
        std::set<size_t> checkComponentsUnique(components[i].begin(), components[i].end());
        assert(components[i].size() == checkComponentsUnique.size());
    }
#endif

    return components;
}



template<class Mesh>
void LaplacianPos(Mesh &poly_m,std::vector<typename Mesh::CoordType> &AvVert)
{
    //cumulate step
    AvVert.clear();
    AvVert.resize(poly_m.vert.size(),typename Mesh::CoordType(0,0,0));
    std::vector<typename Mesh::ScalarType> AvSum(poly_m.vert.size(),0);
    for (size_t i=0;i<poly_m.face.size();i++)
        for (size_t j=0;j<(size_t)poly_m.face[i].VN();j++)
        {
            //get current vertex
            typename Mesh::VertexType *currV=poly_m.face[i].V(j);
            //and its position
            typename Mesh::CoordType currP=currV->P();
            //cumulate over other positions
            typename Mesh::ScalarType W=vcg::PolyArea(poly_m.face[i]);
            //assert(W!=0);
            for (size_t k=0;k<(size_t)poly_m.face[i].VN();k++)
            {
                if (k==j) continue;
                int IndexV=vcg::tri::Index(poly_m,poly_m.face[i].V(k));
                AvVert[IndexV]+=currP*W;
                AvSum[IndexV]+=W;
            }
        }

    //average step
    for (size_t i=0;i<poly_m.vert.size();i++)
    {
        if (AvSum[i]==0)continue;
        AvVert[i]/=AvSum[i];
    }
}

template <class Mesh>
void LaplacianGeodesic(
        Mesh &poly_m,
        int nstep,
        const double maxDistance,
        const double minDumpS,
        std::vector<size_t>& smoothedVertices)
{
    std::vector<typename Mesh::VertexPointer> seedVec;
    for (int i = 0; i < poly_m.vert.size(); i++) {
        if (poly_m.vert[i].IsS()) {
            seedVec.push_back(&poly_m.vert[i]);
        }
    }
    vcg::tri::UpdateQuality<Mesh>::VertexConstant(poly_m, std::numeric_limits<float>::max());
    vcg::tri::EuclideanDistance<Mesh> ed;
    vcg::tri::UpdateTopology<Mesh>::VertexFace(poly_m);
    vcg::tri::Geodesic<Mesh>::Compute(poly_m,seedVec, ed);

    smoothedVertices.clear();
    std::vector<double> DampS(poly_m.vert.size());
    for (int i = 0; i < poly_m.vert.size(); i++) {        
        if (!poly_m.vert[i].IsD()) {
            if (poly_m.vert[i].Q() < maxDistance) {
                smoothedVertices.push_back(i);
                DampS[i] = poly_m.vert[i].Q() / maxDistance;
                assert(DampS[i] >= 0 && DampS[i] <= 1);
                DampS[i] = minDumpS + DampS[i]*(1-minDumpS);
            }
            else {
                DampS[i] = std::numeric_limits<double>::max();
            }
        }
    }

    for (int s=0;s<nstep;s++)
    {
        std::vector< typename Mesh::CoordType> AvVert;
        LaplacianPos(poly_m,AvVert);

        for (size_t i=0;i<poly_m.vert.size();i++)
        {
            if (!poly_m.vert[i].IsD() && DampS[i] <= 1) {
                poly_m.vert[i].P()=poly_m.vert[i].P()*DampS[i]+
                        AvVert[i]*(1-DampS[i]);
            }
        }
    }
}


template <class Mesh>
bool isTriangleMesh(Mesh& mesh) {
    for (size_t i = 0; i < mesh.face.size(); i++) {
        if (mesh.face[i].VN() != 3)
            return false;
    }

    return true;
}

template <class Mesh>
bool isQuadMesh(Mesh& mesh) {
    for (size_t i = 0; i < mesh.face.size(); i++) {
        if (mesh.face[i].VN() != 4)
            return false;
    }

    return true;
}


inline bool findVertexChainPathRecursive(
        const size_t& vCurrentId,
        const size_t& vStartId,
        const std::vector<std::vector<size_t>>& vertexNextMap,
        std::vector<size_t>& nextConfiguration)
{
    if (vCurrentId == vStartId)
        return true;

    for (size_t i = 0; i < vertexNextMap[vCurrentId].size(); i++) {
        nextConfiguration[vCurrentId] = i;
        if (findVertexChainPathRecursive(vertexNextMap[vCurrentId][i], vStartId, vertexNextMap, nextConfiguration)) {
            return true;
        }
    }
    return false;
}

inline std::vector<size_t> findVertexChainPath(
        const size_t& vCurrentId,
        const std::vector<std::vector<size_t>>& vertexNextMap)
{
    std::vector<size_t> nextConfiguration(vertexNextMap.size());

    findVertexChainPathRecursive(vCurrentId, vCurrentId, vertexNextMap, nextConfiguration);

    return nextConfiguration;
}

}
