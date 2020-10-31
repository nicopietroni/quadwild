#include "qr_utils.h"

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

namespace QuadRetopology {
namespace internal {

template <class PolyMeshType>
int removeUnreferencedVertices(PolyMeshType& m, bool onlySelected)
{
    int deleted = 0;

    std::vector<bool> referredVec(m.vert.size(),false);

    for(auto fi = m.face.begin(); fi != m.face.end(); ++fi) {
        if(fi->IsD())
            continue;

        for(int j=0; j < fi->VN(); ++j) {
            referredVec[vcg::tri::Index(m, fi->V(j))]=true;
        }
    }

    for(auto vi=m.vert.begin();vi!=m.vert.end();++vi) {
        if(vi->IsD())
            continue;

        if(onlySelected && !vi->IsS())
            continue;

        if(!referredVec[vcg::tri::Index(m,*vi)]) {
            vcg::tri::Allocator<PolyMeshType>::DeleteVertex(m,*vi);
            ++deleted;
        }
    }

    return deleted;
}

template <class PolyMeshType, class ScalarType>
int clusterVertices(PolyMeshType &mesh, const bool onlySelected, const ScalarType radius)
{
    typedef typename PolyMeshType::CoordType CoordType;
    typedef typename PolyMeshType::VertexPointer VertexPointer;

    size_t nVertices = mesh.vert.size();

    if(nVertices == 0)
        return 0;

    int clustered = 0;

    std::vector<VertexPointer> vertices;

    for(size_t i = 0; i < mesh.vert.size(); i++) {
        if (mesh.vert[i].IsD())
            continue;

        vertices.push_back(&mesh.vert[i]);
    }

    std::sort(vertices.begin(), vertices.end(), [] (VertexPointer const &a, VertexPointer const &b)
    {
        return a->cP() < b->cP();
    });

    for(size_t i = 0; i < vertices.size(); i++) {
        if (vertices[i]->IsD())
            continue;

        const CoordType& p = vertices[i]->cP();

        std::vector<VertexPointer> cluster;

        cluster.push_back(vertices[i]);

        bool done = false;
        for(size_t j = i+1; !done && j < vertices.size(); j++) {
            if (vertices[j]->IsD())
                continue;

            const CoordType& q = vertices[j]->cP();

            double distanceX = std::abs(q.X() - p.X());
            double distanceY = std::abs(q.Y() - p.Y());
            double distanceZ = std::abs(q.Z() - p.Z());

            if (distanceX > radius && distanceY > radius && distanceZ > radius) {
                done = true;
            }
            else if ((p - q).Norm() <= radius) {
                cluster.push_back(vertices[j]);
            }
        }

        if (cluster.size() > 1) {
            CoordType point = cluster[0]->cP();

            if (onlySelected) {
                for(size_t j = 0; j < cluster.size(); j++) {
                    if (!cluster[j]->IsS()) {
                        point = cluster[j]->cP();
                        break;
                    }
                }
            }

            for(size_t j = 0; j < cluster.size(); j++) {
                if (cluster[j]->cP() != point && (!onlySelected || cluster[j]->IsS())) {
                    cluster[j]->P() = point;

                    clustered++;
                }
            }
        }
    }

    return clustered;
}

template <class PolyMeshType>
int removeDegenerateFaces(PolyMeshType& mesh, bool onlySelected, bool alwaysDelete)
{
    typedef typename PolyMeshType::FaceIterator FaceIterator;
    typedef typename PolyMeshType::VertexPointer VertexPointer;
    int count_fd = 0;

    for(FaceIterator fi = mesh.face.begin(); fi != mesh.face.end(); ++fi) {
        if(fi->IsD())
            continue;

        if(onlySelected && !fi->IsS())
            continue;

        int numRemainingVertices = fi->VN();

        for (int j1 = 0; j1 < numRemainingVertices - 1; ++j1) {
            for (int j2 = j1 + 1; j2 < numRemainingVertices; ++j2) {
                if (fi->V(j1) == fi->V(j2)) {
                    for (int k = j2; k < numRemainingVertices - 1; ++k) {
                        fi->V0(k) = fi->V1(k);
                    }

                    numRemainingVertices--;

                    j2--;
                }
            }
        }

        if (numRemainingVertices < fi->VN()) {
            if(alwaysDelete || numRemainingVertices < 3) {
                Allocator<PolyMeshType>::DeleteFace(mesh,*fi);
            }
            else {
                std::vector<VertexPointer> vp(numRemainingVertices);

                for (int j = 0; j < numRemainingVertices; j++) {
                    vp[j] = fi->V(j);
                }

                fi->Dealloc();
                fi->Alloc(numRemainingVertices);

                for (int j = 0; j < numRemainingVertices; j++) {
                    fi->V(j) = vp[j];
                }
            }

            count_fd++;
        }
    }

    return count_fd;
}

template <class PolyMeshType>
int removeDuplicateVertices(PolyMeshType& mesh, const bool onlySelected)
{
    typedef typename PolyMeshType::CoordType CoordType;
    typedef typename PolyMeshType::VertexPointer VertexPointer;
    typedef typename PolyMeshType::VertexPointer VertexIterator;
    typedef typename PolyMeshType::FaceIterator FaceIterator;

    int deleted = 0;

    std::vector<VertexPointer> vertices;

    for(size_t i = 0; i < mesh.vert.size(); i++) {
        if (mesh.vert[i].IsD())
            continue;

        vertices.push_back(&mesh.vert[i]);
    }

    std::sort(vertices.begin(), vertices.end(), [] (VertexPointer const &a, VertexPointer const &b)
    {
        return a->cP() < b->cP();
    });

    std::map<VertexPointer, VertexPointer> vertexMap;
    for(size_t i = 0; i < vertices.size(); i++) {
        if (vertices[i]->IsD())
            continue;

        const CoordType& p = vertices[i]->cP();

        std::vector<VertexPointer> cluster;

        cluster.push_back(vertices[i]);

        bool done = false;
        for(size_t j = i+1; !done && j < vertices.size(); j++) {
            if (vertices[j]->IsD())
                continue;

            const CoordType& q = vertices[j]->cP();

            if (p == q) {
                cluster.push_back(vertices[j]);
            }
            else {
                done = true;
            }
        }

        if (cluster.size() > 1) {
            size_t kept = 0;

            if (onlySelected) {
                for(size_t j = 0; j < cluster.size(); j++) {
                    if (!cluster[j]->IsS()) {
                        kept = j;
                        break;
                    }
                }
            }

            for(size_t j = 0; j < cluster.size(); j++) {
                if (j != kept && (!onlySelected || cluster[j]->IsS())) {
                    Allocator<PolyMeshType>::DeleteVertex(mesh, *cluster[j]);
                    deleted++;

                    vertexMap[cluster[j]] = cluster[kept];
                }
            }
        }
    }

    for(FaceIterator fi = mesh.face.begin(); fi != mesh.face.end(); ++fi) {
        if(fi->IsD())
            continue;

        for (int k = 0; k < fi->VN(); ++k) {
            typename std::map<VertexPointer, VertexPointer>::iterator it = vertexMap.find(fi->V(k));
            if (it != vertexMap.end()) {
                fi->V(k) = it->second;
            }
        }
    }

    return deleted;
}

template <class PolyMeshType>
int removeDoublets(PolyMeshType& mesh, bool onlySelectedVertices)
{
    int count = 0;

    vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(mesh);

    for (size_t i = 0; i < mesh.face.size(); ++i) {
        if (mesh.face[i].IsD() || mesh.face[i].VN() != 4)
            continue;

        typename PolyMeshType::FacePointer face = &mesh.face[i];

        for (int j1 = 0; j1 < face->VN(); ++j1) {
            int j2 = (j1+1) % face->VN();
            int j3 = (j1+2) % face->VN();
            int j4 = (j1+3) % face->VN();

            if (onlySelectedVertices && (!face->V(j1)->IsS() || !face->V(j2)->IsS() || !face->V(j3)->IsS()))
                continue;

            if (vcg::face::IsBorder(*face, j1) || vcg::face::IsBorder(*face, j2))
                continue;

            if (face->FFp(j1) == face->FFp(j2)) {
                typename PolyMeshType::FacePointer adjFace = face->FFp(j2);

                if (adjFace->IsD() || adjFace->VN() != 4)
                    continue;

                int adjJ1 = adjFace->FFi(j1);
                int adjJ2 = (adjJ1+1) % adjFace->VN();
                int adjJ3 = (adjJ1+2) % adjFace->VN();
                int adjJ4 = (adjJ1+3) % adjFace->VN();

                //Remove warnings
                static_cast<void>(j3);
                static_cast<void>(j4);
                static_cast<void>(adjJ2);

                if (adjFace->FFp(adjJ1) != face || adjFace->FFp(adjJ4) != face) //Non manifoldness (it should not happen)
                    continue;

                assert(face->V(j1) == adjFace->V(adjJ2) && face->V(j2) == adjFace->V(adjJ1) && face->V(j3) == adjFace->V(adjJ4));

                face->V(j2) = adjFace->V(adjJ3);

                vcg::tri::Allocator<PolyMeshType>::DeleteFace(mesh, *adjFace);
                count++;
            }
        }
    }
    return count;
}

bool findVertexChainPathRecursive(
        const size_t& vCurrentId,
        const size_t& vStartId,
        const std::vector<std::vector<size_t>>& vertexNextMap,
        std::vector<size_t>& nextConfiguration);


template <class MeshType>
void updateAllMeshAttributes(MeshType &mesh)
{
    vcg::tri::UpdateNormal<MeshType>::PerFaceNormalized(mesh);
    vcg::tri::UpdateNormal<MeshType>::PerVertexNormalized(mesh);
    vcg::tri::UpdateBounding<MeshType>::Box(mesh);
    vcg::tri::UpdateTopology<MeshType>::FaceFace(mesh);
    vcg::tri::UpdateTopology<MeshType>::VertexFace(mesh);
    vcg::tri::UpdateFlags<MeshType>::FaceBorderFromFF(mesh);
    vcg::tri::UpdateFlags<MeshType>::VertexBorderFromNone(mesh);
}

template <class MeshType>
std::vector<std::vector<size_t>> findConnectedComponents(
        const MeshType &mesh)
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

#ifndef NDEBUG
    //Check unique faces in components
    for (size_t i = 0; i < components.size(); i++)
    {
        std::set<size_t> checkComponentsUnique(components[i].begin(), components[i].end());
        assert(components[i].size() == checkComponentsUnique.size());
    }
#endif

    return components;
}



template<class MeshType>
void LaplacianPos(MeshType &poly_m,std::vector<typename MeshType::CoordType> &AvVert)
{
    //cumulate step
    AvVert.clear();
    AvVert.resize(poly_m.vert.size(),typename MeshType::CoordType(0,0,0));
    std::vector<typename MeshType::ScalarType> AvSum(poly_m.vert.size(),0);
    for (size_t i=0;i<poly_m.face.size();i++)
        for (size_t j=0;j<(size_t)poly_m.face[i].VN();j++)
        {
            //get current vertex
            typename MeshType::VertexType *currV=poly_m.face[i].V(j);
            //and its position
            typename MeshType::CoordType currP=currV->P();
            //cumulate over other positions
            typename MeshType::ScalarType W=vcg::PolyArea(poly_m.face[i]);
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

template <class MeshType>
void LaplacianGeodesic(
        MeshType &poly_m,
        int nstep,
        const double maxDistance,
        const double minDumpS,
        std::vector<size_t>& smoothedVertices)
{
    std::vector<typename MeshType::VertexPointer> seedVec;
    for (int i = 0; i < poly_m.vert.size(); i++) {
        if (poly_m.vert[i].IsS()) {
            seedVec.push_back(&poly_m.vert[i]);
        }
    }
    vcg::tri::UpdateQuality<MeshType>::VertexConstant(poly_m, std::numeric_limits<float>::max());
    vcg::tri::EuclideanDistance<MeshType> ed;
    vcg::tri::UpdateTopology<MeshType>::VertexFace(poly_m);
    vcg::tri::Geodesic<MeshType>::Compute(poly_m,seedVec, ed);

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
        std::vector< typename MeshType::CoordType> AvVert;
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


template <class MeshType>
std::vector<int> splitFacesInTriangles(MeshType& mesh) {
    typedef typename MeshType::VertexType VertexType;

    size_t numFaces = mesh.face.size();
    std::vector<int> birthFace(mesh.face.size(), -1);

    for (size_t i = 0; i < numFaces; i++) {
        size_t prevSize = mesh.face.size();

        if (mesh.face[i].VN() == 4) {
            VertexType* v[4];
            v[0] = mesh.face[i].V(0);
            v[1] = mesh.face[i].V(1);
            v[2] = mesh.face[i].V(2);
            v[3] = mesh.face[i].V(3);

            size_t startIndex = 0;
            if ((v[0]->P() - v[2]->P()).Norm() > (v[1]->P() - v[3]->P()).Norm()) {
                startIndex++;
            }

            vcg::tri::Allocator<MeshType>::AddFaces(mesh,1);
            mesh.face.back().Alloc(3);
            mesh.face.back().V(0) = v[startIndex];
            mesh.face.back().V(1) = v[(startIndex + 1)%4];
            mesh.face.back().V(2) = v[(startIndex + 2)%4];

            mesh.face[i].Dealloc();
            mesh.face[i].Alloc(3);
            mesh.face[i].V(0)=v[(startIndex + 2)%4];
            mesh.face[i].V(1)=v[(startIndex + 3)%4];
            mesh.face[i].V(2)=v[(startIndex + 4)%4];
        }
        else if (mesh.face[i].VN() > 4) {
            vcg::PolygonalAlgorithm<MeshType>::Triangulate(mesh, i);
        }

        birthFace[i] = static_cast<int>(i);

        for (size_t j = 0; j < mesh.face.size() - prevSize; j++)
            birthFace.push_back(static_cast<int>(i));
    }

    return birthFace;
}

template<class MeshType>
typename MeshType::ScalarType averageEdgeLength(MeshType& mesh, const std::vector<size_t>& faces) {
    typename MeshType::ScalarType length = 0;
    size_t numEdges = 0;
    for (size_t fId : faces) {
        for (size_t j=0;j<mesh.face[fId].VN();j++)
        {
            size_t index0=vcg::tri::Index(mesh,mesh.face[fId].V0(j));
            size_t index1=vcg::tri::Index(mesh,mesh.face[fId].V1(j));
            typename MeshType::CoordType p0=mesh.vert[index0].P();
            typename MeshType::CoordType p1=mesh.vert[index1].P();
            length += (p0-p1).Norm();
            numEdges++;
        }

    }

    length /= numEdges;

    return length;
}

template<class MeshType>
typename MeshType::ScalarType averageEdgeLength(MeshType& mesh) {
    std::vector<size_t> faces;

    for (size_t i=0;i<mesh.face.size();i++) {
        if (!mesh.face[i].IsD()) {
            faces.push_back(i);
        }
    }

    return averageEdgeLength(mesh, faces);
}

template <class MeshType>
class OrientFaces
{
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;

public:

    static void InvertFace(FaceType &f0)
    {
        std::vector<VertexType*> faceVert;
        for (int i=0;i<f0.VN();i++)
            faceVert.push_back(f0.V(i));
        std::reverse(faceVert.begin(),faceVert.end());
        for (int i=0;i<f0.VN();i++)
            f0.V(i)=faceVert[i];
    }

private:

    static bool IsCoherent(const FaceType &f0,const FaceType &f1,const int IndexE)
    {
        assert(f0.cFFp(IndexE)==&f1);
        assert(&f0!=&f1);
        const VertexType *v0=f0.cV(IndexE);
        const VertexType *v1=f0.cV((IndexE+1)%f0.VN());
        int IndexEopp=f0.cFFi(IndexE);
        assert(f1.cFFp(IndexEopp)==&f0);
        const VertexType *v2=f1.cV(IndexEopp);
        const VertexType *v3=f1.cV((IndexEopp+1)%f1.VN());
        if (v0==v2){assert(v1==v3);return false;}
        assert(v0==v3);
        assert(v1==v2);
        return true;
    }

    static void InvertFaces(MeshType &PolyM,
                            const std::vector<int> &ToInvert)
    {
        for (size_t i=0;i<ToInvert.size();i++)
            InvertFace(PolyM.face[ToInvert[i]]);
    }

    static void PropagateFrom(MeshType &PolyM,int &fI0,
                       std::vector<int> &OrientSet0,
                       std::vector<int> &OrientSet1)
    {
        OrientSet0.clear();
        OrientSet1.clear();

        std::vector<int> CoherentSet(PolyM.face.size(),-1);
        CoherentSet[fI0]=0;

        assert(!PolyM.face[fI0].IsS());

        std::vector<int> exploreStack;
        exploreStack.push_back(fI0);
        do{
            //get from the stack and set as explored
            int currF=exploreStack.back();
            exploreStack.pop_back();
            if(PolyM.face[currF].IsS())continue;//already explored
            PolyM.face[currF].SetS();

            //put in the right coherent set
            int currSet=CoherentSet[currF];
            if (currSet==0)
                OrientSet0.push_back(currF);
            else
                OrientSet1.push_back(currF);
            for (int i=0;i<PolyM.face[currF].VN();i++)
            {
                FaceType *f0=&PolyM.face[currF];
                FaceType *f1=PolyM.face[currF].FFp(i);
                if (f1==f0)continue;//border

                int IndexF1=vcg::tri::Index(PolyM,f1);
                if(PolyM.face[IndexF1].IsS())continue;

                exploreStack.push_back(IndexF1);//add to the stack

                //either is coherent or the opposite
                if (IsCoherent(*f0,*f1,i))
                    CoherentSet[IndexF1]=currSet;
                else
                    CoherentSet[IndexF1]=(currSet+1)%2;
            }
        }while(!exploreStack.empty());
    }

public:

    static void AutoOrientFaces(MeshType &PolyM)
    {
        vcg::tri::UpdateTopology<MeshType>::FaceFace(PolyM);
        vcg::tri::UpdateFlags<MeshType>::FaceClearS(PolyM);
        for (int i=0;i<(int)PolyM.face.size();i++)
        {
            if (PolyM.face[i].IsS())continue;
            std::cout<<"Reoriented face."<<std::endl;
            std::vector<int> OrientSet0,OrientSet1;
            PropagateFrom(PolyM,i,OrientSet0,OrientSet1);
            if (OrientSet0.size()<OrientSet1.size())
                InvertFaces(PolyM,OrientSet0);
            else
                InvertFaces(PolyM,OrientSet1);
        }
    }

};


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


template <class MeshType>
bool isTriangleMesh(MeshType& mesh) {
    for (size_t i = 0; i < mesh.face.size(); i++) {
        if (mesh.face[i].VN() != 3)
            return false;
    }

    return true;
}

template <class MeshType>
bool isQuadMesh(MeshType& mesh) {
    for (size_t i = 0; i < mesh.face.size(); i++) {
        if (mesh.face[i].VN() != 4)
            return false;
    }

    return true;
}

}
}
