#ifndef SMOOTH_MESH_H
#define SMOOTH_MESH_H

#include <vector>
#include <vcg/complex/algorithms/polygonal_algorithms.h>
#include <vcg/complex/algorithms/polygonal_algorithms.h>

/* ----- Triangle mesh ----- */

class BasicVertex;
class BasicEdge;
class BasicFace;

struct MyBasicTypes : public vcg::UsedTypes<
        vcg::Use<BasicVertex>::AsVertexType,
        vcg::Use<BasicEdge>::AsEdgeType,
        vcg::Use<BasicFace>::AsFaceType>{};

class BasicVertex : public vcg::Vertex<
        MyBasicTypes,
        vcg::vertex::Normal3d,
        vcg::vertex::Coord3d,
        vcg::vertex::BitFlags>{};

class BasicEdge : public vcg::Edge<
        MyBasicTypes,
        vcg::edge::VertexRef,
        vcg::edge::BitFlags> {};

class BasicFace : public vcg::Face<
        MyBasicTypes,
        vcg::face::VertexRef,
        vcg::face::BitFlags,
        vcg::face::Mark> {};

class BasicMesh : public vcg::tri::TriMesh<
        std::vector<BasicVertex>,
        std::vector<BasicEdge>,
        std::vector<BasicFace> > {};

template <class TriangleMeshType>
void ExtractEdgeMesh(const TriangleMeshType &triMesh,
                     const std::vector<std::pair<size_t,size_t> > features,
                     BasicMesh &EdgeMesh)
{
    typedef typename TriangleMeshType::FaceType TriFaceType;
    typedef typename TriangleMeshType::ScalarType TriScalarType;
    typedef typename TriangleMeshType::CoordType CoordType;

    EdgeMesh.Clear();
    for (size_t i=0;i<features.size();i++)
    {
        size_t IndexF=features[i].first;
        size_t IndexE=features[i].second;
        CoordType P0=triMesh.face[IndexF].cP0(IndexE);
        CoordType P1=triMesh.face[IndexF].cP1(IndexE);
        vcg::tri::Allocator<BasicMesh>::AddEdge(EdgeMesh,P0,P1);
    }
    vcg::tri::Clean<BasicMesh>::RemoveDuplicateVertex(EdgeMesh);
    vcg::tri::Clean<BasicMesh>::RemoveDuplicateEdge(EdgeMesh);
    vcg::tri::Allocator<BasicMesh>::CompactEveryVector(EdgeMesh);
}

void SelectCorners(BasicMesh &EdgeMesh)
{
    std::vector<size_t> NumE(EdgeMesh.vert.size(),0);
    vcg::tri::UpdateFlags<BasicMesh>::VertexClearS(EdgeMesh);
    for (size_t i=0;i<EdgeMesh.edge.size();i++)
    {
        size_t IndexV0=vcg::tri::Index(EdgeMesh,EdgeMesh.edge[i].V(0));
        size_t IndexV1=vcg::tri::Index(EdgeMesh,EdgeMesh.edge[i].V(1));
        NumE[IndexV0]++;
        NumE[IndexV1]++;
    }
    for (size_t i=0;i<NumE.size();i++)
        if (NumE[i]!=2)EdgeMesh.vert[i].SetS();
}

void ClosestPointEMesh(const typename BasicMesh::CoordType &Pos,
                       const BasicMesh &EdgeMesh,
                       size_t &IndexE,
                       typename BasicMesh::ScalarType &t,
                       typename BasicMesh::ScalarType &MinD,
                       typename BasicMesh::CoordType &Clos)
{
    typedef typename BasicMesh::ScalarType ScalarType;
    typedef typename BasicMesh::CoordType CoordType;
    MinD=std::numeric_limits<ScalarType>::max();
    for (size_t i=0;i<EdgeMesh.edge.size();i++)
    {
        CoordType P0=EdgeMesh.edge[i].V(0)->P();
        CoordType P1=EdgeMesh.edge[i].V(1)->P();
        vcg::Segment3<ScalarType> STest(P0,P1);
        ScalarType testD;
        CoordType ClosTest;
        vcg::SegmentPointDistance(STest,Pos,ClosTest,testD);
        if (testD>MinD)continue;
        Clos=ClosTest;
        MinD=testD;
        IndexE=i;
        t=1-(Clos-P0).Norm()/(P1-P0).Norm();
    }
}

enum VType{Internal,Feature,Corner};

//template <class PolyMeshType>
//void GetVertType(const BasicMesh &EdgeMesh,
//                 const PolyMeshType &polyMesh,
//                 std::vector<VType> &VertTypes)
//{
//    typedef typename BasicMesh::ScalarType ScalarType;
//    typedef typename BasicMesh::CoordType CoordType;

//    VertTypes.clear();
//    VertTypes.resize(polyMesh.vert.size(),Internal);
//    const ScalarType Delta=0.0001;
//    for (size_t i=0;i<polyMesh.vert.size();i++)
//    {
//        size_t IndexE;
//        ScalarType t,minD;
//        CoordType clos;

//        ClosestPoint(polyMesh.vert[i].cP(),EdgeMesh,IndexE,t,minD,clos);
//        if (minD>Delta)continue;
//        BasicVertex *EV0=EdgeMesh.edge[IndexE].V(0);
//        BasicVertex *EV1=EdgeMesh.edge[IndexE].V(1);

////        if ((t>(1-Delta))&&(EV0->IsS()))
////        {
////            VertTypes[i]=Corner;
////            continue;
////        }
////        if ((t<Delta)&&(EV1->IsS()))
////        {
////            VertTypes[i]=Corner;
////            continue;
////        }
//        VertTypes[i]=Feature;
//    }
//}


template <class PolyMeshType,class TriangleMeshType>
void GetVertType(const TriangleMeshType &tri_mesh,
                 const PolyMeshType &poly_mesh,
                 const std::vector<std::pair<size_t,size_t> > FeatureEdges,
                 const std::vector<size_t> &tri_face_partition,
                 const std::vector<size_t> &quad_face_partition,
                 std::vector<std::vector<bool> > &IsSharp,
                 std::vector<VType> &VertTypes)
{
    typedef typename BasicMesh::ScalarType ScalarType;
    typedef typename BasicMesh::CoordType CoordType;


    std::set<std::pair<size_t,size_t> > SharpPart;

    for (size_t i=0;i<FeatureEdges.size();i++)
    {
        size_t IndexF=FeatureEdges[i].first;
        size_t IndexE=FeatureEdges[i].second;
        if (vcg::face::IsBorder(tri_mesh.face[IndexF],IndexE))continue;
        size_t Part0=tri_face_partition[IndexF];

        size_t IndexFOpp=vcg::tri::Index(tri_mesh,tri_mesh.face[IndexF].cFFp(IndexE));
        size_t Part1=tri_face_partition[IndexFOpp];
        assert(Part0!=Part1);
        std::pair<size_t,size_t> key(std::min(Part0,Part1),std::max(Part0,Part1));
        SharpPart.insert(key);
    }

    //then check the edge of the quads
    IsSharp.resize(poly_mesh.face.size());
    for (size_t i=0;i<poly_mesh.face.size();i++)
    {
        size_t Part0=quad_face_partition[i];
        IsSharp[i].resize(poly_mesh.face[i].VN(),false);
        for (size_t j=0;j<poly_mesh.face[i].VN();j++)
        {
            if (vcg::face::IsBorder(poly_mesh.face[i],j))
                IsSharp[i][j]=true;
            else
            {
                size_t IndexFOpp=vcg::tri::Index(poly_mesh,poly_mesh.face[i].cFFp(j));
                size_t Part1=quad_face_partition[IndexFOpp];
                std::pair<size_t,size_t> key(std::min(Part0,Part1),std::max(Part0,Part1));

                if (SharpPart.count(key))
                    IsSharp[i][j]=true;
            }
        }
    }

    //finally set the vertex type
    VertTypes.clear();
    VertTypes.resize(poly_mesh.vert.size(),Internal);
    std::vector<size_t> SharpVal(poly_mesh.vert.size(),0);
    for (size_t i=0;i<poly_mesh.face.size();i++)
    {
        size_t NumV=poly_mesh.face[i].VN();
        for (size_t j=0;j<NumV;j++)
        {
            if (!IsSharp[i][j])continue;
            size_t IndexV0=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V(j));
            size_t IndexV1=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V((j+1)%NumV));
            SharpVal[IndexV0]++;
            SharpVal[IndexV1]++;
        }
    }

    for (size_t i=0;i<SharpVal.size();i++)
    {
        if (SharpVal[i]==0)continue;
        if (SharpVal[i]==4)
            VertTypes[i]=Feature;
        else
            VertTypes[i]=Corner;
    }
}

enum SmoothType{Laplacian,TemplateFit};


template <class PolyMeshType,class TriangleMeshType>
void SmoothWithFeatures(TriangleMeshType &tri_mesh,
                        PolyMeshType &poly_mesh,
                        const std::vector<std::pair<size_t,size_t> > &features,
                        const std::vector<size_t> &tri_face_partition,
                        const std::vector<size_t> &quad_face_partition,
                        SmoothType SType,
                        size_t Steps=10,
                        typename PolyMeshType::ScalarType Damp=0.5)
{
    typedef typename TriangleMeshType::FaceType TriFaceType;
    typedef typename TriangleMeshType::ScalarType TriScalarType;
    typedef typename TriangleMeshType::CoordType CoordType;

    typedef typename PolyMeshType::VertexType PolyVertexType;

    BasicMesh EdgeMesh;
    ExtractEdgeMesh(tri_mesh,features,EdgeMesh);

    std::vector<std::vector<bool> > IsSharp;
    std::vector<VType> VertTypes;
    GetVertType(tri_mesh,poly_mesh,features,tri_face_partition,quad_face_partition,IsSharp,VertTypes);

    vcg::GridStaticPtr<TriFaceType,TriScalarType> TriGrid;
    TriGrid.Set(tri_mesh.face.begin(),tri_mesh.face.end());

    //    vcg::tri::io::ExporterOBJ<BasicMesh>::Save(EdgeMesh,"test_edge.obj",vcg::tri::io::Mask::IOM_EDGEINDEX);
    //    return;
    for (size_t s=0;s<Steps;s++)
    {
        //Smooth only along features
        std::vector<CoordType> TargetPos(poly_mesh.vert.size(),CoordType(0,0,0));
        std::vector<size_t> TargetNum(poly_mesh.vert.size(),0);
        vcg::tri::UpdateSelection<PolyMeshType>::VertexClear(poly_mesh);
        for (size_t i=0;i<poly_mesh.face.size();i++)
        {
            size_t NumV=poly_mesh.face[i].VN();
            for (size_t j=0;j<poly_mesh.face[i].VN();j++)
            {
                if (!IsSharp[i][j])continue;

                PolyVertexType *V0=poly_mesh.face[i].V(j);
                PolyVertexType *V1=poly_mesh.face[i].V((j+1)%NumV);

                size_t IndexV0=vcg::tri::Index(poly_mesh,V0);
                size_t IndexV1=vcg::tri::Index(poly_mesh,V1);

                poly_mesh.vert[IndexV0].SetS();
                poly_mesh.vert[IndexV1].SetS();

                CoordType PosV0=poly_mesh.vert[IndexV0].P();
                CoordType PosV1=poly_mesh.vert[IndexV1].P();
                TargetPos[IndexV0]+=PosV1;
                TargetPos[IndexV1]+=PosV0;
                TargetNum[IndexV0]++;
                TargetNum[IndexV1]++;
            }
        }
        for (size_t i=0;i<poly_mesh.vert.size();i++)
        {
            if (VertTypes[i]!=Feature)continue;
            if (TargetNum[i]==0)continue;

            //then average
            CoordType TargetP=TargetPos[i]/TargetNum[i];
            poly_mesh.vert[i].P()=poly_mesh.vert[i].P()*Damp+TargetP*(1-Damp);

            //and reproject
            size_t IndexE;
            TriScalarType t,MinD;
            CoordType closestPt;
            ClosestPointEMesh(poly_mesh.vert[i].P(),EdgeMesh,IndexE,t,MinD,closestPt);
            poly_mesh.vert[i].P()=closestPt;
        }

        if (SType==Laplacian)
            vcg::PolygonalAlgorithm<PolyMeshType>::Laplacian(poly_mesh,true,1,Damp);
        else
            vcg::PolygonalAlgorithm<PolyMeshType>::SmoothPCA(poly_mesh,1,Damp,true,true,0.5,false);

        //then reproject
        for (size_t i=0;i<poly_mesh.vert.size();i++)
        {
            if (VertTypes[i]!=Internal)continue;

            TriScalarType MaxD=poly_mesh.bbox.Diag();
            TriScalarType MinD;
            CoordType closestPt;
            vcg::tri::GetClosestFaceBase(tri_mesh,TriGrid,poly_mesh.vert[i].P(),MaxD,MinD,closestPt);
            poly_mesh.vert[i].P()=closestPt;
        }
    }
    //    vcg::tri::io::ExporterOBJ<BasicMesh>::Save(EdgeMesh,"test_edge.obj",vcg::tri::io::Mask::IOM_EDGEINDEX);
    //return;
}

#endif // SMOOTH_MESH_H
