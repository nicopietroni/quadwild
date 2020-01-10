#ifndef SMOOTH_MESH_H
#define SMOOTH_MESH_H

#include <vector>
#include <vcg/complex/algorithms/polygonal_algorithms.h>
#include <wrap/io_trimesh/export.h>

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
void GetVertType(const BasicMesh &edge_mesh,
                 const typename PolyMeshType::ScalarType &AvEdge,
                 const TriangleMeshType &tri_mesh,
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

        size_t IndexE;
        BasicMesh::ScalarType t,MinD;
        BasicMesh::CoordType Clos;

        ClosestPointEMesh(poly_mesh.vert[i].cP(),edge_mesh,IndexE,t,MinD,Clos);

        if (MinD>AvEdge/10)continue;

        if (SharpVal[i]==4)
            VertTypes[i]=Feature;
        else
            VertTypes[i]=Corner;
    }
}

//template <class PolyMeshType,class TriangleMeshType>
//void GetCorners(const std::vector<std::pair<size_t,size_t> > &FeatureEdges,
//                const TriangleMeshType &tri_mesh,
//                const PolyMeshType &poly_mesh,
//                const std::vector<size_t> &tri_face_partition,
//                const std::vector<size_t> &quad_face_partition,
//                std::vector<size_t> &cornersIdx)
//{
//    cornersIdx.clear();
//    typedef typename BasicMesh::ScalarType ScalarType;
//    typedef typename BasicMesh::CoordType CoordType;
//    typedef typename vcg::Segment3<typename PolyMeshType::ScalarType> SegmentType;

//    //first find the corners
//    std::vector<size_t> Val(tri_mesh.vert.size(),0);

//    for (size_t i=0;i<FeatureEdges.size();i++)
//    {
//        size_t IndexF=FeatureEdges[i].first;
//        size_t IndexE=FeatureEdges[i].second;
//        int IndV0=vcg::tri::Index(tri_mesh,tri_mesh.face[IndexF].cV0(IndexE));
//        int IndV1=vcg::tri::Index(tri_mesh,tri_mesh.face[IndexF].cV1(IndexE));
//        Val[IndV0]++;
//        Val[IndV1]++;
//    }

//    std::vector<std::vector<size_t> > VertPartitions;
//    VertPartitions.resize(tri_mesh.vert.size());
//    assert(tri_face_partition.size()==tri_mesh.face.size());
//    for (size_t i=0;i<tri_face_partition.size();i++)
//        for (size_t j=0;j<3;j++)
//        {
//            int IndV=vcg::tri::Index(tri_mesh,tri_mesh.face[i].cV(j));
//            VertPartitions[IndV].push_back(tri_face_partition[i]);
//        }

//    for (size_t i=0;i<VertPartitions.size();i++)
//    {
//        std::sort(VertPartitions[i].begin(),VertPartitions[i].end());
//        std::vector<size_t>::iterator it;
//        it = std::unique (VertPartitions[i].begin(),VertPartitions[i].end());
//        VertPartitions[i].resize( std::distance(VertPartitions[i].begin(),it) );
//    }

//    //finally find the partitions for the corners
//    std::set<std::vector<size_t> > CornersPart;
//    for (size_t i=0;i<Val.size();i++)
//    {
//        if (Val[i]<3)continue;
//        assert(VertPartitions[i].size()>2);
//        CornersPart.insert(VertPartitions[i]);
//    }
//    //std::cout<<"THERE ARE "<<CornersPart.size()<<"Partition"<<std::endl;

//    std::vector<std::vector<size_t> > PolyVertPartitions;
//    PolyVertPartitions.resize(poly_mesh.vert.size());
//    assert(quad_face_partition.size()==poly_mesh.face.size());
//    for (size_t i=0;i<quad_face_partition.size();i++)
//        for (size_t j=0;j<poly_mesh.face[i].VN();j++)
//        {
//            int IndV=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V(j));
//            PolyVertPartitions[IndV].push_back(quad_face_partition[i]);
//        }

//    for (size_t i=0;i<PolyVertPartitions.size();i++)
//    {
//        std::sort(PolyVertPartitions[i].begin(),PolyVertPartitions[i].end());
//        std::vector<size_t>::iterator it;
//        it = std::unique (PolyVertPartitions[i].begin(),PolyVertPartitions[i].end());
//        PolyVertPartitions[i].resize( std::distance(PolyVertPartitions[i].begin(),it) );
//    }

//    //vcg::tri::UpdateSelection<PolyMeshType>::VertexClear(poly_mesh);
//    for (size_t i=0;i<PolyVertPartitions.size();i++)
//        if (CornersPart.count(PolyVertPartitions[i])>0)
//            cornersIdx.push_back(i);
//    //            std::cout<<"AAAAAA"<<std::endl;
//    //            poly_mesh.vert[i].SetS();
//}

template <class PolyMeshType,class TriangleMeshType>
void GetCorners(const std::vector<std::pair<size_t,size_t> > &FeatureEdges,
                const std::vector<size_t > &tri_feature_C,
                const TriangleMeshType &tri_mesh,
                const std::vector<size_t > &quad_corner,
                const PolyMeshType &poly_mesh,
                const typename PolyMeshType::ScalarType &AvEdge,
                std::vector<size_t> &cornersIdx)
{
    cornersIdx.clear();
    typedef typename BasicMesh::ScalarType ScalarType;
    typedef typename BasicMesh::CoordType CoordType;
    //    typedef typename vcg::Segment3<typename PolyMeshType::ScalarType> SegmentType;


    //    //first find the one ending in nothind that must be disabled
    //    std::vector<size_t> Val(tri_mesh.vert.size(),0);

    //    for (size_t i=0;i<FeatureEdges.size();i++)
    //    {
    //        size_t IndexF=FeatureEdges[i].first;
    //        size_t IndexE=FeatureEdges[i].second;
    //        int IndV0=vcg::tri::Index(tri_mesh,tri_mesh.face[IndexF].cV0(IndexE));
    //        int IndV1=vcg::tri::Index(tri_mesh,tri_mesh.face[IndexF].cV1(IndexE));
    //        Val[IndV0]++;
    //        Val[IndV1]++;
    //    }
    //    std::set<std::vector<size_t> > CornersEnd;
    //    for (size_t i=0;i<Val.size();i++)
    //    {
    //        if (tri_mesh.vert[i].IsB())continue;
    //        if (Val[i]==1)CornersEnd.insert(Val);
    //    }

    for (size_t i=0;i<tri_feature_C.size();i++)
    {
        CoordType TriCorner=tri_mesh.vert[tri_feature_C[i]].P();
        ScalarType MinD=std::numeric_limits<ScalarType>::max();
        size_t CornerI=0;
        for (size_t j=0;j<quad_corner.size();j++)
        {
            CoordType TestPos=poly_mesh.vert[quad_corner[j]].P();
            ScalarType TestD=(TriCorner-TestPos).Norm();
            if (TestD>MinD)continue;
            MinD=TestD;
            CornerI=quad_corner[j];
        }
        if (MinD>AvEdge)continue;
        cornersIdx.push_back(CornerI);
    }
}

template<class ScalarType>
void Merge(const int &IndexMerge0,
           const int &IndexMerge1,
           std::vector<int> &VertBasis,
           std::vector<std::vector<vcg::Segment3<ScalarType> > > &ProjBasis)
{
    assert(IndexMerge0>=0);
    assert(IndexMerge1>=0);
    assert(IndexMerge0<ProjBasis.size());
    assert(IndexMerge1<ProjBasis.size());
    ProjBasis[IndexMerge0].insert(ProjBasis[IndexMerge0].end(),
                                  ProjBasis[IndexMerge1].begin(),
                                  ProjBasis[IndexMerge1].end());

    ProjBasis[IndexMerge1].clear();
    for (size_t i=0;i<VertBasis.size();i++)
        if (VertBasis[i]==IndexMerge1)
            VertBasis[i]=IndexMerge0;
}

template<class TriangleMeshType>
bool IsMergable(const std::vector<vcg::Segment3<typename TriangleMeshType::ScalarType> > &SegBase0,
                const std::vector<vcg::Segment3<typename TriangleMeshType::ScalarType> > &SegBase1,
                const std::set<typename TriangleMeshType::CoordType> &CornerPos)
{
    typedef typename TriangleMeshType::ScalarType ScalarType;
    typedef typename TriangleMeshType::CoordType CoordType;

    for (size_t i=0;i<SegBase0.size();i++)
    {
        CoordType P0[2];
        P0[0]=SegBase0[i].P0();
        P0[1]=SegBase0[i].P1();
        for (size_t j=0;j<SegBase1.size();j++)
        {
            CoordType P1[2];
            P1[0]=SegBase1[j].P0();
            P1[1]=SegBase1[j].P1();
            if ((P0[0]==P1[0])&&(CornerPos.count(P0[0])==0))
                return true;

            if ((P0[0]==P1[1])&&(CornerPos.count(P0[0])==0))
                return true;

            if ((P0[1]==P1[0])&&(CornerPos.count(P0[1])==0))
                return true;

            if ((P0[1]==P1[1])&&(CornerPos.count(P0[1])==0))
                return true;
        }
    }
    return false;
}

template<class TriangleMeshType>
void MergeProjBasis(const TriangleMeshType &tri_mesh,
                    const std::vector<size_t> &featuresC,
                    std::vector<int> &VertBasis,
                    std::vector<std::vector<vcg::Segment3<typename TriangleMeshType::ScalarType> > > &ProjBasis)
{
    typedef typename TriangleMeshType::ScalarType ScalarType;
    typedef typename TriangleMeshType::CoordType CoordType;

    //set the coordinates of the corners
    std::set<CoordType> CornerPos;
    for (size_t i=0;i<featuresC.size();i++)
        CornerPos.insert(tri_mesh.vert[featuresC[i]].P());

    bool merged=false;
    do
    {
        merged=false;
        if (ProjBasis.size()==1)return;
        int MergeI0=-1;
        int MergeI1=-1;
        for (size_t i=0;i<ProjBasis.size()-1;i++)
        {
            if (ProjBasis[i].size()==0)continue;

            for (size_t j=(i+1);j<ProjBasis.size();j++)
            {
                if (ProjBasis[j].size()==0)continue;
                if (IsMergable<TriangleMeshType>(ProjBasis[i],ProjBasis[j],CornerPos))
                {
                    MergeI0=i;
                    MergeI1=j;
                    merged=true;
                    break;
                }
            }
        }

        if ((MergeI0>=0)&&(MergeI1>=0))
        {
            assert(MergeI0!=MergeI1);
            Merge(MergeI0,MergeI1,VertBasis,ProjBasis);
        }

    }while(merged);
}

template <class PolyMeshType,class TriangleMeshType>
void GetVertProjBasis(const typename PolyMeshType::ScalarType &AvEdge,
                      const std::vector<std::pair<size_t,size_t> > &FeatureEdges,
                      const std::vector<size_t> &featuresC,
                      TriangleMeshType &tri_mesh,
                      const std::vector<size_t > &quad_corner,
                      PolyMeshType &poly_mesh,
                      const std::vector<size_t> &tri_face_partition,
                      const std::vector<size_t> &quad_face_partition,
                      std::vector<int> &VertBasis,
                      std::vector<std::vector<vcg::Segment3<typename PolyMeshType::ScalarType> > > &ProjBasis)
{
    typedef typename BasicMesh::ScalarType ScalarType;
    typedef typename BasicMesh::CoordType CoordType;
    typedef typename vcg::Segment3<typename PolyMeshType::ScalarType> SegmentType;

    VertBasis.clear();
    ProjBasis.clear();

    //std::cout<<"A"<<std::endl;

    std::vector<size_t> cornersIdx;
    GetCorners(FeatureEdges,featuresC,tri_mesh,quad_corner,poly_mesh,AvEdge,cornersIdx);

    //GetCorners(featuresC,tri_mesh,poly_mesh,tri_face_partition,quad_face_partition,cornersIdx);

    //find each segment set per pair of partitions
    typedef std::pair<int,int> PatchPairKey;
    std::map<PatchPairKey,size_t > BasisMap;

    for (size_t i=0;i<FeatureEdges.size();i++)
    {
        size_t IndexF=FeatureEdges[i].first;
        size_t IndexE=FeatureEdges[i].second;
        int Part0=tri_face_partition[IndexF];
        int Part1=-1;//when is adjacent to border
        if (!vcg::face::IsBorder(tri_mesh.face[IndexF],IndexE))
        {
            size_t IndexFOpp=vcg::tri::Index(tri_mesh,tri_mesh.face[IndexF].cFFp(IndexE));
            Part1=tri_face_partition[IndexFOpp];
        }
        assert(Part0!=Part1);
        PatchPairKey key(std::min(Part0,Part1),std::max(Part0,Part1));

        CoordType Pos0=tri_mesh.face[IndexF].cP0(IndexE);
        CoordType Pos1=tri_mesh.face[IndexF].cP1(IndexE);
        SegmentType SBasis(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
        if (BasisMap.count(key)==0)
        {
            ProjBasis.resize(ProjBasis.size()+1);
            ProjBasis.back().push_back(SBasis);
            BasisMap[key]=ProjBasis.size()-1;
        }
        else
        {
            int IndexP=BasisMap[key];
            ProjBasis[IndexP].push_back(SBasis);
        }
    }

    //std::cout<<"B"<<std::endl;

    //then check the edge of the quads, by default =-2, no projection, internal
    VertBasis.resize(poly_mesh.vert.size(),-2);
    //    BasicMesh EdgeMesh;
    //    ExtractEdgeMesh(tri_mesh,features,EdgeMesh);

    for (size_t i=0;i<poly_mesh.face.size();i++)
    {
        int Part0=quad_face_partition[i];

        //get neighbours
        for (size_t j=0;j<poly_mesh.face[i].VN();j++)
        {
            //default border
            int Part1=-1;
            if (!vcg::face::IsBorder(poly_mesh.face[i],j))
            {
                size_t IndexFOpp=vcg::tri::Index(poly_mesh,poly_mesh.face[i].cFFp(j));
                Part1=quad_face_partition[IndexFOpp];
            }
            std::pair<size_t,size_t> key(std::min(Part0,Part1),std::max(Part0,Part1));
            if (BasisMap.count(key)>0)
            {
                int IndV0=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V0(j));
                int IndV1=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V1(j));

                //                ClosestPointEMesh(poly_mesh.vert[IndV0].cP(),edge_mesh,IndexE,t,MinD,Clos);
                //                if (MinD>AvEdge/10)continue;

                VertBasis[IndV0]=BasisMap[key];
                VertBasis[IndV1]=BasisMap[key];
            }
        }
    }

    //std::cout<<"C"<<std::endl;
    for (size_t i=0;i<cornersIdx.size();i++)
        VertBasis[cornersIdx[i]]=-1;

    if (ProjBasis.size()>0)
    MergeProjBasis(tri_mesh,featuresC,VertBasis,ProjBasis);
    //std::cout<<"D"<<std::endl;
}

enum SmoothType{Laplacian,TemplateFit};


//template <class PolyMeshType,class TriangleMeshType>
//void SmoothWithFeatures(TriangleMeshType &tri_mesh,
//                        PolyMeshType &poly_mesh,
//                        const std::vector<std::pair<size_t,size_t> > &features,
//                        const std::vector<size_t> &tri_face_partition,
//                        const std::vector<size_t> &quad_face_partition,
//                        SmoothType SType,
//                        size_t Steps,
//                        typename PolyMeshType::ScalarType Damp,
//                        typename PolyMeshType::ScalarType AvEdge)
//{
//    typedef typename TriangleMeshType::FaceType TriFaceType;
//    typedef typename TriangleMeshType::ScalarType TriScalarType;
//    typedef typename TriangleMeshType::CoordType CoordType;

//    typedef typename PolyMeshType::VertexType PolyVertexType;

//    BasicMesh EdgeMesh;
//    ExtractEdgeMesh(tri_mesh,features,EdgeMesh);

//    std::vector<std::vector<bool> > IsSharp;
//    std::vector<VType> VertTypes;
//    GetVertType(EdgeMesh,AvEdge,tri_mesh,poly_mesh,features,tri_face_partition,quad_face_partition,IsSharp,VertTypes);

//    vcg::GridStaticPtr<TriFaceType,TriScalarType> TriGrid;
//    TriGrid.Set(tri_mesh.face.begin(),tri_mesh.face.end());

//    //    vcg::tri::io::ExporterOBJ<BasicMesh>::Save(EdgeMesh,"test_edge.obj",vcg::tri::io::Mask::IOM_EDGEINDEX);
//    //    return;
//    for (size_t s=0;s<Steps;s++)
//    {
//        //Smooth only along features
//        std::vector<CoordType> TargetPos(poly_mesh.vert.size(),CoordType(0,0,0));
//        std::vector<size_t> TargetNum(poly_mesh.vert.size(),0);
//        vcg::tri::UpdateSelection<PolyMeshType>::VertexClear(poly_mesh);
//        for (size_t i=0;i<poly_mesh.face.size();i++)
//        {
//            size_t NumV=poly_mesh.face[i].VN();
//            for (size_t j=0;j<poly_mesh.face[i].VN();j++)
//            {
//                if (!IsSharp[i][j])continue;

//                PolyVertexType *V0=poly_mesh.face[i].V(j);
//                PolyVertexType *V1=poly_mesh.face[i].V((j+1)%NumV);

//                size_t IndexV0=vcg::tri::Index(poly_mesh,V0);
//                size_t IndexV1=vcg::tri::Index(poly_mesh,V1);

//                poly_mesh.vert[IndexV0].SetS();
//                poly_mesh.vert[IndexV1].SetS();

//                CoordType PosV0=poly_mesh.vert[IndexV0].P();
//                CoordType PosV1=poly_mesh.vert[IndexV1].P();
//                TargetPos[IndexV0]+=PosV1;
//                TargetPos[IndexV1]+=PosV0;
//                TargetNum[IndexV0]++;
//                TargetNum[IndexV1]++;
//            }
//        }
//        for (size_t i=0;i<poly_mesh.vert.size();i++)
//        {
//            if (VertTypes[i]!=Feature)continue;
//            if (TargetNum[i]==0)continue;

//            //then average
//            CoordType TargetP=TargetPos[i]/TargetNum[i];
//            poly_mesh.vert[i].P()=poly_mesh.vert[i].P()*Damp+TargetP*(1-Damp);

//            //and reproject
//            size_t IndexE;
//            TriScalarType t,MinD;
//            CoordType closestPt;
//            ClosestPointEMesh(poly_mesh.vert[i].P(),EdgeMesh,IndexE,t,MinD,closestPt);
//            poly_mesh.vert[i].P()=closestPt;
//        }

//        if (SType==Laplacian)
//            vcg::PolygonalAlgorithm<PolyMeshType>::Laplacian(poly_mesh,true,1,Damp);
//        else
//            vcg::PolygonalAlgorithm<PolyMeshType>::SmoothPCA(poly_mesh,1,Damp,true,true,0.5,false);

//        //then reproject
//        for (size_t i=0;i<poly_mesh.vert.size();i++)
//        {
//            if (VertTypes[i]!=Internal)continue;

//            TriScalarType MaxD=poly_mesh.bbox.Diag();
//            TriScalarType MinD;
//            CoordType closestPt;
//            vcg::tri::GetClosestFaceBase(tri_mesh,TriGrid,poly_mesh.vert[i].P(),MaxD,MinD,closestPt);
//            poly_mesh.vert[i].P()=closestPt;
//        }
//    }
//    //    vcg::tri::io::ExporterOBJ<BasicMesh>::Save(EdgeMesh,"test_edge.obj",vcg::tri::io::Mask::IOM_EDGEINDEX);
//    //return;
//}

//template <class PolyMeshType,class TriangleMeshType>
//void SmoothWithFeatures(TriangleMeshType &tri_mesh,
//                        PolyMeshType &poly_mesh,
//                        const std::vector<std::pair<size_t,size_t> > &features,
//                        const std::vector<size_t> &tri_face_partition,
//                        const std::vector<size_t> &quad_face_partition,
//                        SmoothType SType,
//                        size_t Steps,
//                        typename PolyMeshType::ScalarType Damp,
//                        typename PolyMeshType::ScalarType AvEdge)
//{
//    typedef typename TriangleMeshType::FaceType TriFaceType;
//    typedef typename TriangleMeshType::ScalarType TriScalarType;
//    typedef typename TriangleMeshType::CoordType CoordType;

//    typedef typename PolyMeshType::VertexType PolyVertexType;

//    BasicMesh EdgeMesh;
//    ExtractEdgeMesh(tri_mesh,features,EdgeMesh);

//    std::vector<std::vector<bool> > IsSharp;
//    std::vector<VType> VertTypes;
//    GetVertType(EdgeMesh,AvEdge,tri_mesh,poly_mesh,features,tri_face_partition,quad_face_partition,IsSharp,VertTypes);

//    //select corners
////    vcg::tri::UpdateSelection<PolyMeshType>::VertexClear(poly_mesh);
////    for (size_t i=0;i<VertTypes.size();i++)
////        if (VertTypes[i]==Corner)
////            poly_mesh.vert[i].SetS();

//    vcg::GridStaticPtr<TriFaceType,TriScalarType> TriGrid;
//    TriGrid.Set(tri_mesh.face.begin(),tri_mesh.face.end());

//    //    vcg::tri::io::ExporterOBJ<BasicMesh>::Save(EdgeMesh,"test_edge.obj",vcg::tri::io::Mask::IOM_EDGEINDEX);
//    //    return;
//    for (size_t s=0;s<Steps;s++)
//    {

//        if (SType==Laplacian)
//            vcg::PolygonalAlgorithm<PolyMeshType>::Laplacian(poly_mesh,true,1,Damp);
//        else
//            vcg::PolygonalAlgorithm<PolyMeshType>::SmoothPCA(poly_mesh,1,Damp,true,true,0.5,false);

//        //then reproject
//        for (size_t i=0;i<poly_mesh.vert.size();i++)
//        {
//            //if (VertTypes[i]==Corner)continue;

//            size_t IndexE;
//            TriScalarType t,MinD;
//            CoordType closestPt;
//            if (VertTypes[i]!=Internal)
//            {
//                ClosestPointEMesh(poly_mesh.vert[i].P(),EdgeMesh,IndexE,t,MinD,closestPt);
//                poly_mesh.vert[i].P()=closestPt;
//            }
//            TriScalarType MaxD=poly_mesh.bbox.Diag();
////            TriScalarType MinD;
////            CoordType closestPt;
//            vcg::tri::GetClosestFaceBase(tri_mesh,TriGrid,poly_mesh.vert[i].P(),MaxD,MinD,closestPt);
//            poly_mesh.vert[i].P()=closestPt;
//        }
//    }
//    //    vcg::tri::io::ExporterOBJ<BasicMesh>::Save(EdgeMesh,"test_edge.obj",vcg::tri::io::Mask::IOM_EDGEINDEX);
//    //return;
//}

template <class PolyMeshType>
typename PolyMeshType::CoordType ClosestPointSegSet(const typename BasicMesh::CoordType &Pos,
                                                    std::vector<vcg::Segment3<typename BasicMesh::ScalarType> > &SegSet)
{
    typedef typename PolyMeshType::ScalarType ScalarType;
    typedef typename PolyMeshType::CoordType CoordType;
    ScalarType MinD=std::numeric_limits<ScalarType>::max();
    CoordType Clos;
    for (size_t i=0;i<SegSet.size();i++)
    {
        vcg::Segment3<ScalarType> STest=SegSet[i];
        ScalarType testD;
        CoordType ClosTest;
        vcg::SegmentPointDistance(STest,Pos,ClosTest,testD);
        if (testD>MinD)continue;
        Clos=ClosTest;
        MinD=testD;
    }
    return Clos;
}

template <class PolyMeshType>
void LaplacianEdge(PolyMeshType &poly_mesh,bool FixS,typename PolyMeshType::ScalarType Damp)
{
    typedef typename PolyMeshType::ScalarType ScalarType;
    typedef typename PolyMeshType::CoordType CoordType;
    std::vector<CoordType> SumPos(poly_mesh.vert.size(),CoordType(0,0,0));
    std::vector<size_t> NumPos(poly_mesh.vert.size(),0);
    for (size_t i=0;i<poly_mesh.face.size();i++)
    {
        int sizeP=poly_mesh.face[i].VN();
        for (size_t j=0;j<poly_mesh.face[i].VN();j++)
        {
            CoordType Pos0=poly_mesh.face[i].P(j);
            CoordType Pos1=poly_mesh.face[i].P((j+1)%sizeP);
            size_t VIndex0=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V(j));
            size_t VIndex1=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V((j+1)%sizeP));
            SumPos[VIndex0]+=Pos1;
            SumPos[VIndex1]+=Pos0;
            NumPos[VIndex0]++;
            NumPos[VIndex1]++;
        }
    }

    for (size_t i=0;i<SumPos.size();i++)
    {
        if (NumPos[i]==0)continue;
        SumPos[i]/=NumPos[i];
    }

    for (size_t i=0;i<poly_mesh.vert.size();i++)
    {
        if (FixS && (poly_mesh.vert[i].IsS()))continue;
        poly_mesh.vert[i].P()=poly_mesh.vert[i].P()*Damp+SumPos[i]*(1-Damp);
    }

}

template <class PolyMeshType,class TriangleMeshType>
void SmoothWithFeatures(TriangleMeshType &tri_mesh,
                        PolyMeshType &poly_mesh,
                        const std::vector<std::pair<size_t,size_t> > &features,
                        const std::vector<size_t> &featuresC,
                        const std::vector<size_t> &tri_face_partition,
                        const std::vector<size_t > &quad_corner,
                        const std::vector<size_t> &quad_face_partition,
                        SmoothType SType,
                        size_t Steps,
                        typename PolyMeshType::ScalarType Damp,
                        typename PolyMeshType::ScalarType AvEdge)
{
    typedef typename TriangleMeshType::FaceType TriFaceType;
    typedef typename TriangleMeshType::ScalarType TriScalarType;
    typedef typename TriangleMeshType::CoordType CoordType;

    typedef typename PolyMeshType::VertexType PolyVertexType;

    std::cout<<"*** Getting Projection basis ***"<<std::endl;
    std::vector<int> VertBasis;
    std::vector<std::vector<vcg::Segment3<typename PolyMeshType::ScalarType> > > ProjBasis;
    GetVertProjBasis(AvEdge,features,featuresC,tri_mesh,quad_corner,poly_mesh,tri_face_partition,quad_face_partition,VertBasis,ProjBasis);
    std::cout<<"*** End Getting Projection basis ***"<<std::endl;
    //GetVertProjBasis(AvEdge,features,featuresC,tri_mesh,poly_mesh,tri_face_partition,quad_face_partition,VertBasis,ProjBasis);


    BasicMesh EdgeMesh;
    ExtractEdgeMesh(tri_mesh,features,EdgeMesh);

    vcg::tri::io::ExporterOBJ<BasicMesh>::Save(EdgeMesh,"test_edge.obj",vcg::tri::io::Mask::IOM_EDGEINDEX);

    //select corners
    vcg::tri::UpdateSelection<PolyMeshType>::VertexClear(poly_mesh);
    for (size_t i=0;i<VertBasis.size();i++)
        if (VertBasis[i]==-1)
            poly_mesh.vert[i].SetS();

    vcg::GridStaticPtr<TriFaceType,TriScalarType> TriGrid;
    TriGrid.Set(tri_mesh.face.begin(),tri_mesh.face.end());

    for (size_t s=0;s<Steps;s++)
    {
        if (SType==Laplacian)
            LaplacianEdge(poly_mesh,true,Damp);
        //vcg::PolygonalAlgorithm<PolyMeshType>::Laplacian(poly_mesh,true,1,Damp);
        else
            vcg::PolygonalAlgorithm<PolyMeshType>::SmoothPCA(poly_mesh,1,Damp,true,true,0.2,false);
        //then reproject
        for (size_t i=0;i<poly_mesh.vert.size();i++)
        {
            //if (VertTypes[i]==Corner)continue;


            CoordType closestPt;
            if (VertBasis[i]>=0)
            {
                closestPt=ClosestPointSegSet<PolyMeshType>(poly_mesh.vert[i].P(),ProjBasis[VertBasis[i]]);
                poly_mesh.vert[i].P()=closestPt;
            }
            else
            {
                TriScalarType MinD;
                TriScalarType MaxD=poly_mesh.bbox.Diag();
                vcg::tri::GetClosestFaceBase(tri_mesh,TriGrid,poly_mesh.vert[i].P(),MaxD,MinD,closestPt);
                poly_mesh.vert[i].P()=closestPt;
            }
        }
    }

}
#endif // SMOOTH_MESH_H
