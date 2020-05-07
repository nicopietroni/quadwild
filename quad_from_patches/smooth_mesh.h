#ifndef SMOOTH_MESH_H
#define SMOOTH_MESH_H

#include <vector>
#include <vcg/complex/algorithms/polygonal_algorithms.h>
#include <wrap/io_trimesh/export.h>
#include <wrap/igl/smooth_field.h>
#include <vcg/complex/algorithms/implicit_smooth.h>

///* ----- Triangle mesh ----- */

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

//void SelectCorners(BasicMesh &EdgeMesh)
//{
//    std::vector<size_t> NumE(EdgeMesh.vert.size(),0);
//    vcg::tri::UpdateFlags<BasicMesh>::VertexClearS(EdgeMesh);
//    for (size_t i=0;i<EdgeMesh.edge.size();i++)
//    {
//        size_t IndexV0=vcg::tri::Index(EdgeMesh,EdgeMesh.edge[i].V(0));
//        size_t IndexV1=vcg::tri::Index(EdgeMesh,EdgeMesh.edge[i].V(1));
//        NumE[IndexV0]++;
//        NumE[IndexV1]++;
//    }
//    for (size_t i=0;i<NumE.size();i++)
//        if (NumE[i]!=2)EdgeMesh.vert[i].SetS();
//}

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

//enum VType{Internal,Feature,Corner};



//template <class PolyMeshType,class TriangleMeshType>
//void GetVertType(const BasicMesh &edge_mesh,
//                 const typename PolyMeshType::ScalarType &AvEdge,
//                 const TriangleMeshType &tri_mesh,
//                 const PolyMeshType &poly_mesh,
//                 const std::vector<std::pair<size_t,size_t> > FeatureEdges,
//                 const std::vector<size_t> &tri_face_partition,
//                 const std::vector<size_t> &quad_face_partition,
//                 std::vector<std::vector<bool> > &IsSharp,
//                 std::vector<VType> &VertTypes)
//{
//    typedef typename BasicMesh::ScalarType ScalarType;
//    typedef typename BasicMesh::CoordType CoordType;


//    std::set<std::pair<size_t,size_t> > SharpPart;

//    for (size_t i=0;i<FeatureEdges.size();i++)
//    {
//        size_t IndexF=FeatureEdges[i].first;
//        size_t IndexE=FeatureEdges[i].second;
//        if (vcg::face::IsBorder(tri_mesh.face[IndexF],IndexE))continue;
//        size_t Part0=tri_face_partition[IndexF];

//        size_t IndexFOpp=vcg::tri::Index(tri_mesh,tri_mesh.face[IndexF].cFFp(IndexE));
//        size_t Part1=tri_face_partition[IndexFOpp];
//        assert(Part0!=Part1);
//        std::pair<size_t,size_t> key(std::min(Part0,Part1),std::max(Part0,Part1));
//        SharpPart.insert(key);
//    }

//    //then check the edge of the quads
//    IsSharp.resize(poly_mesh.face.size());
//    for (size_t i=0;i<poly_mesh.face.size();i++)
//    {
//        size_t Part0=quad_face_partition[i];
//        IsSharp[i].resize(poly_mesh.face[i].VN(),false);
//        for (size_t j=0;j<poly_mesh.face[i].VN();j++)
//        {
//            if (vcg::face::IsBorder(poly_mesh.face[i],j))
//                IsSharp[i][j]=true;
//            else
//            {
//                size_t IndexFOpp=vcg::tri::Index(poly_mesh,poly_mesh.face[i].cFFp(j));
//                size_t Part1=quad_face_partition[IndexFOpp];
//                std::pair<size_t,size_t> key(std::min(Part0,Part1),std::max(Part0,Part1));

//                if (SharpPart.count(key))
//                    IsSharp[i][j]=true;
//            }
//        }
//    }

//    //finally set the vertex type
//    VertTypes.clear();
//    VertTypes.resize(poly_mesh.vert.size(),Internal);
//    std::vector<size_t> SharpVal(poly_mesh.vert.size(),0);
//    for (size_t i=0;i<poly_mesh.face.size();i++)
//    {
//        size_t NumV=poly_mesh.face[i].VN();
//        for (size_t j=0;j<NumV;j++)
//        {
//            if (!IsSharp[i][j])continue;
//            size_t IndexV0=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V(j));
//            size_t IndexV1=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V((j+1)%NumV));
//            SharpVal[IndexV0]++;
//            SharpVal[IndexV1]++;
//        }
//    }

//    for (size_t i=0;i<SharpVal.size();i++)
//    {
//        if (SharpVal[i]==0)continue;

//        size_t IndexE;
//        BasicMesh::ScalarType t,MinD;
//        BasicMesh::CoordType Clos;

//        ClosestPointEMesh(poly_mesh.vert[i].cP(),edge_mesh,IndexE,t,MinD,Clos);

//        if (MinD>AvEdge/10)continue;

//        if (SharpVal[i]==4)
//            VertTypes[i]=Feature;
//        else
//            VertTypes[i]=Corner;
//    }
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
    typedef typename PolyMeshType::ScalarType ScalarType;
    typedef typename PolyMeshType::CoordType CoordType;

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

//template<class ScalarType>
//void Merge(const int &IndexMerge0,
//           const int &IndexMerge1,
//           std::vector<int> &VertBasis,
//           std::vector<std::vector<vcg::Segment3<ScalarType> > > &ProjBasis)
//{
//    assert(IndexMerge0>=0);
//    assert(IndexMerge1>=0);
//    assert(IndexMerge0<ProjBasis.size());
//    assert(IndexMerge1<ProjBasis.size());
//    ProjBasis[IndexMerge0].insert(ProjBasis[IndexMerge0].end(),
//                                  ProjBasis[IndexMerge1].begin(),
//                                  ProjBasis[IndexMerge1].end());

//    ProjBasis[IndexMerge1].clear();
//    for (size_t i=0;i<VertBasis.size();i++)
//        if (VertBasis[i]==IndexMerge1)
//            VertBasis[i]=IndexMerge0;
//}

//template<class TriangleMeshType>
//bool IsMergable(const std::vector<vcg::Segment3<typename TriangleMeshType::ScalarType> > &SegBase0,
//                const std::vector<vcg::Segment3<typename TriangleMeshType::ScalarType> > &SegBase1,
//                const std::set<typename TriangleMeshType::CoordType> &CornerPos)
//{
//    typedef typename TriangleMeshType::ScalarType ScalarType;
//    typedef typename TriangleMeshType::CoordType CoordType;

//    for (size_t i=0;i<SegBase0.size();i++)
//    {
//        CoordType P0[2];
//        P0[0]=SegBase0[i].P0();
//        P0[1]=SegBase0[i].P1();
//        for (size_t j=0;j<SegBase1.size();j++)
//        {
//            CoordType P1[2];
//            P1[0]=SegBase1[j].P0();
//            P1[1]=SegBase1[j].P1();
//            if ((P0[0]==P1[0])&&(CornerPos.count(P0[0])==0))
//                return true;

//            if ((P0[0]==P1[1])&&(CornerPos.count(P0[0])==0))
//                return true;

//            if ((P0[1]==P1[0])&&(CornerPos.count(P0[1])==0))
//                return true;

//            if ((P0[1]==P1[1])&&(CornerPos.count(P0[1])==0))
//                return true;
//        }
//    }
//    return false;
//}

//template<class TriangleMeshType>
//void MergeProjBasis(const TriangleMeshType &tri_mesh,
//                    const std::vector<size_t> &featuresC,
//                    std::vector<int> &VertBasis,
//                    std::vector<std::vector<vcg::Segment3<typename TriangleMeshType::ScalarType> > > &ProjBasis)
//{
//    typedef typename TriangleMeshType::ScalarType ScalarType;
//    typedef typename TriangleMeshType::CoordType CoordType;

//    //set the coordinates of the corners
//    std::set<CoordType> CornerPos;
//    for (size_t i=0;i<featuresC.size();i++)
//        CornerPos.insert(tri_mesh.vert[featuresC[i]].P());

//    bool merged=false;
//    do
//    {
//        merged=false;
//        if (ProjBasis.size()==1)return;
//        int MergeI0=-1;
//        int MergeI1=-1;
//        for (size_t i=0;i<ProjBasis.size()-1;i++)
//        {
//            if (ProjBasis[i].size()==0)continue;

//            for (size_t j=(i+1);j<ProjBasis.size();j++)
//            {
//                if (ProjBasis[j].size()==0)continue;
//                if (IsMergable<TriangleMeshType>(ProjBasis[i],ProjBasis[j],CornerPos))
//                {
//                    MergeI0=i;
//                    MergeI1=j;
//                    merged=true;
//                    break;
//                }
//            }
//        }

//        if ((MergeI0>=0)&&(MergeI1>=0))
//        {
//            assert(MergeI0!=MergeI1);
//            Merge(MergeI0,MergeI1,VertBasis,ProjBasis);
//        }

//    }while(merged);
//}


//template <class PolyMeshType,class TriangleMeshType>
//void GetVertProjBasis(const typename PolyMeshType::ScalarType &AvEdge,
//                      const std::vector<std::pair<size_t,size_t> > &FeatureEdges,
//                      const std::vector<size_t> &featuresC,
//                      TriangleMeshType &tri_mesh,
//                      const std::vector<size_t > &quad_corner,
//                      PolyMeshType &poly_mesh,
//                      const std::vector<size_t> &tri_face_partition,
//                      const std::vector<size_t> &quad_face_partition,
//                      std::vector<int> &VertBasis,
//                      std::vector<std::vector<vcg::Segment3<typename PolyMeshType::ScalarType> > > &ProjBasis)
//{
//    typedef typename PolyMeshType::ScalarType ScalarType;
//    typedef typename PolyMeshType::CoordType CoordType;
//    typedef typename vcg::Segment3<typename PolyMeshType::ScalarType> SegmentType;

//    VertBasis.clear();
//    ProjBasis.clear();

//    //std::cout<<"A"<<std::endl;

//    std::vector<size_t> cornersIdx;
//    GetCorners(FeatureEdges,featuresC,tri_mesh,quad_corner,poly_mesh,AvEdge,cornersIdx);

//    //GetCorners(featuresC,tri_mesh,poly_mesh,tri_face_partition,quad_face_partition,cornersIdx);

//    //find each segment set per pair of partitions
//    typedef std::pair<int,int> PatchPairKey;
//    std::map<PatchPairKey,size_t > BasisMap;

//    for (size_t i=0;i<FeatureEdges.size();i++)
//    {
//        size_t IndexF=FeatureEdges[i].first;
//        size_t IndexE=FeatureEdges[i].second;
//        int Part0=tri_face_partition[IndexF];
//        int Part1=-1;//when is adjacent to border
//        if (!vcg::face::IsBorder(tri_mesh.face[IndexF],IndexE))
//        {
//            size_t IndexFOpp=vcg::tri::Index(tri_mesh,tri_mesh.face[IndexF].cFFp(IndexE));
//            Part1=tri_face_partition[IndexFOpp];
//        }
//        assert(Part0!=Part1);
//        PatchPairKey key(std::min(Part0,Part1),std::max(Part0,Part1));

//        CoordType Pos0=tri_mesh.face[IndexF].cP0(IndexE);
//        CoordType Pos1=tri_mesh.face[IndexF].cP1(IndexE);
//        SegmentType SBasis(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
//        if (BasisMap.count(key)==0)
//        {
//            ProjBasis.resize(ProjBasis.size()+1);
//            ProjBasis.back().push_back(SBasis);
//            BasisMap[key]=ProjBasis.size()-1;
//        }
//        else
//        {
//            int IndexP=BasisMap[key];
//            ProjBasis[IndexP].push_back(SBasis);
//        }
//    }

//    //std::cout<<"B"<<std::endl;

//    //then check the edge of the quads, by default =-2, no projection, internal
//    VertBasis.resize(poly_mesh.vert.size(),-2);
//    //    BasicMesh EdgeMesh;
//    //    ExtractEdgeMesh(tri_mesh,features,EdgeMesh);

//    for (size_t i=0;i<poly_mesh.face.size();i++)
//    {
//        int Part0=quad_face_partition[i];

//        //get neighbours
//        for (size_t j=0;j<poly_mesh.face[i].VN();j++)
//        {
//            //default border
//            int Part1=-1;
//            if (!vcg::face::IsBorder(poly_mesh.face[i],j))
//            {
//                size_t IndexFOpp=vcg::tri::Index(poly_mesh,poly_mesh.face[i].cFFp(j));
//                Part1=quad_face_partition[IndexFOpp];
//            }
//            std::pair<size_t,size_t> key(std::min(Part0,Part1),std::max(Part0,Part1));
//            if (BasisMap.count(key)>0)
//            {
//                int IndV0=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V0(j));
//                int IndV1=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V1(j));

//                //                ClosestPointEMesh(poly_mesh.vert[IndV0].cP(),edge_mesh,IndexE,t,MinD,Clos);
//                //                if (MinD>AvEdge/10)continue;

//                VertBasis[IndV0]=BasisMap[key];
//                VertBasis[IndV1]=BasisMap[key];
//            }
//        }
//    }

//    //std::cout<<"C"<<std::endl;
//    for (size_t i=0;i<cornersIdx.size();i++)
//        VertBasis[cornersIdx[i]]=-1;

//    if (ProjBasis.size()>0)
//        MergeProjBasis(tri_mesh,featuresC,VertBasis,ProjBasis);
//    //std::cout<<"D"<<std::endl;
//}

enum ProjType{ProjNone,ProjSuface,ProjSharp,ProjCorner};
struct ProjectionBase
{
    std::vector<ProjType> VertProjType;
    //edge sequences as pair of vertices
    std::vector<std::vector<std::pair<size_t,size_t> > > SharpEdge;
    //indexes of vertices
    std::vector<size_t> Corner;
    void Clear()
    {VertProjType.clear();SharpEdge.clear();Corner.clear();}
};

enum SmoothType{Laplacian,TemplateFit};

template <class PolyMeshType>
void LaplacianEdgePos(const PolyMeshType &poly_mesh,
                      const ProjectionBase &PolyProjBase,
                      const typename PolyMeshType::ScalarType Damp,
                      std::vector<typename PolyMeshType::CoordType> &TargetPos)
{
    typedef typename PolyMeshType::ScalarType ScalarType;
    typedef typename PolyMeshType::CoordType CoordType;

    TargetPos = std::vector<CoordType>(poly_mesh.vert.size(),CoordType(0,0,0));
    std::vector<size_t> NumPos(poly_mesh.vert.size(),0);
    for (size_t i=0;i<poly_mesh.face.size();i++)
    {
        int sizeP=poly_mesh.face[i].VN();
        for (size_t j=0;j<poly_mesh.face[i].VN();j++)
        {
            CoordType Pos0=poly_mesh.face[i].cP(j);
            CoordType Pos1=poly_mesh.face[i].cP((j+1)%sizeP);
            size_t VIndex0=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V(j));
            size_t VIndex1=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V((j+1)%sizeP));
            TargetPos[VIndex0]+=Pos1;
            TargetPos[VIndex1]+=Pos0;
            NumPos[VIndex0]++;
            NumPos[VIndex1]++;
        }
    }

    for (size_t i=0;i<TargetPos.size();i++)
    {
        if (NumPos[i]==0)
            TargetPos[i]=poly_mesh.vert[i].cP();
        else
            TargetPos[i]/=NumPos[i];
    }

    for (size_t i=0;i<poly_mesh.vert.size();i++)
    {
        //if (FixS && (poly_mesh.vert[i].IsS()))continue;
        if (PolyProjBase.VertProjType[i]==ProjCorner)
            TargetPos[i]=poly_mesh.vert[i].cP();
        else
            TargetPos[i]=poly_mesh.vert[i].cP()*Damp+TargetPos[i]*(1-Damp);
    }
}

template <class PolyMeshType>
void TemplatePolyPos(PolyMeshType &poly_mesh,
                     const ProjectionBase &PolyProjBase,
                     const typename PolyMeshType::ScalarType Damp,
                     std::vector<typename PolyMeshType::CoordType> &TargetPos)
{
    typedef typename PolyMeshType::ScalarType ScalarType;
    typedef typename PolyMeshType::CoordType CoordType;

    //save old position
    vcg::tri::UpdateFlags<PolyMeshType>::VertexClearS(poly_mesh);
    std::vector<CoordType> OldPos;
    for (size_t i=0;i<poly_mesh.vert.size();i++)
    {
        if (PolyProjBase.VertProjType[i]==ProjCorner)
            poly_mesh.vert[i].SetS();
        OldPos.push_back(poly_mesh.vert[i].P());
    }

    //smooth
    vcg::PolygonalAlgorithm<PolyMeshType>::SmoothPCA(poly_mesh,10,Damp,true);
    TargetPos.clear();
    for (size_t i=0;i<poly_mesh.vert.size();i++)
    {
        TargetPos.push_back(poly_mesh.vert[i].P());
        poly_mesh.vert[i].P()=OldPos[i];
    }
    vcg::tri::UpdateFlags<PolyMeshType>::VertexClearS(poly_mesh);
    //    TargetPos = std::vector<CoordType>(poly_mesh.vert.size(),CoordType(0,0,0));
    //    std::vector<size_t> NumPos(poly_mesh.vert.size(),0);
    //    for (size_t i=0;i<poly_mesh.face.size();i++)
    //    {
    //        int sizeP=poly_mesh.face[i].VN();
    //        for (size_t j=0;j<poly_mesh.face[i].VN();j++)
    //        {
    //            CoordType Pos0=poly_mesh.face[i].cP(j);
    //            CoordType Pos1=poly_mesh.face[i].cP((j+1)%sizeP);
    //            size_t VIndex0=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V(j));
    //            size_t VIndex1=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V((j+1)%sizeP));
    //            TargetPos[VIndex0]+=Pos1;
    //            TargetPos[VIndex1]+=Pos0;
    //            NumPos[VIndex0]++;
    //            NumPos[VIndex1]++;
    //        }
    //    }

    //    for (size_t i=0;i<TargetPos.size();i++)
    //    {
    //        if (NumPos[i]==0)
    //            TargetPos[i]=poly_mesh.vert[i].cP();
    //        else
    //            TargetPos[i]/=NumPos[i];
    //    }

    //    for (size_t i=0;i<poly_mesh.vert.size();i++)
    //    {
    //        //if (FixS && (poly_mesh.vert[i].IsS()))continue;
    //        if (PolyProjBase.VertProjType[i]==ProjCorner)
    //            TargetPos[i]=poly_mesh.vert[i].cP();
    //        else
    //            TargetPos[i]=poly_mesh.vert[i].cP()*Damp+TargetPos[i]*(1-Damp);
    //    }
}

template <class PolyMeshType,class TriMeshType>
void GetQuadInterpW(const typename TriMeshType::FaceType &FTris,
                    const typename PolyMeshType::FaceType &FPoly,
                    const typename TriMeshType::CoordType &ClosestPt,
                    std::vector<typename TriMeshType::ScalarType> &VertWeigths)
{
    typedef typename PolyMeshType::ScalarType ScalarType;
    typedef typename PolyMeshType::CoordType CoordType;
    typedef typename TriMeshType::FaceType TriFaceType;

    VertWeigths.clear();
    VertWeigths=std::vector<ScalarType>(4,0);

    //get the barycentric coordinates
    CoordType bary;
    //bool Interpolated=
    vcg::InterpolationParameters(FTris,ClosestPt,bary);

    //then check which are the index of the vertices
    int IndexV[3]={-1,-1,-1};
    for (size_t i=0;i<3;i++)
    {
        for (size_t j=0;j<4;j++)
        {
            if (FTris.cV(i)->cQ()==FPoly.cV(j)->cQ())
            {
                IndexV[i]=j;
                break;
            }
        }
    }

    //then check there in only and only one -1
    assert((IndexV[0]==-1)||(IndexV[1]==-1)||(IndexV[2]==-1));

    //   //and all different
    //   std::cout<<" "<<FTris.cV(0)->cQ()<<" "<<FTris.cV(1)->Q()<<" "<<FTris.cV(2)->cQ()<<std::endl;

    //   std::cout<<" "<<FTris.cV(0)->P().X()<<","<<
    //                   FTris.cV(0)->P().Y()<<","<<
    //                   FTris.cV(0)->P().Z()<<"-"<<
    //                   FTris.cV(1)->P().X()<<","<<
    //                   FTris.cV(1)->P().Y()<<","<<
    //                   FTris.cV(1)->P().Z()<<"-"<<
    //                   FTris.cV(2)->P().X()<<","<<
    //                   FTris.cV(2)->P().Y()<<","<<
    //                   FTris.cV(2)->P().Z()<<std::endl;


    //   std::cout<<" "<<FPoly.cV(0)->cQ()<<" "<<FPoly.cV(1)->cQ()<<
    //              " "<<FPoly.cV(2)->Q()<< " "<<FPoly.cV(3)->Q()<<std::endl;

    //   std::cout<<" "<<FPoly.cV(0)->P().X()<<","<<
    //                   FPoly.cV(0)->P().Y()<<","<<
    //                   FPoly.cV(0)->P().Z()<<"-"<<
    //                   FPoly.cV(1)->P().X()<<","<<
    //                   FPoly.cV(1)->P().Y()<<","<<
    //                   FPoly.cV(1)->P().Z()<<"-"<<
    //                   FPoly.cV(2)->P().X()<<","<<
    //                   FPoly.cV(2)->P().Y()<<","<<
    //                   FPoly.cV(2)->P().Z()<<"-"<<
    //              FPoly.cV(3)->P().X()<<","<<
    //              FPoly.cV(3)->P().Y()<<","<<
    //              FPoly.cV(3)->P().Z()<<std::endl;

    //   std::cout<<" "<<IndexV[0]<<" "<<IndexV[1]<<" "<<IndexV[2]<<std::endl;
    assert((IndexV[0]!=IndexV[1])&&(IndexV[1]!=IndexV[2])&&(IndexV[0]!=IndexV[2]));

    for (size_t i=0;i<3;i++)
        bary.V(i)=std::max(bary.V(i),(ScalarType)0);
    //then set the weights
    for (size_t i=0;i<3;i++)
    {
        assert(bary.V(i)>=0);
        //distribute on all
        if (IndexV[i]==-1)
        {
            VertWeigths[0]+=bary.V(i)/4;
            VertWeigths[1]+=bary.V(i)/4;
            VertWeigths[2]+=bary.V(i)/4;
            VertWeigths[3]+=bary.V(i)/4;
        }
        else
        {
            assert(IndexV[i]>=0);
            assert(IndexV[i]<4);
            VertWeigths[IndexV[i]] =bary.V(i);
        }
    }

    //finally normalize
    ScalarType Sum=0;
    for (size_t i=0;i<4;i++)
        Sum+=VertWeigths[i];

    for (size_t i=0;i<4;i++)
    {
        assert(VertWeigths[i]>=0);
        VertWeigths[i]/=Sum;
    }
}

template <class PolyMeshType,class TriMeshType>
void ProjectStepPositions(PolyMeshType &PolyM,TriMeshType &TriM,
                          vcg::GridStaticPtr<typename TriMeshType::FaceType,typename TriMeshType::ScalarType> TriGrid,
                          const ProjectionBase &TriProjBase,
                          const BasicMesh &EdgeM,
                          const ProjectionBase &PolyProjBase,
                          std::vector<typename PolyMeshType::CoordType> &TargetPos)
{
    typedef typename PolyMeshType::ScalarType ScalarType;
    typedef typename PolyMeshType::CoordType CoordType;
    typedef typename TriMeshType::FaceType TriFaceType;
    TargetPos.clear();

    std::vector<int> SharpProj(PolyM.vert.size(),-1);
    for (size_t i=0;i<PolyProjBase.SharpEdge.size();i++)
        for (size_t j=0;j<PolyProjBase.SharpEdge[i].size();j++)
        {
            size_t IndexV0=PolyProjBase.SharpEdge[i][j].first;
            size_t IndexV1=PolyProjBase.SharpEdge[i][j].second;
            SharpProj[IndexV0]=i;
            SharpProj[IndexV1]=i;
        }

    for (size_t i=0;i<PolyM.vert.size();i++)
    {
        if ((PolyProjBase.VertProjType[i]==ProjCorner)||
                (PolyProjBase.VertProjType[i]==ProjNone))
        {
            TargetPos.push_back(PolyM.vert[i].P());
            //            CoordType TestPos=PolyM.vert[i].P();
            //            int IndexSh=SharpProj[i];
            //            assert(IndexSh>=0);
            //            assert(IndexSh<TriProjBase.SharpEdge.size());
            //            size_t ClosestSeg;
            //            CoordType ClosestPos;
            //            ClosestPointEdgeSet(TestPos,TriM,TriProjBase.SharpEdge[IndexSh],ClosestSeg,ClosestPos);
            //            TargetPos.push_back(ClosestPos);
        }
        if (PolyProjBase.VertProjType[i]==ProjSharp)
        {
            CoordType TestPos=PolyM.vert[i].P();
                        size_t IndexE;
                        ScalarType t,MinD;
                        CoordType Clos;
                        ClosestPointEMesh(TestPos,EdgeM,IndexE,t,MinD,Clos);
                        TargetPos.push_back(Clos);
            //            //std::cout<<"TEST"<<std::endl;
            //
//            int IndexSh=SharpProj[i];
//            assert(IndexSh>=0);
//            assert(IndexSh<TriProjBase.SharpEdge.size());
//            size_t ClosestSeg;
//            CoordType ClosestPos;
//            ClosestPointEdgeSet(TestPos,TriM,TriProjBase.SharpEdge[IndexSh],ClosestSeg,ClosestPos);
//            TargetPos.push_back(ClosestPos);
        }
        if (PolyProjBase.VertProjType[i]==ProjSuface)
        {
            CoordType TestPos=PolyM.vert[i].P();
            CoordType closestPt;
            ScalarType MaxD=PolyM.bbox.Diag();
            ScalarType MinD;
            TriFaceType *f=vcg::tri::GetClosestFaceBase(TriM,TriGrid,TestPos,MaxD,MinD,closestPt);
            assert(f!=NULL);
            TargetPos.push_back(closestPt);
        }
    }
}

//template <class PolyMeshType,class TriMeshType>
//void GetFlippedVert(const PolyMeshType &poly_mesh,
//                    const TriMeshType &tri_mesh,
//                    vcg::GridStaticPtr<typename TriMeshType::FaceType,typename TriMeshType::ScalarType> TriGrid,
//                    std::vector<size_t> &FlippedV)
//{
//    typedef typename PolyMeshType::ScalarType ScalarType;
//    typedef typename PolyMeshType::CoordType CoordType;
//    typedef typename TriMeshType::FaceType TriFaceType;

//    FlippedV.clear();

//    std::vector<CoordType> TargetN;
//    for (size_t i=0;i<poly_mesh.vert.size();i++)
//    {
//        CoordType TestPos=poly_mesh.vert[i].P();
//        CoordType closestPt;
//        ScalarType MaxD=PolyM.bbox.Diag();
//        ScalarType MinD;
//        TriFaceType *f=vcg::tri::GetClosestFaceBase(TriM,TriGrid,TestPos,MaxD,MinD,closestPt);
//        assert(f!=NULL);
//        TargetPos.push_back(closestPt);
//    }
//}

template <class PolyMeshType>
void GetFlippedFaces(const PolyMeshType &poly_mesh,
                     std::vector<size_t> &FlippedF)
{
    typedef typename PolyMeshType::ScalarType ScalarType;
    typedef typename PolyMeshType::CoordType CoordType;
    //    vcg::PolygonalAlgorithm::
    FlippedF.clear();

    for (size_t i=0;i<poly_mesh.face.size();i++)
    {
        std::vector<CoordType> VertN;
        size_t sizeF=poly_mesh.face[i].VN();
        for (size_t j=0;j<poly_mesh.face[i].VN();j++)
        {
            CoordType Pos0=poly_mesh.face[i].V((j+sizeF-1)%sizeF)->P();
            CoordType Pos1=poly_mesh.face[i].V(j)->P();
            CoordType Pos2=poly_mesh.face[i].V((j+1)%sizeF)->P();
            CoordType Dir0=Pos0-Pos1;
            CoordType Dir1=Pos2-Pos1;
            Dir0.Normalize();
            Dir1.Normalize();
            CoordType Norm=Dir0^Dir1;
            Norm.Normalize();
            VertN.push_back(Norm);
        }
        bool isFlipped=false;
        for (size_t j=0;j<VertN.size();j++)
        {
            CoordType Norm0=VertN[j];
            CoordType Norm1=VertN[(j+1)%VertN.size()];
            if ((Norm0*Norm1)<0)isFlipped=true;
        }
        if (isFlipped)
            FlippedF.push_back(i);
    }

}

template <class PolyMeshType>
void IsFlippedFaces(const PolyMeshType &poly_mesh,
                     std::vector<bool> &IsFlipped)
{
    std::vector<size_t> FlippedF;
    IsFlipped=std::vector<bool>(poly_mesh.face.size(),false);
    GetFlippedFaces(poly_mesh,FlippedF);
    for (size_t i=0;i<FlippedF.size();i++)
        IsFlipped[FlippedF[i]]=true;
}

template <class PolyMeshType,class TriMeshType>
void GetProjectionBasis(const TriMeshType &tri_mesh,
                        const std::vector<std::pair<size_t,size_t> > &FeatureTris,
                        const std::vector<size_t> &FeatureTrisC,
                        const std::vector<size_t> &tri_face_partition,
                        PolyMeshType &poly_mesh,
                        const std::vector<size_t > &poly_corner,
                        const std::vector<size_t> &poly_face_partition,
                        const typename PolyMeshType::ScalarType &AvEdge,
                        ProjectionBase &PBaseTris,
                        ProjectionBase &PBasePoly)
{
    //    typedef typename PolyMeshType::ScalarType ScalarType;
    //    typedef typename PolyMeshType::CoordType CoordType;
    //    typedef typename vcg::Segment3<typename PolyMeshType::ScalarType> SegmentType;

    PBaseTris.Clear();
    PBasePoly.Clear();

    PBaseTris.VertProjType.resize(tri_mesh.vert.size(),ProjSuface);
    PBasePoly.VertProjType.resize(poly_mesh.vert.size(),ProjSuface);


    //get for each corner of the tri mesh the corresponding on the poly mesh
    PBaseTris.Corner=FeatureTrisC;
    GetCorners(FeatureTris,PBaseTris.Corner,tri_mesh,poly_corner,poly_mesh,AvEdge,PBasePoly.Corner);


    //find each segment set per pair of partitions
    typedef std::pair<int,int> PatchPairKey;
    std::map<PatchPairKey,size_t > BasisMap;
    std::set<std::pair<size_t,size_t> > InsertedTriEdges;
    for (size_t i=0;i<FeatureTris.size();i++)
    {
        size_t IndexF=FeatureTris[i].first;
        size_t IndexE=FeatureTris[i].second;
        int Part0=tri_face_partition[IndexF];
        int Part1=-1;//when is adjacent to border
        if (!vcg::face::IsBorder(tri_mesh.face[IndexF],IndexE))
        {
            size_t IndexFOpp=vcg::tri::Index(tri_mesh,tri_mesh.face[IndexF].cFFp(IndexE));
            Part1=tri_face_partition[IndexFOpp];
        }

        //see if there is such configuration of adjacent patches
        assert(Part0!=Part1);
        PatchPairKey keyPatch(std::min(Part0,Part1),std::max(Part0,Part1));


        //        CoordType Pos0=tri_mesh.face[IndexF].cP0(IndexE);
        //        CoordType Pos1=tri_mesh.face[IndexF].cP1(IndexE);
        //        SegmentType SBasis(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
        size_t IndexV0=vcg::tri::Index(tri_mesh,tri_mesh.face[IndexF].cV0(IndexE));
        size_t IndexV1=vcg::tri::Index(tri_mesh,tri_mesh.face[IndexF].cV1(IndexE));
        std::pair<size_t,size_t> EdgeKey=std::pair<size_t,size_t>(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));

        PBaseTris.VertProjType[IndexV0]=ProjSharp;
        PBaseTris.VertProjType[IndexV1]=ProjSharp;

        //check if already added
        if (InsertedTriEdges.count(EdgeKey)>0)continue;
        InsertedTriEdges.insert(EdgeKey);

        if (BasisMap.count(keyPatch)==0)
        {
            //            ProjBasis.resize(ProjBasis.size()+1);
            //            ProjBasis.back().push_back(SBasis);
            PBaseTris.SharpEdge.resize(PBaseTris.SharpEdge.size()+1);
            PBaseTris.SharpEdge.back().push_back(EdgeKey);
            BasisMap[keyPatch]=PBaseTris.SharpEdge.size()-1;
        }
        else
        {
            int IndexSharp=BasisMap[keyPatch];
            PBaseTris.SharpEdge[IndexSharp].push_back(EdgeKey);
            //ProjBasis[IndexP].push_back(SBasis);
        }
    }

    //    //std::cout<<"B"<<std::endl;

    //    //then check the edge of the quads, by default =-2, no projection, internal
    //    VertBasis.resize(poly_mesh.vert.size(),-2);
    //    //    BasicMesh EdgeMesh;
    //    //    ExtractEdgeMesh(tri_mesh,features,EdgeMesh);

    PBasePoly.SharpEdge.resize(PBaseTris.SharpEdge.size());

    std::set<std::pair<size_t,size_t> > InsertedPolyEdges;
    for (size_t i=0;i<poly_mesh.face.size();i++)
    {
        int Part0=poly_face_partition[i];

        //get neighbours
        for (size_t j=0;j<poly_mesh.face[i].VN();j++)
        {
            //default border
            int Part1=-1;
            //std::cout<<"i j "<<i<<","<<j<<std::endl;

            //            int sizeF=poly_mesh.face[i].VN();
            //            size_t IndexVa0=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V(j));
            //            size_t IndexVa1=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V((j+1)%sizeF));
            //            std::pair<size_t,size_t> edgeKey(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
            //std::cout<<"Va0 Va1 "<<IndexVa0<<","<<IndexVa1<<std::endl;

            if (!vcg::face::IsBorder(poly_mesh.face[i],j))
            {
                size_t IndexFOpp=vcg::tri::Index(poly_mesh,poly_mesh.face[i].cFFp(j));
                Part1=poly_face_partition[IndexFOpp];
            }
            std::pair<int,int> keyPatch(std::min(Part0,Part1),std::max(Part0,Part1));
            //check if this is part of a sharp feature
            if (BasisMap.count(keyPatch)>0)
            {
                if (Part0==Part1)std::cout<<"WARNING: Internal SHARP FEATURE of a patch"<<std::endl;
                int sizeF=poly_mesh.face[i].VN();
                size_t IndexV0=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V(j));
                size_t IndexV1=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V((j+1)%sizeF));
                std::pair<size_t,size_t> edgeKey(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
                //                std::cout<<"V0 V1 "<<IndexV0<<","<<IndexV1<<std::endl;

                //                std::cout<<"WTest 0"<<std::endl;

                if (InsertedPolyEdges.count(edgeKey)>0)continue;

                //                std::cout<<"WTest 1"<<std::endl;

                InsertedPolyEdges.insert(edgeKey);

                //  ClosestPointEMesh(poly_mesh.vert[IndV0].cP(),edge_mesh,IndexE,t,MinD,Clos);
                //  if (MinD>AvEdge/10)continue;

                size_t IndexSharp=BasisMap[keyPatch];
                assert(IndexSharp<PBasePoly.SharpEdge.size());
                assert(IndexSharp<PBaseTris.SharpEdge.size());

                PBasePoly.SharpEdge[IndexSharp].push_back(edgeKey);
                PBasePoly.VertProjType[IndexV0]=ProjSharp;
                PBasePoly.VertProjType[IndexV1]=ProjSharp;

                //                VertBasis[IndV0]=BasisMap[key];
                //                VertBasis[IndV1]=BasisMap[key];
            }
        }
    }

    //    //std::cout<<"C"<<std::endl;
    //    for (size_t i=0;i<cornersIdx.size();i++)
    //        VertBasis[cornersIdx[i]]=-1;

    //    if (ProjBasis.size()>0)
    //        MergeProjBasis(tri_mesh,featuresC,VertBasis,ProjBasis);
    //    //std::cout<<"D"<<std::endl;

    //finally set the corners
    PBaseTris.Corner=FeatureTrisC;
    GetCorners(FeatureTris,PBaseTris.Corner,tri_mesh,poly_corner,poly_mesh,AvEdge,PBasePoly.Corner);
    for (size_t i=0;i<PBaseTris.Corner.size();i++)
    {
        size_t IndexC=PBaseTris.Corner[i];
        PBaseTris.VertProjType[IndexC]=ProjCorner;
    }
    for (size_t i=0;i<PBasePoly.Corner.size();i++)
    {
        size_t IndexC=PBasePoly.Corner[i];
        PBasePoly.VertProjType[IndexC]=ProjCorner;
    }

//    //make vertices equivalent
//    for (size_t i=0;i<PBasePoly.Corner.size();i++)
//    {
//        size_t IndexCPoly=PBasePoly.Corner[i];
//        size_t IndexCTri=PBaseTris.Corner[i];
//        poly_mesh.vert[IndexCPoly].P()=tri_mesh.vert[IndexCTri].P();
//    }
}

template <class PolyMeshType,class TriMeshType>
void InitPolyTrisMesh(PolyMeshType &PolyM,TriMeshType &poly_tris)
{
    //set the index of original vertices
    for (size_t i=0;i<PolyM.vert.size();i++)
        PolyM.vert[i].Q()=i;

    //and the original faces
    for (size_t i=0;i<PolyM.face.size();i++)
        PolyM.face[i].Q()=i;

    vcg::PolygonalAlgorithm<PolyMeshType>::TriangulateToTriMesh(PolyM,poly_tris);
    assert(poly_tris.vert.size()==(PolyM.vert.size()+PolyM.face.size()));
    assert(poly_tris.face.size()==(4*PolyM.face.size()));

    //then update values on the mesh
    vcg::tri::UpdateNormal<TriMeshType>::PerFaceNormalized(poly_tris);
    vcg::tri::UpdateNormal<TriMeshType>::PerVertexNormalized(poly_tris);
    vcg::tri::UpdateBounding<TriMeshType>::Box(poly_tris);

    //check the index of the polygonal face
    for (size_t i=0;i<poly_tris.vert.size();i++)
    {
        assert(poly_tris.face[i].Q()>=0);
        assert(poly_tris.face[i].Q()<PolyM.face.size());
    }

    //set the index of the original vertices
    for (size_t i=0;i<poly_tris.vert.size();i++)
    {
        poly_tris.vert[i].Q()=-1;
        if (i<PolyM.vert.size())
            poly_tris.vert[i].Q()=i;
    }
}

template <class PolyMeshType,class TriMeshType>
void GetMovingPointOnSurface(PolyMeshType &PolyM,TriMeshType &poly_tris,
                             const typename TriMeshType::CoordType &TestPos,
                             vcg::GridStaticPtr<typename TriMeshType::FaceType,typename TriMeshType::ScalarType> &Poly_tri_Grid,
                             std::vector<std::vector<typename TriMeshType::CoordType> > &VertMove,
                             std::vector<std::vector<typename TriMeshType::ScalarType> > &VertWeight)
{
    typedef typename PolyMeshType::ScalarType ScalarType;
    typedef typename PolyMeshType::CoordType CoordType;
    typedef typename TriMeshType::FaceType TriFaceType;

    assert(VertMove.size()==PolyM.vert.size());
    assert(VertWeight.size()==PolyM.vert.size());

    CoordType closestPt;
    ScalarType MaxD=PolyM.bbox.Diag();
    ScalarType MinD;
    TriFaceType *f=vcg::tri::GetClosestFaceBase(poly_tris,Poly_tri_Grid,TestPos,MaxD,MinD,closestPt);
    //assert(f!=NULL);

    //retrieve the original face
    int IndexTriF=vcg::tri::Index(poly_tris,f);
    int IndexPolyF=poly_tris.face[IndexTriF].Q();

    assert(IndexPolyF>=0);
    assert(IndexPolyF<PolyM.face.size());

    std::vector<typename TriMeshType::ScalarType> VertWeigths;
    //std::cout<<"A"<<std::endl;
    GetQuadInterpW<PolyMeshType,TriMeshType>(poly_tris.face[IndexTriF],
                                             PolyM.face[IndexPolyF],
                                             closestPt,VertWeigths);
    //std::cout<<"B"<<std::endl;
    CoordType MoveVect=TestPos-closestPt;

    for (size_t j=0;j<4;j++)
    {
        int IndexV=vcg::tri::Index(PolyM,PolyM.face[IndexPolyF].V(j));
        VertMove[IndexV].push_back(MoveVect*VertWeigths[j]);
        VertWeight[IndexV].push_back(VertWeigths[j]);
    }
}

//template <class PolyMeshType,class TriMeshType>
//void BackProjectStepPositions(PolyMeshType &PolyM,TriMeshType &TriM,
//                              const ProjectionBase &TriProjBase,
//                              const ProjectionBase &PolyProjBase,
//                              std::vector<typename PolyMeshType::CoordType> &TargetPos)
//{
//    typedef typename PolyMeshType::ScalarType ScalarType;
//    typedef typename PolyMeshType::CoordType CoordType;
//    typedef typename TriMeshType::FaceType TriFaceType;

//    assert(TriProjBase.VertProjType.size()==TriM.vert.size());
//    assert(PolyProjBase.VertProjType.size()==PolyM.vert.size());
//    assert(TriProjBase.SharpEdge.size()==PolyProjBase.SharpEdge.size());
//    assert(PolyProjBase.Corner.size()==PolyProjBase.Corner.size());

//    //transform the polygonal mesh into triangle one
//    TriMeshType poly_tris;

//    //    //set the index of original vertices
//    //    for (size_t i=0;i<PolyM.vert.size();i++)
//    //        PolyM.vert[i].Q()=i;

//    //    //and the original faces
//    //    for (size_t i=0;i<PolyM.face.size();i++)
//    //        PolyM.face[i].Q()=i;

//    //    vcg::PolygonalAlgorithm<PolyMeshType>::TriangulateToTriMesh(PolyM,poly_tris);
//    //    assert(poly_tris.vert.size()==(PolyM.vert.size()+PolyM.face.size()));
//    //    assert(poly_tris.face.size()==(4*PolyM.face.size()));

//    //    //then update values on the mesh
//    //    vcg::tri::UpdateNormal<TriMeshType>::PerFaceNormalized(poly_tris);
//    //    vcg::tri::UpdateNormal<TriMeshType>::PerVertexNormalized(poly_tris);
//    //    vcg::tri::UpdateBounding<TriMeshType>::Box(poly_tris);

//    //    //check the index of the polygonal face
//    //    for (size_t i=0;i<poly_tris.vert.size();i++)
//    //    {
//    //       assert(poly_tris.face[i].Q()>=0);
//    //       assert(poly_tris.face[i].Q()<PolyM.face.size());
//    //    }

//    //    //set the index of the original vertices
//    //    for (size_t i=0;i<poly_tris.vert.size();i++)
//    //    {
//    //        poly_tris.vert[i].Q()=-1;
//    //        if (i<PolyM.vert.size())
//    //            poly_tris.vert[i].Q()=i;
//    //    }

//    InitPolyTrisMesh(PolyM,poly_tris);


//    std::vector<std::vector<CoordType> > VertMove(PolyM.vert.size());
//    std::vector<std::vector<ScalarType> > VertWeight(PolyM.vert.size());

//    //first set the internal
//    //then initialize the grid
//    vcg::GridStaticPtr<TriFaceType,ScalarType> TriGrid;
//    TriGrid.Set(poly_tris.face.begin(),poly_tris.face.end());

//    for (size_t i=0;i<TriM.vert.size();i++)
//    {
//        if (TriProjBase.VertProjType[i]!=ProjSuface)continue;

//        CoordType TestPos=TriM.vert[i].P();
//        GetMovingPointOnSurface(PolyM,poly_tris,TestPos,TriGrid,VertMove,VertWeight);

//        //        CoordType closestPt;
//        //        ScalarType MaxD=PolyM.bbox.Diag();
//        //        ScalarType MinD;
//        //        TriFaceType *f=vcg::tri::GetClosestFaceBase(poly_tris,TriGrid,TestPos,MaxD,MinD,closestPt);
//        //        assert(f!=NULL);

//        //        //retrieve the original face
//        //        int IndexTriF=vcg::tri::Index(poly_tris,f);
//        //        int IndexPolyF=poly_tris.face[IndexTriF].Q();

//        //        assert(IndexPolyF>=0);
//        //        assert(IndexPolyF<PolyM.face.size());

//        //        std::vector<typename TriMeshType::ScalarType> VertWeigths;
//        //        //std::cout<<"A"<<std::endl;
//        //        GetQuadInterpW<PolyMeshType,TriMeshType>(poly_tris.face[IndexTriF],
//        //                                                 PolyM.face[IndexPolyF],
//        //                                                 closestPt,VertWeigths);
//        //        //std::cout<<"B"<<std::endl;
//        //        CoordType MoveVect=TestPos-closestPt;

//        //        for (size_t j=0;j<4;j++)
//        //        {
//        //            int IndexV=vcg::tri::Index(PolyM,PolyM.face[IndexPolyF].V(j));
//        //            VertMove[IndexV].push_back(MoveVect*VertWeigths[j]);
//        //            VertWeight[IndexV].push_back(VertWeigths[j]);
//        //        }

//    }

//    //normalize the weight and average the direction
//    std::vector<CoordType> TargetMov(PolyM.vert.size(),CoordType(0,0,0));
//    for (size_t i=0;i<VertWeight.size();i++)
//    {
//        ScalarType SumW=0;
//        for (size_t j=0;j<VertWeight[i].size();j++)
//            SumW+=VertWeight[i][j];

//        if (SumW==0)continue;

//        for (size_t j=0;j<VertWeight[i].size();j++)
//            VertWeight[i][j]/=SumW;

//        for (size_t j=0;j<VertMove[i].size();j++)
//        {
//            TargetMov[i]+=VertMove[i][j]*VertWeight[i][j];
//        }
//    }

//    TargetPos.clear();
//    for (size_t i=0;i<PolyM.vert.size();i++)
//    {
//        if (PolyProjBase.VertProjType[i]==ProjCorner)
//            TargetPos.push_back(PolyM.vert[i].P());
//        else
//            TargetPos.push_back(PolyM.vert[i].P()+TargetMov[i]);
//    }
//}

template <class PolyMeshType,class TriMeshType>
void BackProjectStepPositions(PolyMeshType &PolyM,TriMeshType &TriM,
                              const ProjectionBase &TriProjBase,
                              const ProjectionBase &PolyProjBase,
                              std::vector<typename PolyMeshType::CoordType> &TargetPos)
{
    typedef typename PolyMeshType::ScalarType ScalarType;
    typedef typename PolyMeshType::CoordType CoordType;
    typedef typename TriMeshType::FaceType TriFaceType;

    assert(TriProjBase.VertProjType.size()==TriM.vert.size());
    assert(PolyProjBase.VertProjType.size()==PolyM.vert.size());
    assert(TriProjBase.SharpEdge.size()==PolyProjBase.SharpEdge.size());
    assert(PolyProjBase.Corner.size()==PolyProjBase.Corner.size());

    //transform the polygonal mesh into triangle one
    TriMeshType poly_tris;

    InitPolyTrisMesh(PolyM,poly_tris);

    std::vector<std::vector<CoordType> > VertMove(PolyM.vert.size());
    std::vector<std::vector<ScalarType> > VertWeight(PolyM.vert.size());

    //first set the internal
    //then initialize the grid
    vcg::GridStaticPtr<TriFaceType,ScalarType> TriGrid;
    TriGrid.Set(poly_tris.face.begin(),poly_tris.face.end());

    for (size_t i=0;i<TriM.vert.size();i++)
    {
        if (TriProjBase.VertProjType[i]!=ProjSuface)continue;

        CoordType TestPos=TriM.vert[i].P();
        GetMovingPointOnSurface(PolyM,poly_tris,TestPos,TriGrid,VertMove,VertWeight);
    }

    //normalize the weight and average the direction
    std::vector<CoordType> TargetMov(PolyM.vert.size(),CoordType(0,0,0));
    for (size_t i=0;i<VertWeight.size();i++)
    {
        ScalarType SumW=0;
        for (size_t j=0;j<VertWeight[i].size();j++)
            SumW+=VertWeight[i][j];

        if (SumW==0)continue;

        for (size_t j=0;j<VertWeight[i].size();j++)
            VertWeight[i][j]/=SumW;

        for (size_t j=0;j<VertMove[i].size();j++)
        {
            TargetMov[i]+=VertMove[i][j]*VertWeight[i][j];
        }
    }

    TargetPos.clear();
    for (size_t i=0;i<PolyM.vert.size();i++)
    {
        if ((PolyProjBase.VertProjType[i]==ProjCorner)||
                (PolyProjBase.VertProjType[i]==ProjNone))
            TargetPos.push_back(PolyM.vert[i].P());
        else
            TargetPos.push_back(PolyM.vert[i].P()+TargetMov[i]);
    }
}

template <class PolyMeshType,class TriMeshType>
void MultiCostraintSmoothStep(PolyMeshType &PolyM,TriMeshType &TriM,
                              const BasicMesh &EdgeM,
                              vcg::GridStaticPtr<typename TriMeshType::FaceType,typename TriMeshType::ScalarType> &TriGrid,
                              ProjectionBase &TriProjBase,ProjectionBase &PolyProjBase,
                              SmoothType SType,
                              const typename PolyMeshType::ScalarType Damp,
                              size_t back_proj_steps)
{
    typedef typename PolyMeshType::CoordType CoordType;

    std::vector<CoordType> TargetPosProj;
    std::vector<CoordType> TargetPosBackProj;
    std::vector<CoordType> TargetPosSmooth;

//    //save all old positions
//    std::vector<bool> IsFlipped0;
//    IsFlippedFaces(PolyM,IsFlipped0);

//    std::vector<CoordType> OldPos;
//    for (size_t i=0;i<PolyM.vert.size();i++)
//        OldPos.push_back(PolyM.vert[i].P());

    if (SType==Laplacian)
        LaplacianEdgePos(PolyM,PolyProjBase,Damp,TargetPosSmooth);
    else
        TemplatePolyPos(PolyM,PolyProjBase,Damp,TargetPosSmooth);

    //then blend between the two
    for (size_t i=0;i<PolyM.vert.size();i++)
        PolyM.vert[i].P()=TargetPosSmooth[i];

    for (size_t i=0;i<back_proj_steps;i++)
    {
        BackProjectStepPositions(PolyM,TriM,TriProjBase,PolyProjBase,TargetPosBackProj);
        for (size_t i=0;i<PolyM.vert.size();i++)
            PolyM.vert[i].P()=TargetPosBackProj[i];
    }


    ProjectStepPositions(PolyM,TriM,TriGrid,TriProjBase,EdgeM,PolyProjBase,TargetPosProj);

    //ProjectStepPositions(PolyM,TriM,TriGrid,TriProjBase,PolyProjBase,TargetPosProj);
    for (size_t i=0;i<PolyM.vert.size();i++)
        PolyM.vert[i].P()=TargetPosProj[i];

//    //check the flip
//    std::vector<bool> IsFlipped1;
//    IsFlippedFaces(PolyM,IsFlipped1);
//    for (size_t i=0;i<PolyM.face.size();i++)
//    {
//        if ((!IsFlipped0[i])&&(IsFlipped1[i]))
//        for (size_t j=0;j<PolyM.face[i].VN();j++)
//        {
//            size_t IndexV=vcg::tri::Index(PolyM,PolyM.face[i].V(j));
//            PolyM.vert[IndexV].P()=OldPos[IndexV];
//        }
//    }
    //        //then blend between the two
    //        for (size_t i=0;i<PolyM.vert.size();i++)
    //            PolyM.vert[i].P()=TargetPosSmooth[i]*0.5+TargetPosProj[i]*0.5;
}

template <class PolyMeshType>
void DeCostrainFace(const PolyMeshType &PolyM,
                    const size_t &IndexF,
                    ProjectionBase &PolyProjBase)
{
    bool HasModified=false;
    for (size_t i=0;i<PolyM.face[IndexF].VN();i++)
    {
        size_t IndexV=vcg::tri::Index(PolyM,PolyM.face[IndexF].V(i));
        if (PolyProjBase.VertProjType[IndexV]==ProjCorner)
        {
            PolyProjBase.VertProjType[IndexV]=ProjSharp;
            HasModified=true;
        }
    }
    if (HasModified)return;
    for (size_t i=0;i<PolyM.face[IndexF].VN();i++)
    {
        size_t IndexV=vcg::tri::Index(PolyM,PolyM.face[IndexF].V(i));
        if (PolyProjBase.VertProjType[IndexV]==ProjSharp)
        {
            PolyProjBase.VertProjType[IndexV]=ProjSuface;
            HasModified=true;
        }
    }
    if (HasModified)return;
    for (size_t i=0;i<PolyM.face[IndexF].VN();i++)
    {
        size_t IndexV=vcg::tri::Index(PolyM,PolyM.face[IndexF].V(i));
        PolyProjBase.VertProjType[IndexV]=ProjNone;
    }

}

template <class PolyMeshType>
void DeCostrainFace(const PolyMeshType &PolyM,
                    const std::vector<size_t> &IndexF,
                    ProjectionBase &PolyProjBase)
{
    for (size_t i=0;i<IndexF.size();i++)
        DeCostrainFace(PolyM,IndexF[i],PolyProjBase);
}

template <class PolyMeshType,class TriMeshType>
void MultiCostraintSmooth(PolyMeshType &PolyM,
                          TriMeshType &TriM,
                          const std::vector<std::pair<size_t,size_t> > &features,
                          const std::vector<size_t> &featuresC,
                          const std::vector<size_t> &tri_face_partition,
                          const std::vector<size_t > &quad_corner,
                          const std::vector<size_t> &quad_face_partition,
                          //SmoothType SType,
                          const typename PolyMeshType::ScalarType Damp,
                          const typename PolyMeshType::ScalarType AvEdge,
                          size_t step_num,
                          size_t back_proj_steps)
{
    typedef typename PolyMeshType::ScalarType ScalarType;
    typedef typename PolyMeshType::CoordType CoordType;
    typedef typename TriMeshType::FaceType TriFaceType;

    std::cout<<"*** Getting Projection basis ***"<<std::endl;
    ProjectionBase TriProjBase,PolyProjBase;

    GetProjectionBasis(TriM,features,featuresC,tri_face_partition,PolyM,
                       quad_corner,quad_face_partition,AvEdge,
                       TriProjBase,PolyProjBase);

    //    std::vector<int> VertBasis;
    //    std::vector<std::vector<vcg::Segment3<typename PolyMeshType::ScalarType> > > ProjBasis;
    //    GetVertProjBasis(AvEdge,features,featuresC,TriM,quad_corner,PolyM,tri_face_partition,quad_face_partition,VertBasis,ProjBasis);
    //    std::cout<<"*** End Getting Projection basis ***"<<std::endl;

    //    //select corners
    //    vcg::tri::UpdateSelection<PolyMeshType>::VertexClear(PolyM);
    //    for (size_t i=0;i<VertBasis.size();i++)
    //        if (VertBasis[i]==-1)
    //            PolyM.vert[i].SetS();

    BasicMesh EdgeM;
    ExtractEdgeMesh(TriM,features,EdgeM);


    vcg::GridStaticPtr<TriFaceType,ScalarType> TriGrid;
    TriGrid.Set(TriM.face.begin(),TriM.face.end());
    for (size_t s=0;s<step_num/2;s++)
        MultiCostraintSmoothStep<PolyMeshType,TriMeshType>(PolyM,TriM,EdgeM,TriGrid,TriProjBase,PolyProjBase,Laplacian,Damp,back_proj_steps);
//    for (size_t s=0;s<step_num/2;s++)
//        MultiCostraintSmoothStep<PolyMeshType,TriMeshType>(PolyM,TriM,EdgeM,TriGrid,TriProjBase,PolyProjBase,TemplateFit,Damp,back_proj_steps);


    std::vector<size_t> FlippedF;
    GetFlippedFaces(PolyM,FlippedF);
    for (size_t i=0;i<FlippedF.size();i++)
        PolyM.face[FlippedF[i]].C()=vcg::Color4b(255,0,0,255);//true;
//    int MaxSteps=20;
//    int currSteps=0;
//    std::cout<<"There are  "<<FlippedF.size()<<" flipped faces "<<std::endl;
//    while ((FlippedF.size()>0)&&(currSteps<MaxSteps))
//    {
//        DeCostrainFace(PolyM,FlippedF,PolyProjBase);
//        MultiCostraintSmoothStep<PolyMeshType,TriMeshType>(PolyM,TriM,EdgeM,TriGrid,TriProjBase,PolyProjBase,Laplacian,Damp,back_proj_steps);
//        GetFlippedFaces(PolyM,FlippedF);

//        currSteps++;
//        std::cout<<"There are  "<<FlippedF.size()<<" flipped faces "<<std::endl;
//    }

}



template <class PolyMeshType> typename PolyMeshType::CoordType ClosestPointSegSet(const typename PolyMeshType::CoordType &Pos,
                                                                                  std::vector<vcg::Segment3<typename PolyMeshType::ScalarType> > &SegSet)
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
void ClosestPointEdgeSet(const typename PolyMeshType::CoordType &Pos,
                         const std::vector<vcg::Segment3<typename PolyMeshType::ScalarType> > &SegSet,
                         size_t &ClosestSeg,typename PolyMeshType::CoordType &Clos)
{
    typedef typename PolyMeshType::ScalarType ScalarType;
    typedef typename PolyMeshType::CoordType CoordType;
    ScalarType MinD=std::numeric_limits<ScalarType>::max();
    ClosestSeg=0;
    for (size_t i=0;i<SegSet.size();i++)
    {
        vcg::Segment3<ScalarType> STest=SegSet[i];
        ScalarType testD;
        CoordType ClosTest;
        vcg::SegmentPointDistance(STest,Pos,ClosTest,testD);
        if (testD>MinD)continue;
        ClosestSeg=i;
        Clos=ClosTest;
        MinD=testD;
    }
}

template <class MeshType>
void ClosestPointEdgeSet(const typename MeshType::CoordType &Pos,
                         const MeshType &sampleMesh,
                         const std::vector<std::pair<size_t,size_t> > &EdgeSet,
                         size_t &ClosestSeg,
                         typename MeshType::CoordType &Clos)
{
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;

    std::vector<vcg::Segment3<typename MeshType::ScalarType> > SegSet;
    for (size_t i=0;i<EdgeSet.size();i++)
    {
        size_t IndexV0=EdgeSet[i].first;
        size_t IndexV1=EdgeSet[i].second;
        assert(IndexV0!=IndexV1);
        assert(IndexV0<sampleMesh.vert.size());
        assert(IndexV1<sampleMesh.vert.size());
        CoordType Pos0=sampleMesh.vert[IndexV0].cP();
        CoordType Pos1=sampleMesh.vert[IndexV1].cP();
        vcg::Segment3<ScalarType> seg(Pos0,Pos1);
        SegSet.push_back(seg);
    }

    ClosestPointEdgeSet<MeshType>(Pos,SegSet,ClosestSeg,Clos);
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

//void InitProejct(TriangleMeshType &tri_mesh,
//                        PolyMeshType &poly_mesh,
//                        const std::vector<std::pair<size_t,size_t> > &features,
//                        const std::vector<size_t> &featuresC,
//                        const std::vector<size_t> &tri_face_partition,
//                        const std::vector<size_t > &quad_corner,
//                        const std::vector<size_t> &quad_face_partition,
//                        SmoothType SType,
//                        size_t Steps,
//                        typename PolyMeshType::ScalarType Damp,
//                        typename PolyMeshType::ScalarType AvEdge)
//{
////InitByCurvature(MeshType & mesh,
////                                unsigned Nring,
////                                bool UpdateFaces=true)
//}

//template<class TriangleMeshType>
//void ProjectMovingDirection(vcg::GridStaticPtr<typename TriangleMeshType::FaceType,
//                                               typename TriangleMeshType::ScalarType> &TriGrid,
//                            TriangleMeshType &tri_mesh,
//                            const typename TriangleMeshType::CoordType &SamplePos,
//                            typename TriangleMeshType::CoordType &MovingDir)
//{
//    typedef typename TriangleMeshType::FaceType TriFaceType;
//    typedef typename TriangleMeshType::ScalarType TriScalarType;
//    typedef typename TriangleMeshType::CoordType CoordType;

//    //get the closest face point
//    CoordType closestPt,Normf,bary;
//    TriScalarType MaxD=tri_mesh.bbox.Diag();
//    TriScalarType MinD;
//    TriFaceType *f=vcg::tri::GetClosestFaceBase(tri_mesh,TriGrid,SamplePos,MaxD,MinD,closestPt,Normf,bary);

//    //project moving direction cannot move along normal
//    MovingDir-=Normf*(MovingDir*Normf);
//}

template<class TriangleMeshType>
void ProjectMovingDirectionStep(vcg::GridStaticPtr<typename TriangleMeshType::FaceType,
                                typename TriangleMeshType::ScalarType> &TriGrid,
                                TriangleMeshType &tri_mesh,
                                const typename TriangleMeshType::CoordType &SamplePos,
                                typename TriangleMeshType::CoordType &MovingDir)
{
    typedef typename TriangleMeshType::FaceType TriFaceType;
    typedef typename TriangleMeshType::ScalarType TriScalarType;
    typedef typename TriangleMeshType::CoordType CoordType;

    //get the closest face point
    CoordType closestPt,Normf,bary;
    TriScalarType MaxD=tri_mesh.bbox.Diag();
    TriScalarType MinD;
    TriFaceType *f=vcg::tri::GetClosestFaceBase(tri_mesh,TriGrid,SamplePos,MaxD,MinD,closestPt,Normf,bary);

    //project moving direction cannot move along normal
    MovingDir-=Normf*(MovingDir*Normf);

    //then get the two directions
    CoordType PD1=f->PD1();
    CoordType PD2=f->PD2();
    PD1.Normalize();
    PD2.Normalize();

    TriScalarType Mag1=fabs(f->K1());//0.9
    TriScalarType Mag2=fabs(f->K2());//0.1
    //    TriScalarType SumMag=Mag1+Mag2;
    //    Mag1/=SumMag;
    //    Mag2/=SumMag;
    //    Mag1=1-Mag1;//0.1
    //    Mag2=1-Mag2;//0.9

    CoordType Dir=MovingDir;
    Dir.Normalize();
    TriScalarType Dot1=fabs(Dir*PD1);//0.1
    TriScalarType Dot2=fabs(Dir*PD2);//0.9

    //TriScalarType InterpDir=Dot1/(Dot1+Dot2);//0.9
    TriScalarType InterpMag=Mag1*Dot1+Mag2*Dot2;//0.1*0.1+0.9*0.9
    MovingDir*=InterpMag;
}

template<class TriangleMeshType>
void ProjectMovingDirection(vcg::GridStaticPtr<typename TriangleMeshType::FaceType,
                            typename TriangleMeshType::ScalarType> &TriGrid,
                            TriangleMeshType &tri_mesh,
                            const typename TriangleMeshType::CoordType &SamplePos,
                            typename TriangleMeshType::CoordType &MovingDir,
                            size_t Steps=10)
{
    typedef typename TriangleMeshType::FaceType TriFaceType;
    typedef typename TriangleMeshType::ScalarType TriScalarType;
    typedef typename TriangleMeshType::CoordType CoordType;

    CoordType CurrStep=MovingDir/Steps;
    CoordType CurrPos=SamplePos;
    for (size_t i=0;i<Steps;i++)
    {
        //make a step
        ProjectMovingDirectionStep(TriGrid,tri_mesh,CurrPos,CurrStep);
        //update current positions
        CurrPos+=CurrStep;
    }

    //the return the final position
    MovingDir=CurrPos-SamplePos;
}

//template <class PolyMeshType,class TriangleMeshType>
//void SmoothWithFeatures(TriangleMeshType &tri_mesh,
//                        PolyMeshType &poly_mesh,
//                        const std::vector<std::pair<size_t,size_t> > &features,
//                        const std::vector<size_t> &featuresC,
//                        const std::vector<size_t> &tri_face_partition,
//                        const std::vector<size_t > &quad_corner,
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

//    std::cout<<"*** Getting Projection basis ***"<<std::endl;
//    std::vector<int> VertBasis;
//    std::vector<std::vector<vcg::Segment3<typename PolyMeshType::ScalarType> > > ProjBasis;
//    GetVertProjBasis(AvEdge,features,featuresC,tri_mesh,quad_corner,poly_mesh,tri_face_partition,
//                     quad_face_partition,VertBasis,ProjBasis);
//    std::cout<<"*** End Getting Projection basis ***"<<std::endl;


//    //    BasicMesh EdgeMesh;
//    //    ExtractEdgeMesh(tri_mesh,features,EdgeMesh);

//    //vcg::tri::io::ExporterOBJ<BasicMesh>::Save(EdgeMesh,"test_edge.obj",vcg::tri::io::Mask::IOM_EDGEINDEX);

//    //select corners
//    vcg::tri::UpdateSelection<PolyMeshType>::VertexClear(poly_mesh);
//    for (size_t i=0;i<VertBasis.size();i++)
//        if (VertBasis[i]==-1)
//            poly_mesh.vert[i].SetS();

//    vcg::GridStaticPtr<TriFaceType,TriScalarType> TriGrid;
//    TriGrid.Set(tri_mesh.face.begin(),tri_mesh.face.end());

//    //    for (size_t i=0;i<poly_mesh.vert.size();i++)
//    //    {

//    //        TriScalarType MinD;
//    //        TriScalarType MaxD=poly_mesh.bbox.Diag();
//    //        CoordType closestPt;
//    //        vcg::tri::GetClosestFaceBase(tri_mesh,TriGrid,poly_mesh.vert[i].P(),MaxD,MinD,closestPt);
//    //        poly_mesh.vert[i].P()=closestPt;

//    //    }

//    //    return;

//    for (size_t s=0;s<Steps;s++)
//    {
//        //save the initial position
//        std::vector<CoordType> Pos0;
//        for (size_t i=0;i<poly_mesh.vert.size();i++)
//            Pos0.push_back(poly_mesh.vert[i].P());

//        if (SType==Laplacian)
//            LaplacianEdge(poly_mesh,true,Damp);
//        //vcg::PolygonalAlgorithm<PolyMeshType>::Laplacian(poly_mesh,true,1,Damp);
//        else
//            vcg::PolygonalAlgorithm<PolyMeshType>::SmoothPCA(poly_mesh,1,Damp,true,true,0.2,false);



//        //move in tangent space
//        for (size_t i=0;i<poly_mesh.vert.size();i++)
//        {
//            CoordType Pos1=poly_mesh.vert[i].P();
//            CoordType MovDir=Pos1-Pos0[i];
//            ProjectMovingDirection<TriangleMeshType>(TriGrid,tri_mesh,Pos0[i],MovDir);
//            poly_mesh.vert[i].P()=Pos0[i]+MovDir;

//            //if (VertTypes[i]==Corner)continue;
//            CoordType closestPt;
//            if (VertBasis[i]>=0)
//            {
//                closestPt=ClosestPointSegSet<PolyMeshType>(poly_mesh.vert[i].P(),ProjBasis[VertBasis[i]]);
//                poly_mesh.vert[i].P()=closestPt;
//            }
//            else
//            {
//                TriScalarType MinD;
//                TriScalarType MaxD=poly_mesh.bbox.Diag();
//                vcg::tri::GetClosestFaceBase(tri_mesh,TriGrid,poly_mesh.vert[i].P(),MaxD,MinD,closestPt);
//                poly_mesh.vert[i].P()=closestPt;
//            }
//            //          }
//        }

//    }
//}


//template <class PolyMeshType,class TriangleMeshType>
//void SmoothSubdivide(TriangleMeshType &tri_mesh,
//                     PolyMeshType &poly_mesh,
//                     const std::vector<std::pair<size_t,size_t> > &features,
//                     const std::vector<size_t> &featuresC,
//                     const std::vector<size_t> &tri_face_partition,
//                     const std::vector<size_t > &quad_corner,
//                     const std::vector<size_t> &quad_face_partition,
//                     SmoothType SType,
//                     size_t Steps,
//                     typename PolyMeshType::ScalarType Damp,
//                     typename PolyMeshType::ScalarType AvEdge)
//{
//    typedef typename TriangleMeshType::FaceType TriFaceType;
//    typedef typename TriangleMeshType::ScalarType TriScalarType;
//    typedef typename TriangleMeshType::CoordType CoordType;

//    vcg::tri::FieldSmoother<TriangleMeshType>::InitByCurvature(tri_mesh,5);
//    vcg::tri::CrossField<TriangleMeshType>::SetFaceCrossVectorFromVert(tri_mesh);
//    std::vector<TriScalarType> KVal;
//    for (size_t i=0;i<tri_mesh.vert.size();i++)
//    {
//        tri_mesh.vert[i].K1()/=fabs(tri_mesh.vert[i].K1());
//        tri_mesh.vert[i].K2()/=fabs(tri_mesh.vert[i].K2());
//    }

//    for (size_t i=0;i<tri_mesh.vert.size();i++)
//    {
//        KVal.push_back(fabs(tri_mesh.vert[i].K1()));
//        KVal.push_back(fabs(tri_mesh.vert[i].K2()));
//        //std::cout<<"PD1: "<<tri_mesh.vert[i].K1()<<"PD2: "<<tri_mesh.vert[i].K2()<<std::endl;
//    }
//    std::sort(KVal.begin(),KVal.end());
//    int Index=KVal.size()*0.8;
//    TriScalarType MaxV=KVal[Index];
//    for (size_t i=0;i<tri_mesh.vert.size();i++)
//    {
//        tri_mesh.vert[i].K1()/=MaxV;
//        tri_mesh.vert[i].K2()/=MaxV;
//        if (tri_mesh.vert[i].K1()>1)tri_mesh.vert[i].K1()=1;
//        if (tri_mesh.vert[i].K2()>1)tri_mesh.vert[i].K2()=1;
//        tri_mesh.vert[i].K1()=1-tri_mesh.vert[i].K1();
//        tri_mesh.vert[i].K2()=1-tri_mesh.vert[i].K2();
//    }

//    SmoothWithFeatures(tri_mesh,poly_mesh,features,featuresC,tri_face_partition,
//                       quad_corner,quad_face_partition,SType,Steps,Damp,AvEdge);
//}


template <class PolyMeshType,class TriangleMeshType>
void SmoothSubdivide(TriangleMeshType &tri_mesh,
                     PolyMeshType &poly_mesh,
                     const std::vector<std::pair<size_t,size_t> > &trimesh_features,
                     const std::vector<size_t> &trimesh_corners,
                     const std::vector<size_t> &tri_face_partition,
                     const std::vector<size_t > &quad_corner,
                     const std::vector<size_t> &quad_face_partition,
                     //SmoothType SType,
                     size_t Steps,
                     typename PolyMeshType::ScalarType Damp,
                     typename PolyMeshType::ScalarType AvEdge)
{
    typedef typename TriangleMeshType::FaceType TriFaceType;
    typedef typename TriangleMeshType::ScalarType TriScalarType;
    typedef typename TriangleMeshType::CoordType CoordType;

    //    std::vector<CoordType> TargetPos;
    //    BackProjectStepPositions(poly_mesh,tri_mesh,TargetPos);
    //    for (size_t i=0;i<poly_mesh.vert.size();i++)
    //        poly_mesh.vert[i].P()=TargetPos[i];
    //BackProjectOnMesh(poly_mesh,tri_mesh);
    MultiCostraintSmooth(poly_mesh,tri_mesh,trimesh_features,
                         trimesh_corners,tri_face_partition,
                         quad_corner,quad_face_partition,Damp,AvEdge,Steps,2);

    //    return;

    //    //make a copy of the mesh
    //    PolyMeshType poly_swap;
    //    vcg::tri::Append<PolyMeshType,PolyMeshType>::Mesh(poly_swap,poly_mesh);

    //    //subdivide
    //    //vcg::PolygonalAlgorithm<PolyMeshType>::SubdivideStep(poly_swap);
    //    //vcg::PolygonalAlgorithm<PolyMeshType>::SubdivideStep(poly_swap);

    //    //then triangulate
    //    TriangleMeshType poly_tris;
    //    vcg::PolygonalAlgorithm<PolyMeshType>::TriangulateToTriMesh(poly_swap,poly_tris);

    //    //set a grid
    //    vcg::GridStaticPtr<TriFaceType,TriScalarType> TriGrid;
    //    TriGrid.Set(tri_mesh.face.begin(),tri_mesh.face.end());

    //    for (size_t i=0;i<Steps;i++)
    //    {
    //        //smooth
    //        LaplacianEdge(poly_swap,false,Damp);
    //        //the reproject
    //        for (size_t j=0;j<poly_swap.vert.size();j++)
    //        {
    //            TriScalarType MinD;
    //            TriScalarType MaxD=poly_mesh.bbox.Diag();
    //            CoordType closestPt;
    //            vcg::tri::GetClosestFaceBase(tri_mesh,TriGrid,poly_swap.vert[j].P(),MaxD,MinD,closestPt);
    //            poly_swap.vert[j].P()=closestPt;
    //        }
    //    }

    //    for (size_t i=0;i<poly_mesh.vert.size();i++)
    //        poly_mesh.vert[i].P()=poly_swap.vert[i].P();

    //    //    SmoothWithFeatures(tri_mesh,poly_mesh,features,featuresC,tri_face_partition,
    //    //                       quad_corner,quad_face_partition,SType,Steps,Damp,AvEdge);
}

#endif // SMOOTH_MESH_H
