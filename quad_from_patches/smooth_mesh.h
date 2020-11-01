#ifndef SMOOTH_MESH_H
#define SMOOTH_MESH_H

#include <vector>
#include <vcg/complex/algorithms/polygonal_algorithms.h>
#include <wrap/io_trimesh/export.h>
#include <vcg/complex/algorithms/implicit_smooth.h>

//#include "field_smoother.h"

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
//    typedef typename TriangleMeshType::FaceType TriFaceType;
//    typedef typename TriangleMeshType::ScalarType TriScalarType;
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


template <class PolyMeshType,class TriangleMeshType>
void GetCorners(//const std::vector<std::pair<size_t,size_t> > &FeatureEdges,
                const std::vector<size_t > &tri_feature_C,
                const TriangleMeshType &tri_mesh,
                const std::vector<size_t > &quad_corner,
                const PolyMeshType &poly_mesh,
                //const typename PolyMeshType::ScalarType &AvEdge,
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
        //if (MinD>AvEdge)continue;
        cornersIdx.push_back(CornerI);
    }
}


enum ProjType{ProjNone,ProjSuface,ProjSharp,ProjCorner};
struct ProjectionBase
{
    std::vector<ProjType> VertProjType;
    //edge sequences as pair of vertices
    std::vector<std::vector<std::pair<size_t,size_t> > > SharpEdge;
    //indexes of vertices
    //std::vector<size_t> Corner;
    void Clear()
    {VertProjType.clear();SharpEdge.clear();}//Corner.clear();}
};

enum SmoothType{Laplacian,TemplateFit};

template <class PolyMeshType>
void LaplacianEdgePos(const PolyMeshType &poly_mesh,
                      const ProjectionBase &PolyProjBase,
                      const typename PolyMeshType::ScalarType Damp,
                      std::vector<typename PolyMeshType::CoordType> &TargetPos,
                      const std::vector<bool> &BlockedF,
                      const std::vector<bool> &BlockedV)
{
    //typedef typename PolyMeshType::ScalarType ScalarType;
    typedef typename PolyMeshType::CoordType CoordType;
    typedef typename PolyMeshType::FaceType FaceType;

    TargetPos = std::vector<CoordType>(poly_mesh.vert.size(),CoordType(0,0,0));
    std::vector<size_t> NumPos(poly_mesh.vert.size(),0);
    for (size_t i=0;i<poly_mesh.face.size();i++)
    {
        int sizeP=poly_mesh.face[i].VN();
        for (size_t j=0;j<poly_mesh.face[i].VN();j++)
        {
            FaceType *oppF=poly_mesh.face[i].FFp(j);
            size_t oppI=vcg::tri::Index(poly_mesh,oppF);
            bool AddContribute=true;
            //only once for edge
            //if (oppF<&poly_mesh.face[i])AddContribute=false;
            if (BlockedF[i])AddContribute=false;
            if (BlockedF[oppI])AddContribute=false;
            if (!AddContribute)continue;
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
        if ((PolyProjBase.VertProjType[i]==ProjCorner)||(BlockedV[i]))
            //if (BlockedV[i])
            TargetPos[i]=poly_mesh.vert[i].cP();
        else
            TargetPos[i]=poly_mesh.vert[i].cP()*Damp+TargetPos[i]*(1-Damp);
    }
}

template <class PolyMeshType>
void LaplacianPos(const PolyMeshType &poly_mesh,
                  const typename PolyMeshType::ScalarType Damp,
                  std::vector<typename PolyMeshType::CoordType> &TargetPos,
                  bool only_sel_contribute)
{
    //typedef typename PolyMeshType::ScalarType ScalarType;
    typedef typename PolyMeshType::CoordType CoordType;
    typedef typename PolyMeshType::FaceType FaceType;

    TargetPos = std::vector<CoordType>(poly_mesh.vert.size(),CoordType(0,0,0));
    std::vector<size_t> NumPos(poly_mesh.vert.size(),0);
    for (size_t i=0;i<poly_mesh.face.size();i++)
    {
        int sizeP=poly_mesh.face[i].VN();
        for (int j=0;j<poly_mesh.face[i].VN();j++)
        {
            FaceType *oppF=poly_mesh.face[i].FFp(j);
            //size_t oppI=vcg::tri::Index(poly_mesh,oppF);
            bool AddContribute=true;
            //only once for edge
            if (oppF<&poly_mesh.face[i])AddContribute=false;
            //if (BlockedF[i])AddContribute=false;
            //if (BlockedF[oppI])AddContribute=false;
            if (!AddContribute)continue;
            CoordType Pos0=poly_mesh.face[i].cP(j);
            CoordType Pos1=poly_mesh.face[i].cP((j+1)%sizeP);
            size_t VIndex0=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V(j));
            size_t VIndex1=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V((j+1)%sizeP));
            bool IsS0=poly_mesh.vert[VIndex0].IsS();
            bool IsS1=poly_mesh.vert[VIndex1].IsS();

            if ((!only_sel_contribute)||((only_sel_contribute)&&(IsS1)))
            {
                NumPos[VIndex0]++;
                TargetPos[VIndex0]+=Pos1;
            }

            if ((!only_sel_contribute)||((only_sel_contribute)&&(IsS0)))
            {
                NumPos[VIndex1]++;
                TargetPos[VIndex1]+=Pos0;
            }
        }
    }

    for (size_t i=0;i<TargetPos.size();i++)
    {
        if (NumPos[i]<=1)
            TargetPos[i]=poly_mesh.vert[i].cP();
        else
        {
            CoordType AvgPos=TargetPos[i]/NumPos[i];
            TargetPos[i]=poly_mesh.vert[i].cP()*Damp + AvgPos*(1-Damp);
        }
    }
}

template <class PolyMeshType>
void TemplatePos(PolyMeshType &poly_mesh,
                 const ProjectionBase &PolyProjBase,
                 const typename PolyMeshType::ScalarType Damp,
                 std::vector<typename PolyMeshType::CoordType> &TargetPos)
{
    //typedef typename PolyMeshType::ScalarType ScalarType;
    typedef typename PolyMeshType::CoordType CoordType;

    TargetPos.clear();

    //save old position
    vcg::tri::UpdateFlags<PolyMeshType>::VertexClearS(poly_mesh);
    std::vector<CoordType> OldPos;
    for (size_t i=0;i<poly_mesh.vert.size();i++)
    {
        if ((PolyProjBase.VertProjType[i]==ProjCorner)||
                (PolyProjBase.VertProjType[i]==ProjSharp))
            poly_mesh.vert[i].SetS();

        OldPos.push_back(poly_mesh.vert[i].P());
    }

    //smooth
    vcg::PolygonalAlgorithm<PolyMeshType>::SmoothPCA(poly_mesh,1,Damp,true,true,0.3,true,true);

    vcg::tri::UpdateFlags<PolyMeshType>::VertexClearS(poly_mesh);

    for (size_t i=0;i<poly_mesh.vert.size();i++)
    {
        TargetPos.push_back(poly_mesh.vert[i].P());
        poly_mesh.vert[i].P()=OldPos[i];
    }
}


template <class PolyMeshType,class TriMeshType>
void GetQuadInterpW(const typename TriMeshType::FaceType &FTris,
                    const typename PolyMeshType::FaceType &FPoly,
                    const typename TriMeshType::CoordType &ClosestPt,
                    std::vector<typename TriMeshType::ScalarType> &VertWeigths)
{
    typedef typename PolyMeshType::ScalarType ScalarType;
    typedef typename PolyMeshType::CoordType CoordType;
    //typedef typename TriMeshType::FaceType TriFaceType;

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
                          std::vector<typename PolyMeshType::CoordType> &TargetPos,
                          const bool ProjEdgeM)
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
            if (ProjEdgeM)
            {
                size_t IndexE;
                ScalarType t,MinD;
                CoordType Clos;
                ClosestPointEMesh(TestPos,EdgeM,IndexE,t,MinD,Clos);
                TargetPos.push_back(Clos);
            }
            else
            {
                int IndexSh=SharpProj[i];
                assert(IndexSh>=0);
                assert(IndexSh<TriProjBase.SharpEdge.size());
                size_t ClosestSeg;
                CoordType ClosestPos;
                ClosestPointEdgeSet(TestPos,TriM,TriProjBase.SharpEdge[IndexSh],ClosestSeg,ClosestPos);
                TargetPos.push_back(ClosestPos);
            }
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


template <class PolyMeshType,class TriMeshType>
void GetProjectionBasis(const BasicMesh &edge_mesh,
                         TriMeshType &tri_mesh,
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
    typedef typename TriMeshType::ScalarType ScalarType;
    typedef typename TriMeshType::CoordType CoordType;

    PBaseTris.Clear();
    PBasePoly.Clear();

    PBaseTris.VertProjType.resize(tri_mesh.vert.size(),ProjSuface);
    PBasePoly.VertProjType.resize(poly_mesh.vert.size(),ProjSuface);


    //get for each corner of the tri mesh the corresponding on the poly mesh
    //PBaseTris.Corner=FeatureTrisC;
    std::vector<size_t> FeaturePolyC;
    GetCorners(FeatureTrisC,tri_mesh,poly_corner,poly_mesh,FeaturePolyC);
    //GetCorners(FeatureTris,FeatureTrisC,tri_mesh,poly_corner,poly_mesh,AvEdge,FeaturePolyC);


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
        size_t IndexFOpp;
        if (!vcg::face::IsBorder(tri_mesh.face[IndexF],IndexE))
        {
            IndexFOpp=vcg::tri::Index(tri_mesh,tri_mesh.face[IndexF].cFFp(IndexE));
            Part1=tri_face_partition[IndexFOpp];
        }

        //see if there is such configuration of adjacent patches
        if (Part0==Part1)
        {
            TriMeshType testM;
            typename TriMeshType::FaceType *f0,*f1;
            f0=&tri_mesh.face[IndexF];
            f1=&tri_mesh.face[IndexFOpp];
            vcg::tri::Allocator<TriMeshType>::AddFace(testM,f0->P(0),f0->P(1),f0->P(2));
            vcg::tri::Allocator<TriMeshType>::AddFace(testM,f1->P(0),f1->P(1),f1->P(2));
            vcg::tri::io::ExporterPLY<TriMeshType>::Save(testM,"test_problem.ply");
            assert(0);
        }
        assert(Part0!=Part1);
        PatchPairKey keyPatch(std::min(Part0,Part1),std::max(Part0,Part1));

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
            PBaseTris.SharpEdge.resize(PBaseTris.SharpEdge.size()+1);
            PBaseTris.SharpEdge.back().push_back(EdgeKey);
            BasisMap[keyPatch]=PBaseTris.SharpEdge.size()-1;
        }
        else
        {
            int IndexSharp=BasisMap[keyPatch];
            PBaseTris.SharpEdge[IndexSharp].push_back(EdgeKey);
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
        for (int j=0;j<poly_mesh.face[i].VN();j++)
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

                ScalarType MinD0,MinD1,t;
                CoordType Clos;
                size_t IndexE;
                ClosestPointEMesh(poly_mesh.vert[IndexV0].cP(),edge_mesh,IndexE,t,MinD0,Clos);
                ClosestPointEMesh(poly_mesh.vert[IndexV1].cP(),edge_mesh,IndexE,t,MinD1,Clos);

                size_t IndexSharp=BasisMap[keyPatch];
                assert(IndexSharp<PBasePoly.SharpEdge.size());
                assert(IndexSharp<PBaseTris.SharpEdge.size());

                bool closeI0=(MinD0<AvEdge/10);
                bool closeI1=(MinD1<AvEdge/10);

                if (closeI0 && closeI1)
                    PBasePoly.SharpEdge[IndexSharp].push_back(edgeKey);

                if (closeI0)
                    PBasePoly.VertProjType[IndexV0]=ProjSharp;

                if (closeI1)
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
    //PBaseTris.Corner=FeatureTrisC;
    //GetCorners(FeatureTris,PBaseTris.Corner,tri_mesh,poly_corner,poly_mesh,AvEdge,PBasePoly.Corner);
    for (size_t i=0;i<FeatureTrisC.size();i++)
    {
        size_t IndexC=FeatureTrisC[i];
        PBaseTris.VertProjType[IndexC]=ProjCorner;
    }

    for (size_t i=0;i<FeaturePolyC.size();i++)
    {
        size_t IndexC=FeaturePolyC[i];
        PBasePoly.VertProjType[IndexC]=ProjCorner;
    }
    vcg::tri::UpdateFlags<PolyMeshType>::VertexClearS(poly_mesh);
    for (size_t i=0;i<poly_mesh.vert.size();i++)
    {
        if ((PBasePoly.VertProjType[i]==ProjCorner)||
                (PBasePoly.VertProjType[i]==ProjSharp))
            poly_mesh.vert[i].SetS();
    }

    //    //also add others depending on sharp ones
    //    std::vector<size_t> FeatureCount(poly_mesh.face.size(),0);
    //    for (size_t i=0;i<poly_mesh.face.size();i++)
    //        for (size_t j=0;j<poly_mesh.face[i].VN();j++)
    //        {
    //            size_t IndexV0=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V0(j));
    //            size_t IndexV1=vcg::tri::Index(poly_mesh,poly_mesh.face[i].V0(j));
    //            if ()
    //        }
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
    //assert(PolyProjBase.Corner.size()==PolyProjBase.Corner.size());

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
void SmoothSharpFeatures(PolyMeshType &PolyM,ProjectionBase &PolyProjBase,
                         const BasicMesh &EdgeM,
                         const typename PolyMeshType::ScalarType Damp,
                         std::vector<bool> &BlockedV)
{
    typedef typename PolyMeshType::ScalarType ScalarType;
    typedef typename PolyMeshType::CoordType CoordType;
    typedef typename TriMeshType::FaceType TriFaceType;

    //select only the one on sharp features
    vcg::tri::UpdateFlags<PolyMeshType>::VertexClearS(PolyM);
    for (size_t i=0;i<PolyProjBase.VertProjType.size();i++)
        if ((PolyProjBase.VertProjType[i]==ProjSharp)||(PolyProjBase.VertProjType[i]==ProjCorner))
            PolyM.vert[i].SetS();

    //then do one laplacian step
    std::vector<typename PolyMeshType::CoordType> TargetPos;
    LaplacianPos(PolyM,Damp,TargetPos,true);

    //set value
    for (size_t i=0;i<PolyM.vert.size();i++)
    {
        if (PolyProjBase.VertProjType[i]!=ProjSharp)continue;
        if (BlockedV[i])continue;
        PolyM.vert[i].P()=TargetPos[i];
    }


    for (size_t i=0;i<PolyM.vert.size();i++)
    {
        if (PolyProjBase.VertProjType[i]!=ProjSharp)continue;
        if (BlockedV[i])continue;
        size_t IndexE;
        ScalarType t,MinD;
        CoordType Clos;
        ClosestPointEMesh(PolyM.vert[i].P(),EdgeM,IndexE,t,MinD,Clos);
        //TargetPos.push_back(Clos);
        PolyM.vert[i].P()=Clos;
    }
}

template <class PolyMeshType,class TriMeshType>
void SmoothInternal(PolyMeshType &PolyM,TriMeshType &TriM,
                    vcg::GridStaticPtr<typename TriMeshType::FaceType,
                    typename TriMeshType::ScalarType> &TriGrid,
                    ProjectionBase &TriProjBase,ProjectionBase &PolyProjBase,
                    size_t back_proj_steps,
                    SmoothType SType,
                    const typename PolyMeshType::ScalarType Damp,
                    std::vector<bool> &BlockedV,
                    bool UseSharp)
{
    typedef typename PolyMeshType::CoordType CoordType;
    typedef typename PolyMeshType::ScalarType ScalarType;
    typedef typename TriMeshType::FaceType TriFaceType;

    //std::vector<CoordType> TargetPosProj;
    std::vector<CoordType> TargetPosBackProj;
    std::vector<CoordType> TargetPosSmooth;

    //select only the one on borders
    vcg::tri::UpdateFlags<PolyMeshType>::VertexClearS(PolyM);

    if (!UseSharp)
    {
        vcg::tri::UpdateFlags<PolyMeshType>::VertexSetS(PolyM);
        for (size_t i=0;i<PolyProjBase.VertProjType.size();i++)
        {
            if ((PolyProjBase.VertProjType[i]==ProjSharp)||(PolyProjBase.VertProjType[i]==ProjCorner))
                PolyM.vert[i].ClearS();
        }
        LaplacianPos(PolyM,Damp,TargetPosSmooth,true);
        vcg::tri::UpdateFlags<PolyMeshType>::VertexClearS(PolyM);
    }
    else
    {
        //then do one laplacian step
        //std::vector<typename PolyMeshType::CoordType> TargetPos;
        if (SType==Laplacian)
            LaplacianPos(PolyM,Damp,TargetPosSmooth,false);
        else
            TemplatePos(PolyM,PolyProjBase,Damp,TargetPosSmooth);
    }

    //smooth
    for (size_t i=0;i<PolyM.vert.size();i++)
    {
        if (BlockedV[i])continue;
        if (PolyProjBase.VertProjType[i]==ProjSuface)
            PolyM.vert[i].P()=TargetPosSmooth[i];
    }

    //back projection
    for (size_t i=0;i<back_proj_steps;i++)
    {
        BackProjectStepPositions(PolyM,TriM,TriProjBase,PolyProjBase,TargetPosBackProj);
        for (size_t i=0;i<PolyM.vert.size();i++)
        {
            if (BlockedV[i])continue;
            if (PolyProjBase.VertProjType[i]==ProjSuface)
                PolyM.vert[i].P()=TargetPosBackProj[i];
        }
    }

    for (size_t i=0;i<PolyM.vert.size();i++)
    {
        if (BlockedV[i])continue;
        if (PolyProjBase.VertProjType[i]==ProjSuface)
        {
            CoordType TestPos=PolyM.vert[i].P();
            CoordType closestPt;
            ScalarType MaxD=PolyM.bbox.Diag();
            ScalarType MinD;
            TriFaceType *f=vcg::tri::GetClosestFaceBase(TriM,TriGrid,TestPos,MaxD,MinD,closestPt);
            assert(f!=NULL);
            PolyM.vert[i].P()=closestPt;
        }
    }
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

    BasicMesh EdgeM;
    ExtractEdgeMesh(TriM,features,EdgeM);

    GetProjectionBasis(EdgeM,TriM,features,featuresC,tri_face_partition,PolyM,
                       quad_corner,quad_face_partition,AvEdge,
                       TriProjBase,PolyProjBase);

    std::vector<bool> BlockedV(PolyM.vert.size(),false);

    vcg::GridStaticPtr<TriFaceType,ScalarType> TriGrid;
    TriGrid.Set(TriM.face.begin(),TriM.face.end());

    for (size_t s=0;s<step_num;s++)
    {
        SmoothSharpFeatures<PolyMeshType,TriMeshType>(PolyM,PolyProjBase,EdgeM,Damp,BlockedV);
        SmoothInternal<PolyMeshType,TriMeshType>(PolyM,TriM,TriGrid,
                                                 TriProjBase,PolyProjBase,
                                                 back_proj_steps,
                                                 TemplateFit,Damp,
                                                 BlockedV,true);
    }

}


#endif // SMOOTH_MESH_H
