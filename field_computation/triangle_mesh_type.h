#ifndef MY_TRI_MESH_TYPE
#define MY_TRI_MESH_TYPE

//#define MINDOT -0.99

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/simplex/face/topology.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <wrap/gl/trimesh.h>
//
#include <wrap/igl/smooth_field.h>
#include <wrap/io_trimesh/export_field.h>
#include <iostream>
#include <fstream>
#include <vcg/complex/algorithms/attribute_seam.h>
#include <vcg/complex/algorithms/crease_cut.h>

class PolyFace;
class PolyVertex;

struct PUsedTypes: public vcg::UsedTypes<vcg::Use<PolyVertex>  ::AsVertexType,
        vcg::Use<PolyFace>	::AsFaceType>{};

class PolyVertex:public vcg::Vertex<	PUsedTypes,
        vcg::vertex::Coord3d,
        vcg::vertex::Normal3d//,
        /*vcg::vertex::Mark,
                                                                                                                                        vcg::vertex::BitFlags,
                                                                                                                                        vcg::vertex::Qualityd,
                                                                                                                                        vcg::vertex::TexCoord2d*/>{} ;

class PolyFace:public vcg::Face<
        PUsedTypes
        ,vcg::face::PolyInfo // this is necessary  if you use component in vcg/simplex/face/component_polygon.h
        // It says "this class is a polygon and the memory for its components (e.g. pointer to its vertices
        // will be allocated dynamically")
        ,vcg::face::PFVAdj	 // Pointer to the vertices (just like FVAdj )
        //,vcg::face::PFFAdj	 // Pointer to edge-adjacent face (just like FFAdj )
        ,vcg::face::BitFlags // bit flags
        ,vcg::face::Normal3d // normal
        /*,vcg::face::Color4b  // color
                                                                                                                                        ,vcg::face::Qualityd      // face quality.
                                                                                                                                        ,vcg::face::BitFlags                                                                                                                             ,vcg::face::Mark*/
        ,vcg::face::CurvatureDirf> {
};

class PMesh: public
        vcg::tri::TriMesh<
        std::vector<PolyVertex>,	// the vector of vertices
        std::vector<PolyFace >     // the vector of faces
        >
{
public:

    void TriangulateQuadBySplit(size_t IndexF)
    {

        size_t sizeV=face[IndexF].VN();
        assert(sizeV==4);

        //then reupdate the faces
        VertexType * v0=face[IndexF].V(0);
        VertexType * v1=face[IndexF].V(1);
        VertexType * v2=face[IndexF].V(2);
        VertexType * v3=face[IndexF].V(3);

        face[IndexF].Dealloc();
        face[IndexF].Alloc(3);
        face[IndexF].V(0)=v0;
        face[IndexF].V(1)=v1;
        face[IndexF].V(2)=v3;

        vcg::tri::Allocator<PMesh>::AddFaces(*this,1);

        face.back().Alloc(3);
        face.back().V(0)=v1;
        face.back().V(1)=v2;
        face.back().V(2)=v3;
    }

    void TriangulateQuadBySplit()
    {
        for (size_t i=0;i<face.size();i++)
        {
            if(face[i].VN()!=4)continue;
            TriangulateQuadBySplit(i);
        }
    }
};

class MyTriFace;
class MyTriEdge;
class MyTriVertex;

enum FeatureKind{ETConcave,ETConvex,ETNone};

struct TriUsedTypes: public vcg::UsedTypes<vcg::Use<MyTriVertex>::AsVertexType,
        vcg::Use<MyTriEdge>::AsEdgeType,
        vcg::Use<MyTriFace>::AsFaceType>{};

class MyTriVertex:public vcg::Vertex<TriUsedTypes,
        vcg::vertex::Coord3d,
        vcg::vertex::Normal3d,
        vcg::vertex::VFAdj,
        vcg::vertex::BitFlags,
        vcg::vertex::CurvatureDird,
        vcg::vertex::Qualityd,
        vcg::vertex::Mark>
{

};

class MyTriEdge : public vcg::Edge<
        TriUsedTypes,
        vcg::edge::VertexRef,
        vcg::edge::BitFlags> {};


class MyTriFace:public vcg::Face<TriUsedTypes,
        vcg::face::VertexRef,
        vcg::face::VFAdj,
        vcg::face::FFAdj,
        vcg::face::BitFlags,
        vcg::face::Normal3d,
        vcg::face::CurvatureDird,
        vcg::face::Qualityd,
        vcg::face::WedgeTexCoord2d,
        vcg::face::Mark>
{

public:

    size_t IndexOriginal;

    void ImportData(const MyTriFace  & left )
    {
        IndexOriginal=left.IndexOriginal;
        vcg::Face<TriUsedTypes,
                vcg::face::VertexRef,
                vcg::face::VFAdj,
                vcg::face::FFAdj,
                vcg::face::BitFlags,
                vcg::face::Normal3d,
                vcg::face::CurvatureDird,
                vcg::face::Qualityd,
                vcg::face::WedgeTexCoord2d,
                vcg::face::Mark>::ImportData(left);
    }

    FeatureKind FKind[3];
};

enum GoemPrecondition{NOVertManifold,NOFaceManifold,
                      DegenerateFace,DegenerateVertex,
                      UnreferencedVert,
                      IsOk};


class MyTriMesh: public vcg::tri::TriMesh< std::vector<MyTriVertex>,
        std::vector<MyTriEdge>,
        std::vector<MyTriFace > >
{
    typedef std::pair<CoordType,CoordType> CoordPair;
    std::set< CoordPair > FeaturesCoord;
public:

    ScalarType LimitConcave;

private:

    // Basic subdivision class
    struct SplitLev : public   std::unary_function<vcg::face::Pos<FaceType> ,CoordType >
    {
        std::map<CoordPair,CoordType> *SplitOps;

        void operator()(VertexType &nv,vcg::face::Pos<FaceType>  ep)
        {
            VertexType* v0=ep.f->V0(ep.z);
            VertexType* v1=ep.f->V1(ep.z);

            assert(v0!=v1);

            CoordType Pos0=v0->P();
            CoordType Pos1=v1->P();

            CoordPair CoordK(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
            assert(SplitOps->count(CoordK)>0);
            nv.P()=(*SplitOps)[CoordK];
        }

        vcg::TexCoord2<ScalarType> WedgeInterp(vcg::TexCoord2<ScalarType> &t0,
                                               vcg::TexCoord2<ScalarType> &t1)
        {
            return (vcg::TexCoord2<ScalarType>(0,0));
        }

        SplitLev(std::map<CoordPair,CoordType> *_SplitOps){SplitOps=_SplitOps;}
    };

    class EdgePred
    {

        std::map<CoordPair,CoordType> *SplitOps;

    public:

        bool operator()(vcg::face::Pos<FaceType> ep) const
        {
            VertexType* v0=ep.f->V0(ep.z);
            VertexType* v1=ep.f->V1(ep.z);

            assert(v0!=v1);

            CoordType Pos0=v0->P();
            CoordType Pos1=v1->P();

            CoordPair CoordK(std::min(Pos0,Pos1),std::max(Pos0,Pos1));

            return (SplitOps->count(CoordK)>0);
        }

        EdgePred(std::map<CoordPair,CoordType> *_SplitOps){SplitOps=_SplitOps;}
    };

    void InitEdgeType()
    {
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (IsConcaveEdge(face[i],j))
                    face[i].FKind[j]=ETConcave;
                else
                    face[i].FKind[j]=ETConvex;
            }
    }

    void InitFeatureCoordsTable()
    {
        FeaturesCoord.clear();
        for (size_t i=0;i<face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;
                CoordPair PosEdge(std::min(face[i].P0(j),face[i].P1(j)),
                                  std::max(face[i].P0(j),face[i].P1(j)));
                FeaturesCoord.insert(PosEdge);
            }
        }
    }

    void SetFeatureFromTable()
    {
        for (size_t i=0;i<face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                face[i].ClearFaceEdgeS(j);
                CoordPair PosEdge(std::min(face[i].P0(j),face[i].P1(j)),
                                  std::max(face[i].P0(j),face[i].P1(j)));
                if(FeaturesCoord.count(PosEdge)==0)continue;
                face[i].SetFaceEdgeS(j);
            }
        }
    }

    bool RefineInternalFacesStepFromEdgeSel()
    {
        InitFeatureCoordsTable();
        std::vector<int> to_refine_face;
        for (size_t i=0;i<face.size();i++)
        {
            //find the number of edges
            int Num=0;
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;
                Num++;
            }
            if (Num==3)
                to_refine_face.push_back(i);
        }
        if (to_refine_face.size()==0)return false;

        std::cout<<"Performing "<<to_refine_face.size()<< " face refinement ops"<<std::endl;
        for (size_t j=0;j<to_refine_face.size();j++)
        {
            int IndexF=to_refine_face[j];
            CoordType PD1=face[IndexF].PD1();
            CoordType PD2=face[IndexF].PD2();
            CoordType NewPos=(face[IndexF].P(0)+
                              face[IndexF].P(1)+
                              face[IndexF].P(2))/3;
            vcg::tri::Allocator<MeshType>::AddVertex(*this,NewPos);
            VertexType *V0=face[IndexF].V(0);
            VertexType *V1=face[IndexF].V(1);
            VertexType *V2=face[IndexF].V(2);
            VertexType *V3=&vert.back();
            face[IndexF].V(2)=V3;
            vcg::tri::Allocator<MeshType>::AddFace(*this,V1,V2,V3);
            face.back().PD1()=PD1;
            face.back().PD2()=PD2;
            vcg::tri::Allocator<MeshType>::AddFace(*this,V2,V0,V3);
            face.back().PD1()=PD1;
            face.back().PD2()=PD2;
        }
        UpdateDataStructures();
        SetFeatureFromTable();
        return true;
    }

    bool IsConcaveEdge(const FaceType &f0,int IndexE)
    {
        FaceType *f1=f0.cFFp(IndexE);
        if (f1==&f0)return false;
        CoordType N0=f0.cN();
        CoordType N1=f1->cN();
        CoordType EdgeDir=f0.cP1(IndexE)-f0.cP0(IndexE);
        EdgeDir.Normalize();
        CoordType Cross=N0^N1;
        return ((Cross*EdgeDir)<LimitConcave);
    }

    bool SplitAdjacentTerminalVertices()
    {
        InitFeatureCoordsTable();
        std::vector<size_t> PerVertConcaveEdge(vert.size(),0);
        std::vector<size_t> PerVertConvexEdge(vert.size(),0);
        //count concave vs convex
        for (size_t i=0;i<face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;
                size_t IndexV0=vcg::tri::Index(*this,face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(*this,face[i].V1(j));
                //only on one side
                if (IndexV0>IndexV1)continue;
                if (IsConcaveEdge(face[i],j))
                {
                    PerVertConcaveEdge[IndexV0]++;
                    PerVertConcaveEdge[IndexV1]++;
                }
                else
                {
                    PerVertConvexEdge[IndexV0]++;
                    PerVertConvexEdge[IndexV1]++;
                }
            }
        }

        //count concave vs convex
        std::map<CoordPair,CoordType> ToBeSplitted;
        std::set<CoordPair> NewFeatureEdges;
        for (size_t i=0;i<face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;
                size_t IndexV0=vcg::tri::Index(*this,face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(*this,face[i].V1(j));
                size_t ConcaveEV0=PerVertConcaveEdge[IndexV0];
                size_t ConcaveEV1=PerVertConcaveEdge[IndexV1];
                size_t ConvexEV0=PerVertConvexEdge[IndexV0];
                size_t ConvexEV1=PerVertConvexEdge[IndexV1];
                size_t NumEV0=ConcaveEV0+ConvexEV0;
                size_t NumEV1=ConcaveEV1+ConvexEV1;
                bool IsCornerV0=false;
                bool IsCornerV1=false;

                if (NumEV0==1)IsCornerV0=true;
                if (NumEV0>2)IsCornerV0=true;
                if ((ConcaveEV0>0)&&(ConvexEV0>0))IsCornerV0=true;

                if (NumEV1==1)IsCornerV1=true;
                if (NumEV1>2)IsCornerV1=true;
                if ((ConcaveEV1>0)&&(ConvexEV1>0))IsCornerV1=true;

                //                std::cout<<"ConcaveEV0 "<<ConcaveEV0<<std::endl;
                //                std::cout<<"ConcaveEV1 "<<ConcaveEV1<<std::endl;
                //                std::cout<<"ConvexEV0 "<<ConvexEV0<<std::endl;
                //                std::cout<<"ConvexEV1 "<<ConvexEV1<<std::endl;
                if (IsCornerV0 && IsCornerV1)
                {
                    CoordType P0=face[i].P0(j);
                    CoordType P1=face[i].P1(j);
                    CoordPair Key(std::min(P0,P1),std::max(P0,P1));
                    CoordType Mid=(P0+P1)/2;
                    ToBeSplitted[Key]=Mid;
                    CoordPair newEdge0(std::min(P0,Mid),std::max(P0,Mid));
                    CoordPair newEdge1(std::min(P1,Mid),std::max(P1,Mid));
                    NewFeatureEdges.insert(newEdge0);
                    NewFeatureEdges.insert(newEdge1);
                }
            }
        }


        std::cout<<"Performing "<<ToBeSplitted.size()<< " split ops"<<std::endl;
        if (ToBeSplitted.size()==0)return false;

        SplitLev splMd(&ToBeSplitted);
        EdgePred eP(&ToBeSplitted);

        //do the final split
        bool done=vcg::tri::RefineE<MeshType,SplitLev,EdgePred>(*this,splMd,eP);

        //set old features
        SetFeatureFromTable();

        //and the new ones
        for (size_t i=0;i<face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                CoordType P0=face[i].P0(j);
                CoordType P1=face[i].P1(j);
                CoordPair Key(std::min(P0,P1),std::max(P0,P1));
                if (NewFeatureEdges.count(Key)==0)continue;
                face[i].SetFaceEdgeS(j);
            }
        }
        return done;
    }

    bool SplitAdjacentEdgeSharpFromEdgeSel()
    {
        InitFeatureCoordsTable();
        vcg::tri::UpdateSelection<MeshType>::VertexClear(*this);
        //InitFaceEdgeSelFromFeatureSeq();

        std::set<std::pair<CoordType,CoordType> > EdgePos;

        for (size_t i=0;i<face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;
                int VIndex0=vcg::tri::Index(*this,face[i].V0(j));
                int VIndex1=vcg::tri::Index(*this,face[i].V1(j));
                CoordType P0=vert[VIndex0].P();
                CoordType P1=vert[VIndex1].P();
                vert[VIndex0].SetS();
                vert[VIndex1].SetS();
                EdgePos.insert(std::pair<CoordType,CoordType>(std::min(P0,P1),std::max(P0,P1)));
            }
        }

        //then save the edges to be splitted
        std::map<CoordPair,CoordType> ToBeSplitted;
        for (size_t i=0;i<face.size();i++)
        {
            //find the number of edges
            int Num=0;
            for (size_t j=0;j<3;j++)
            {
                int VIndex0=vcg::tri::Index(*this,face[i].V0(j));
                int VIndex1=vcg::tri::Index(*this,face[i].V1(j));
                if ((!vert[VIndex0].IsS())||(!vert[VIndex1].IsS()))continue;
                CoordType P0=vert[VIndex0].P();
                CoordType P1=vert[VIndex1].P();
                std::pair<CoordType,CoordType> Key(std::min(P0,P1),std::max(P0,P1));
                if (EdgePos.count(Key)==1){Num++;continue;}

                ToBeSplitted[Key]=(P0+P1)/2;
            }
            assert(Num<=2);//this should be already solved
        }
        std::cout<<"Performing "<<ToBeSplitted.size()<< " split ops"<<std::endl;

        SplitLev splMd(&ToBeSplitted);
        EdgePred eP(&ToBeSplitted);

        //do the final split
        bool done=vcg::tri::RefineE<MeshType,SplitLev,EdgePred>(*this,splMd,eP);

        UpdateDataStructures();
        SetFeatureFromTable();
        return done;
    }

    typedef FieldSmoother<MyTriMesh> FieldSmootherType;

public:

    void SmoothField(FieldSmootherType::SmoothParam FieldParam)
    {
        InitFeatureCoordsTable();

        FieldParam.sharp_thr=0.0;
        FieldParam.curv_thr=0.0;
        FieldParam.AddConstr.clear();
        for (size_t i=0;i<face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;
                size_t IndexV0=vcg::tri::Index(*this,face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(*this,face[i].V1(j));
                //only on one side
                if (IndexV0>IndexV1)continue;
                CoordType P0=face[i].V0(j)->P();
                CoordType P1=face[i].V1(j)->P();
                CoordType Dir=P1-P0;
                Dir.Normalize();
                FieldParam.AddConstr.push_back(std::pair<int,CoordType>(i,Dir));

                typename MeshType::FaceType *FOpp=face[i].FFp(j);
                if (FOpp==&face[i])continue;

                int IndexF=vcg::tri::Index(*this,*FOpp);
                FieldParam.AddConstr.push_back(std::pair<int,CoordType>(IndexF,Dir));
            }
        }
        //if there is no alpha curvature and no constraint then set curv_thr
        if ((FieldParam.alpha_curv==0)&&(FieldParam.AddConstr.size()==0))
        {
            FieldParam.curv_thr=0.1;
            FieldParam.align_borders=true;
        }
        std::cout<<"..Alpha.."<<FieldParam.alpha_curv<<std::endl;
        std::cout<<"..Smoothing.."<<std::endl;
        FieldSmootherType::SmoothDirections(*this,FieldParam);
        std::cout<<"..Global Orientation.."<<std::endl;
        vcg::tri::CrossField<MeshType>::OrientDirectionFaceCoherently(*this);
        std::cout<<"..Update Singularities.."<<std::endl;
        vcg::tri::CrossField<MeshType>::UpdateSingularByCross(*this);
        std::cout<<"..Update Features.."<<std::endl;
        UpdateDataStructures();
        SetFeatureFromTable();
    }

    void RefineIfNeeded()
    {

        UpdateDataStructures();
        bool has_refined=false;
        do
        {
            has_refined=false;
            has_refined|=RefineInternalFacesStepFromEdgeSel();
            SolveGeometricIssues();
            has_refined|=SplitAdjacentEdgeSharpFromEdgeSel();
            SolveGeometricIssues();
            //has_refined|=SplitAdjacentTerminalVertices();
            //has_refined|=SplitEdgeSharpSharingVerticesFromEdgeSel();
        }while (has_refined);
        InitEdgeType();
    }

    bool LoadTriMesh(const std::string &filename,bool &allQuad)
    {
        allQuad=false;
        Clear();
        if(filename.empty()) return false;
        int position0=filename.find(".ply");
        int position1=filename.find(".obj");
        int position2=filename.find(".off");

        if (position0!=-1)
        {
            int err=vcg::tri::io::ImporterPLY<MyTriMesh>::Open(*this,filename.c_str());
            if (err!=vcg::ply::E_NOERROR)return false;
            return true;
        }
        if (position1!=-1)
        {
            PMesh pmesh;
            int Mask;
            vcg::tri::io::ImporterOBJ<PMesh>::LoadMask(filename.c_str(), Mask);
            int err=vcg::tri::io::ImporterOBJ<PMesh>::Open(pmesh,filename.c_str(),Mask);
            if ((err!=0)&&(err!=5))return false;
            //check if all quad
            allQuad=true;
            for (size_t i=0;i<pmesh.face.size();i++)
            {
                if (pmesh.face[i].VN()==4)continue;
                allQuad=false;
            }
            if (allQuad)
            {
                for (size_t i=0;i<pmesh.face.size();i++)
                {
                    CoordType PD1[2];
                    PD1[0]=(pmesh.face[i].V(0)->P()-pmesh.face[i].V(1)->P());
                    PD1[1]=(pmesh.face[i].V(3)->P()-pmesh.face[i].V(2)->P());
                    PD1[0].Normalize();
                    PD1[1].Normalize();
                    CoordType PD2[2];
                    PD2[0]=(pmesh.face[i].V(2)->P()-pmesh.face[i].V(1)->P());
                    PD2[1]=(pmesh.face[i].V(3)->P()-pmesh.face[i].V(0)->P());
                    PD2[0].Normalize();
                    PD2[1].Normalize();
                    pmesh.face[i].PD1()=PD1[0]+PD1[1];
                    pmesh.face[i].PD2()=PD2[0]+PD2[1];
                    pmesh.face[i].PD1().Normalize();
                    pmesh.face[i].PD2().Normalize();
                }
                size_t size=pmesh.fn;
                //vcg::PolygonalAlgorithm<PMesh>::Triangulate(pmesh,false);

                pmesh.TriangulateQuadBySplit();
                for (size_t i=0;i<size;i++)
                {
                    pmesh.face[i+size].PD1()=pmesh.face[i].PD1();
                    pmesh.face[i+size].PD2()=pmesh.face[i].PD2();
                }
                //then copy the field
                Clear();
                vcg::tri::Allocator<MyTriMesh>::AddVertices(*this,pmesh.vn);
                vcg::tri::Allocator<MyTriMesh>::AddFaces(*this,pmesh.fn);

                for (size_t i=0;i<pmesh.vert.size();i++)
                    vert[i].P()=pmesh.vert[i].P();

                for (size_t i=0;i<pmesh.face.size();i++)
                {
                    size_t IndexV0=vcg::tri::Index(pmesh,pmesh.face[i].V(0));
                    size_t IndexV1=vcg::tri::Index(pmesh,pmesh.face[i].V(1));
                    size_t IndexV2=vcg::tri::Index(pmesh,pmesh.face[i].V(2));
                    face[i].V(0)=&vert[IndexV0];
                    face[i].V(1)=&vert[IndexV1];
                    face[i].V(2)=&vert[IndexV2];
                    face[i].PD1()=pmesh.face[i].PD1();
                    face[i].PD2()=pmesh.face[i].PD2();
                }
                UpdateDataStructures();
                for (size_t i=0;i<face.size();i++)
                {
                    face[i].PD1()-=face[i].N()*(face[i].PD1()*face[i].N());
                    face[i].PD2()-=face[i].N()*(face[i].PD2()*face[i].N());
                    face[i].PD1().Normalize();
                    face[i].PD2().Normalize();
                    CoordType Avg=face[i].PD1()+face[i].PD2();
                    Avg.Normalize();
                    CoordType Avg1=face[i].N()^Avg;
                    Avg1.Normalize();
                    face[i].PD1()=Avg+Avg1;
                    face[i].PD1().Normalize();
                    face[i].PD2()=face[i].N()^face[i].PD1();
                }
                vcg::tri::CrossField<MyTriMesh>::OrientDirectionFaceCoherently(*this);
                vcg::tri::CrossField<MyTriMesh>::UpdateSingularByCross(*this);
                return true;
            }
            else
            {
                int mask;
                vcg::tri::io::ImporterOBJ<MyTriMesh>::LoadMask(filename.c_str(),mask);
                int err=vcg::tri::io::ImporterOBJ<MyTriMesh>::Open(*this,filename.c_str(),mask);

                if ((err!=0)&&(err!=5))return false;
                return true;
            }
        }
        if (position2!=-1)
        {
            int err=vcg::tri::io::ImporterOFF<MyTriMesh>::Open(*this,filename.c_str());
            if (err!=0)return false;
            return true;
        }
        return false;
    }

    bool SaveSharpFeatures(const std::string &filename)
    {
        if(filename.empty()) return false;
        ofstream myfile;
        myfile.open (filename.c_str());
        size_t num=0;
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;
                num++;
            }
        myfile <<num<<std::endl;
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;
                if (face[i].FKind[j]==ETConcave)
                    myfile <<"0,"<< i <<","<<j<<std::endl;
                else
                    myfile <<"1,"<< i <<","<<j<<std::endl;
            }
        myfile.close();
        return true;
    }

    bool LoadSharpFeatures(const std::string &filename)
    {
        std::cout<<"Loading Sharp Features"<<std::endl;
        for (size_t i=0;i<face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                face[i].FKind[j]=ETNone;
                face[i].ClearFaceEdgeS(j);
            }
        }

        FILE *f=fopen(filename.c_str(),"rt");
        if (f==NULL)return false;
        int Num;
        fscanf(f,"%d/n",&Num);

        for (size_t i=0;i<Num;i++)
        {
            int TypeSh,IndexF,IndexE;
            fscanf(f,"%d,%d,%d,/n",&TypeSh,&IndexF,&IndexE);
            assert((TypeSh==0)||(TypeSh==1));
            assert((IndexE>=0)&&(IndexE<3));
            assert((IndexF>=0)&&(IndexF<face.size()));
            if (TypeSh==0)
                face[IndexF].FKind[IndexE]=ETConcave;
            else
                face[IndexF].FKind[IndexE]=ETConvex;

            face[IndexF].SetFaceEdgeS(IndexE);

            if (!vcg::face::IsBorder(face[IndexF],IndexE))
            {
                FaceType *Fopp=face[IndexF].FFp(IndexE);
                int IOpp=face[IndexF].FFi(IndexE);
                Fopp->SetFaceEdgeS(IOpp);
                if (TypeSh==0)
                    Fopp->FKind[IOpp]=ETConcave;
                else
                    Fopp->FKind[IOpp]=ETConvex;
            }
        }
        fclose(f);
        return true;
    }

    bool LoadSharpFeaturesFL(const std::string &filename)
    {
        std::cout<<"Loading Sharp Features FL Format"<<std::endl;
        for (size_t i=0;i<face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                face[i].FKind[j]=ETNone;
                face[i].ClearFaceEdgeS(j);
            }
        }

        FILE *f=fopen(filename.c_str(),"rt");
        if (f==NULL)return false;
        int Num;
        fscanf(f,"%d/n",&Num);

        std::map<std::pair<size_t,size_t> , std::pair<size_t,size_t> > FaceEdgeMap;
        for (size_t i=0;i<face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                std::pair<size_t,size_t> FaceEdge(i,j);
                size_t IndexV0=vcg::tri::Index(*this,face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(*this,face[i].V1(j));
                std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),
                                             std::max(IndexV0,IndexV1));
                FaceEdgeMap[key]=FaceEdge;
            }
        }

        for (size_t i=0;i<Num;i++)
        {
            int SizeSh;
            fscanf(f,"%d/n",&SizeSh);
            std::vector<int> CurrSh(SizeSh,-1);
            for (size_t j=0;j<SizeSh;j++)
                fscanf(f,"%d",&CurrSh[j]);

            for (size_t j=0;j<CurrSh.size()-1;j++)
            {
                size_t IndexV0=CurrSh[j];
                size_t IndexV1=CurrSh[j+1];
                //std::cout<<"Adding Sharp:"<<IndexV0<<";"<<IndexV1<<std::endl;
                std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),
                                             std::max(IndexV0,IndexV1));
                assert(FaceEdgeMap.count(key)>0);
                size_t IndexF=FaceEdgeMap[key].first;
                size_t IndexE=FaceEdgeMap[key].second;
                face[IndexF].SetFaceEdgeS(IndexE);

                if (IsConcaveEdge(face[IndexF],IndexE))
                    face[IndexF].FKind[IndexE]=ETConcave;
                else
                    face[IndexF].FKind[IndexE]=ETConvex;

                FeatureKind  FKind=face[IndexF].FKind[IndexE];

                if (!vcg::face::IsBorder(face[IndexF],IndexE))
                {
                    FaceType *Fopp=face[IndexF].FFp(IndexE);
                    int IOpp=face[IndexF].FFi(IndexE);
                    Fopp->SetFaceEdgeS(IOpp);
                    Fopp->FKind[IOpp]=FKind;
                }
            }
        }
        fclose(f);
        return true;
    }

    bool SaveTriMesh(const std::string &filename)
    {
        if(filename.empty()) return false;
        int position0=filename.find(".ply");
        int position1=filename.find(".obj");
        int position2=filename.find(".off");

        if (position0!=-1)
        {
            int err=vcg::tri::io::ExporterPLY<MyTriMesh>::Save(*this,filename.c_str());
            if (err!=vcg::ply::E_NOERROR)return false;
            return true;
        }
        if (position1!=-1)
        {
            int mask=0;
            int err=vcg::tri::io::ExporterOBJ<MyTriMesh>::Save(*this,filename.c_str(),mask);
            if ((err!=0)&&(err!=5))return false;
            return true;
        }
        if (position2!=-1)
        {
            int err=vcg::tri::io::ExporterOFF<MyTriMesh>::Save(*this,filename.c_str());
            if (err!=0)return false;
            return true;
        }
        return false;
    }

    bool SaveOrigFace(const std::string &filename)
    {
        if(filename.empty()) return false;
        ofstream myfile;
        myfile.open (filename.c_str());

        for (size_t i=0;i<face.size();i++)
            myfile <<face[i].IndexOriginal <<std::endl;

        myfile.close();
        return true;
    }

    bool SaveField(const std::string &filename)
    {
        if(filename.empty()) return false;
        vcg::tri::io::ExporterFIELD<MyTriMesh>::Save4ROSY(*this,filename.c_str());
        return true;
    }

    GoemPrecondition CheckPreconditions()
    {
        int Num=vcg::tri::Clean<MyTriMesh>::CountNonManifoldVertexFF(*this);
        if (Num>0)return NOVertManifold;
        Num=vcg::tri::Clean<MyTriMesh>::CountNonManifoldEdgeFF(*this);
        if (Num>0)return NOVertManifold;
        Num=vcg::tri::Clean<MyTriMesh>::RemoveDegenerateFace(*this);
        if (Num>0)return DegenerateFace;
        Num=vcg::tri::Clean<MyTriMesh>::RemoveDegenerateVertex(*this);
        if (Num>0)return DegenerateVertex;
        Num=vcg::tri::Clean<MyTriMesh>::RemoveUnreferencedVertex(*this);
        if (Num>0)return UnreferencedVert;
        return IsOk;
    }

    //VCG UPDATING STRUCTURES
    void UpdateDataStructures()
    {
        vcg::tri::Clean<MyTriMesh>::RemoveUnreferencedVertex(*this);
        vcg::tri::Allocator<MyTriMesh>::CompactEveryVector(*this);

        vcg::tri::UpdateBounding<MyTriMesh>::Box(*this);
        vcg::tri::UpdateNormal<MyTriMesh>::PerVertexNormalizedPerFace(*this);
        vcg::tri::UpdateNormal<MyTriMesh>::PerFaceNormalized(*this);
        vcg::tri::UpdateTopology<MyTriMesh>::FaceFace(*this);
        vcg::tri::UpdateTopology<MyTriMesh>::VertexFace(*this);
        vcg::tri::UpdateFlags<MyTriMesh>::VertexBorderFromFaceAdj(*this);
    }

    void InitSharpFeatures(ScalarType SharpAngleDegree)
    {
        if (SharpAngleDegree<=0)return;
        vcg::tri::UpdateFlags<MeshType>::FaceEdgeSelCrease(*this,vcg::math::ToRad(SharpAngleDegree));
        InitEdgeType();
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (vcg::face::IsBorder(face[i],j))
                {
                    face[i].SetFaceEdgeS(j);
                    face[i].FKind[j]=ETConvex;
                }
                if (!vcg::face::IsManifold(face[i],j))
                {
                    face[i].SetFaceEdgeS(j);
                    face[i].FKind[j]=ETConvex;
                }
            }
        std::cout<<"There is "<<SharpLenght()<<" sharp lenght"<<std::endl;
    }

    void GLDrawSharpEdges()
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.9999);
        glLineWidth(5);
        glBegin(GL_LINES);
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;

                if (face[i].FKind[j]==ETConcave)
                    vcg::glColor(vcg::Color4b(0,0,255,255));
                else
                    vcg::glColor(vcg::Color4b(255,0,0,255));

                CoordType Pos0=face[i].P0(j);
                CoordType Pos1=face[i].P1(j);
                vcg::glVertex(Pos0);
                vcg::glVertex(Pos1);
            }
        glEnd();
        glPopAttrib();
    }

    ScalarType SharpLenght()
    {
        ScalarType LSharp=0;
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (face[i].IsFaceEdgeS(j))
                    LSharp+=(face[i].P0(j)-face[i].P1(j)).Norm();
            }
        return (LSharp/2);
    }

    ScalarType Area()
    {
        ScalarType CurrA=0;
        for (size_t i=0;i<face.size();i++)
            CurrA+=vcg::DoubleArea(face[i]);
        return (CurrA/2);
    }

    void SetFeatureValence()
    {
        vcg::tri::UpdateQuality<MyTriMesh>::VertexConstant(*this,0);
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;
                face[i].V0(j)->Q()+=1;
                face[i].V1(j)->Q()+=1;
            }

    }

    void ErodeFeaturesStep()
    {
        SetFeatureValence();

        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;
                if (vcg::face::IsBorder(face[i],j))continue;
                ScalarType Len=(face[i].P0(j)-face[i].P1(j)).Norm();
                if (Len>bbox.Diag()*0.05)continue;

                if ((face[i].V0(j)->Q()==2)||(face[i].V1(j)->Q()==2))
                    face[i].ClearFaceEdgeS(j);
            }
    }

    void DilateFeaturesStep(std::vector<std::pair<size_t,size_t> > &OrigFeatures)
    {
        SetFeatureValence();

        for (size_t i=0;i<OrigFeatures.size();i++)
        {
            size_t IndexF=OrigFeatures[i].first;
            size_t IndexE=OrigFeatures[i].second;
            if ((face[IndexF].V0(IndexE)->Q()==2)&&
                    (!face[IndexF].V0(IndexE)->IsS()))
                face[IndexF].SetFaceEdgeS(IndexE);

            if ((face[IndexF].V1(IndexE)->Q()==2)&&
                    (!face[IndexF].V1(IndexE)->IsS()))
                face[IndexF].SetFaceEdgeS(IndexE);
        }
    }

    void PrintSharpInfo()
    {
        vcg::tri::UpdateQuality<MeshType>::VertexConstant(*this,0);
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;
                if (!face[i].V0(j)->IsB())
                    face[i].V0(j)->Q()+=1;
                if (!face[i].V1(j)->IsB())
                    face[i].V1(j)->Q()+=1;
            }
        size_t Valence1=0;
        for (size_t i=0;i<vert.size();i++)
            if (vert[i].Q()==2)Valence1++;
        std::cout<<"EndPoints "<<Valence1<<std::endl;
    }

    void ErodeDilate(size_t StepNum)
    {
        vcg::tri::UpdateFlags<MyTriMesh>::VertexClearS(*this);
        SetFeatureValence();
        for (size_t i=0;i<vert.size();i++)
            if ((vert[i].Q()>4)||((vert[i].IsB())&&(vert[i].Q()>2)))
                vert[i].SetS();

        std::vector<std::pair<size_t,size_t> > OrigFeatures;

        //save the features
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;
                OrigFeatures.push_back(std::pair<size_t,size_t>(i,j));
            }

        for (size_t s=0;s<StepNum;s++)
            ErodeFeaturesStep();
        for (size_t s=0;s<StepNum;s++)
            DilateFeaturesStep(OrigFeatures);

        PrintSharpInfo();
    }


    inline static void ExtractVertex(const MyTriMesh & srcMesh,
                                     const FaceType & f,
                                     int whichWedge,
                                     const MyTriMesh & dstMesh,
                                     VertexType & v)
    {
        (void)srcMesh;
        (void)dstMesh;

        v.P() = f.cP(whichWedge);
        //v.T() = f.cWT(whichWedge);
    }

    std::set<std::pair<CoordType,CoordType> > ToSplit;

    //    void DilateToSplitStep()
    //    {
    //        vcg::tri::UpdateSelection<MyTriMesh>::VertexClear(*this);
    //        vcg::tri::UpdateQuality<MyTriMesh>::VertexConstant(*this,0);
    //        for (size_t i=0;i<face.size();i++)
    //            for (size_t j=0;j<face[i].VN();j++)
    //            {

    //            }
    //    }

    inline static bool CompareVertex(MyTriMesh & srcMesh,
                                     const VertexType & vA,
                                     const VertexType & vB)
    {
        //(void)srcMesh;

        std::pair<CoordType,CoordType> Key(std::min(vA.cP(),vB.cP()),
                                           std::max(vA.cP(),vB.cP()));
        return (srcMesh.ToSplit.count(Key)==0);
        //return (vA.cT() == vB.cT());
    }

    bool IsFold(const FaceType &f,
                const size_t &IndexE,
                ScalarType MinDot=-0.99)
    {
        if (vcg::face::IsBorder(f,IndexE))return false;
        FaceType *fOpp=f.cFFp(IndexE);
        //int IOpp=f.cFFi(IndexE);
        CoordType N0=f.cN();
        CoordType N1=fOpp->N();
        N0.Normalize();
        N1.Normalize();
        if ((N0*N1)>MinDot)return false;
        return true;
    }

    bool SplitFolds(ScalarType MinDot=-0.99)
    {

        //then save the edges to be splitted
        std::map<CoordPair,CoordType> ToBeSplitted;
        for (size_t i=0;i<face.size();i++)
        {
            //find the number of edges
            for (size_t j=0;j<3;j++)
            {
                if (!IsFold(face[i],j,MinDot))continue;
                int VIndex0=vcg::tri::Index(*this,face[i].V0(j));
                int VIndex1=vcg::tri::Index(*this,face[i].V1(j));
                CoordType P0=vert[VIndex0].P();
                CoordType P1=vert[VIndex1].P();
                std::pair<CoordType,CoordType> Key(std::min(P0,P1),std::max(P0,P1));
                ToBeSplitted[Key]=(P0+P1)/2;
            }
        }
        std::cout<<"Performing "<<ToBeSplitted.size()<< " fold split ops"<<std::endl;

        SplitLev splMd(&ToBeSplitted);
        EdgePred eP(&ToBeSplitted);

        //do the final split

        size_t NumV0=vert.size();
        bool done=vcg::tri::RefineE<MeshType,SplitLev,EdgePred>(*this,splMd,eP);
        size_t NumV1=vert.size();
        std::cout<<"Added "<<NumV1-NumV0<<" vertices"<<std::endl;


        UpdateDataStructures();
        SetFeatureFromTable();
        return done;
    }

    void RemoveFolds(ScalarType MinDot=-0.99)
    {

        InitFeatureCoordsTable();

        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
                face[i].ClearFaceEdgeS(j);

        ScalarType AvgEdge=0;
        size_t Num=0;
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<face[i].VN();j++)
            {
                AvgEdge+=(face[i].P0(j)-face[i].P1(j)).Norm();
                Num++;
                if (vcg::face::IsBorder(face[i],j))continue;
                if (!IsFold(face[i],j,MinDot))continue;
                face[i].SetFaceEdgeS(j);
            }
        AvgEdge/=Num;

        size_t NumV0=vert.size();
        vcg::tri::CutMeshAlongSelectedFaceEdges<MeshType>(*this);
        UpdateDataStructures();
        size_t NumV1=vert.size();
        std::cout<<"Added "<<NumV1-NumV0<<" vertices"<<std::endl;
        for (size_t i=NumV0;i<NumV1;i++)
            vert[i].P()+=vert[i].N()*AvgEdge*0.00001;

        SetFeatureFromTable();
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<face[i].VN();j++)
            {
                if (!vcg::face::IsBorder(face[i],j))continue;
                face[i].SetFaceEdgeS(j);
            }
        //        ToSplit.clear();
        //        //int RemovedNum=0;
        //        for (size_t i=0;i<face.size();i++)
        //            for (size_t j=0;j<face[i].VN();j++)
        //            {
        //                if (vcg::face::IsBorder(face[i],j))continue;
        //                FaceType *fOpp=face[i].FFp(j);
        //                int IOpp=face[i].FFi(j);
        //                CoordType N0=face[i].N();
        //                CoordType N1=fOpp->N();
        //                if ((N0*N1)>MinDot)continue;
        ////                if (!IsFold(face[i],j,MinDot))continue;
        //                if (OnlySharp)
        //                {
        //                    //RemovedNum++;
        //                    face[i].ClearFaceEdgeS(j);
        //                    fOpp->ClearFaceEdgeS(IOpp);
        //                }
        //                else
        //                {
        //                    std::pair<CoordType,CoordType> Key(std::min(face[i].P0(j),face[i].P1(j)),
        //                                                       std::max(face[i].P0(j),face[i].P1(j)));
        //                    //RemovedNum++;
        //                    ToSplit.insert(Key);
        //                }
        //            }

        //        if (!OnlySharp)
        //        {
        //            size_t NumV0=vert.size();
        //            vcg::tri::AttributeSeam::SplitVertex(*this, ExtractVertex, CompareVertex);
        //            vcg::tri::Allocator<TriMesh>::CompactEveryVector(*this);
        //            size_t NumV1=vert.size();
        //            std::cout<<"Added "<<NumV1-NumV0<<" vertices"<<std::endl;
        //            UpdateDataStructures();
        //        }

        //        std::cout<<"Removed "<<ToSplit.size()<<" folded sharp features"<<std::endl;
        //        PrintSharpInfo();
    }



    bool SufficientFeatures(ScalarType SharpFactor)
    {
        ScalarType sqrtCurrA=sqrt(Area());
        ScalarType SharpL=SharpLenght();
        std::cout<<"Sqrt Area "<<sqrtCurrA<<std::endl;
        std::cout<<"Sharp Lenght "<<SharpL<<std::endl;
        ScalarType Ratio=SharpL/sqrtCurrA;
        std::cout<<"Ratio "<<Ratio<<std::endl;
        return(Ratio>SharpFactor);
    }

    size_t NumDuplicatedV()
    {
        std::set<CoordType> Pos;
        vcg::tri::UpdateSelection<MeshType>::VertexClear(*this);
        size_t numDupl=0;
        for (size_t i=0;i<vert.size();i++)
        {
            if (Pos.count(vert[i].P())>0)
            {
                vert[i].SetS();
                numDupl++;
            }
            Pos.insert(vert[i].P());
        }
        return numDupl;
    }

    void Perturb(VertexType &v,ScalarType Magnitudo)
    {
        ScalarType eps=std::numeric_limits<ScalarType>::epsilon()*Magnitudo;
        //take a random direction
        size_t granularity=10000;
        int IntX=(rand()%granularity)-granularity/2;
        int IntY=(rand()%granularity)-granularity/2;
        int IntZ=(rand()%granularity)-granularity/2;
        CoordType Dir=CoordType(IntX,IntY,IntZ);
        Dir.Normalize();
        Dir*=eps;
        //std::cout<<Dir.X()<<";"<<Dir.Y()<<";"<<Dir.Z()<<std::endl;
        v.P()+=Dir;
    }

    bool RepositionDuplicatedV()
    {
        size_t NumD=NumDuplicatedV();
        if (NumD==0)return false;
        //int dilate_step=0;
        ScalarType Magnitudo=2;
        do
        {
            std::cout<<"Repositioning "<<NumD<<" duplicated vertices"<<std::endl;


            //            dilate_step++;
            //            for (size_t i=0;i<dilate_step;i++)
            //            {
            //                vcg::tri::UpdateSelection<MeshType>::FaceFromVertexLoose(*this);
            //                vcg::tri::UpdateSelection<MeshType>::VertexFromFaceLoose(*this);
            //            }
            for (size_t i=0;i<vert.size();i++)
                if (vert[i].IsS())Perturb(vert[i],Magnitudo);

            Magnitudo*=2;
            //vcg::tri::Smooth<MeshType>::VertexCoordLaplacian(*this,1,true);
            NumD=NumDuplicatedV();
        }
        while(NumD>0);
        vcg::tri::UpdateBounding<MeshType>::Box(*this);
        vcg::tri::UpdateNormal<MeshType>::PerVertexNormalizedPerFace(*this);
        vcg::tri::UpdateNormal<MeshType>::PerFaceNormalized(*this);
        return true;
    }

    bool RemoveZeroAreaF()
    {
        //        int nonManifV=0;
        //        int degF=0;
        int zeroAFace=0;
        bool modified=false;
        ScalarType Magnitudo=2;
        do{
            modified=false;
            for (size_t i=0;i<face.size();i++)
            {
                if (vcg::DoubleArea(face[i])>0)continue;
                Perturb(*face[i].V(0),Magnitudo);
                Perturb(*face[i].V(1),Magnitudo);
                Perturb(*face[i].V(2),Magnitudo);
                modified=true;
                zeroAFace++;
            }
            Magnitudo*=2;
        }while (modified);
        vcg::tri::Allocator<MeshType>::CompactEveryVector(*this);
        std::cout<<"Adjusted "<<zeroAFace<<" zero area faces"<<std::endl;
        //        std::cout<<"Removed "<<degF<<" degenerate faces"<<std::endl;
        //        std::cout<<"Removed "<<zeroAFace<<" nonManifV "<<std::endl;
        UpdateDataStructures();
        return modified;
    }

    void RemoveSmallComponents(size_t min_size=100)
    {
        std::vector< std::pair<int, typename MeshType::FacePointer> > CCV;
        vcg::tri::Clean<MeshType>::ConnectedComponents((*this), CCV);
        vcg::tri::ConnectedComponentIterator<MeshType> ci;
        //vcg::tri::UpdateFlags<MeshType>::FaceClearS(*this);
        for(unsigned int i=0;i<CCV.size();++i)
        {
            if (CCV[i].first>min_size)continue;
            //std::vector<typename MeshType::FacePointer> FPV;
            for(ci.start((*this),CCV[i].second);!ci.completed();++ci)
                vcg::tri::Allocator<MeshType>::DeleteFace((*this),(*(*ci)));
        }
        vcg::tri::Clean<MeshType>::RemoveUnreferencedVertex(*this);
        vcg::tri::Allocator<MeshType>::CompactEveryVector(*this);
        UpdateDataStructures();
    }

    void SolveGeometricIssues()
    {
        srand(0);
        bool HasRepositioned=false;
        bool HasSolvedZeroF=false;
        do{
            HasRepositioned=RepositionDuplicatedV();
            HasSolvedZeroF=RemoveZeroAreaF();
        }while (HasRepositioned || HasSolvedZeroF);

    }

    void SetSharp(FaceType &f,int IndexE)
    {
        f.SetFaceEdgeS(IndexE);
        if (IsConcaveEdge(f,IndexE))
            f.FKind[IndexE]=ETConcave;
        else
            f.FKind[IndexE]=ETConvex;

        if (vcg::face::IsBorder(f,IndexE))return;

        MyTriMesh::FaceType *fopp=f.FFp(IndexE);
        int eOpp=f.FFi(IndexE);

        fopp->SetFaceEdgeS(eOpp);
        fopp->FKind[eOpp]= f.FKind[IndexE];
    }

    void ClearSharp(FaceType &f,int IndexE)
    {
        f.ClearFaceEdgeS(IndexE);
        f.FKind[IndexE]=ETNone;
        if (vcg::face::IsBorder(f,IndexE))return;

        MyTriMesh::FaceType *fopp=f.FFp(IndexE);
        int eOpp=f.FFi(IndexE);

        fopp->ClearFaceEdgeS(eOpp);
        fopp->FKind[eOpp]= ETNone;
    }

};


#endif
