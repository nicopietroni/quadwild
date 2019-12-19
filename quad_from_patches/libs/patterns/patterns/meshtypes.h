#ifndef MESHTYPES_H
#define MESHTYPES_H
#include <vcg/complex/complex.h>
//For using Face::EdgePlane
#include<vcg/simplex/face/component_ep.h>
#include <vcg/complex/algorithms/update/halfedge_indexed.h>
using namespace vcg;
class CFace;
class CVertex;
class CEdge;
struct MyUsedTypes : public UsedTypes<	Use<CVertex>::AsVertexType,	Use<CEdge>::AsEdgeType,Use<CFace>::AsFaceType>{};

/// compositing wanted proprieties
class CVertex : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f,vcg::vertex::VFAdj, vcg::vertex::Color4b, vertex::TexCoord2f,vcg::vertex::BitFlags,vcg::vertex::VEAdj, vertex::Qualityf,vertex::Mark>{};
class CEdge   : public Edge<MyUsedTypes, edge::VertexRef,edge::BitFlags,edge::EEAdj,edge::EFAdj,edge::VEAdj>{};
class CFace   : public vcg::Face<  MyUsedTypes, vcg::face::VertexRef, vcg::face::FFAdj,vcg::face::VFAdj, vcg::face::Normal3f, vcg::face::EdgePlane, vcg::face::Color4b, vcg::face::BitFlags,face::Mark,face::WedgeTexCoord2f> {};
class CMesh   : public vcg::tri::TriMesh< std::vector<CVertex>, std::vector<CFace>,std::vector<CEdge>> {};


/* Definition of a mesh of polygons that also supports half-edges
*/
class PFace;
class PVertex;
class PHEdge;
class PEdge;

struct PUsedTypes: public vcg::UsedTypes<vcg::Use<PVertex>  ::AsVertexType,
                                            vcg::Use<PEdge>	::AsEdgeType,
                                            //vcg::Use<PHEdge>::AsHEdgeType,
                                            vcg::Use<PFace>	::AsFaceType
                                            >{};

//class DummyEdge: public vcg::Edge<PolyUsedTypes>{};
class PVertex:public vcg::Vertex<	PUsedTypes,
                                        vcg::vertex::Coord3f,
                                        vcg::vertex::Normal3f,
                                        vcg::vertex::Mark,
                                        vcg::vertex::BitFlags,
                                        vcg::vertex::VHAdj,
                                        vcg::vertex::Qualityf,
                                        vcg::vertex::VFAdj>{} ;

class PEdge : public Edge<PUsedTypes>{};
class PHEdge : public HEdge< PUsedTypes, hedge::BitFlags,
    //hedge::HFAdj,		// pointer to the face
    //hedge::HOppAdj,	// pointer to the opposite edge
    //hedge::HVAdj,		// pointer to the vertex
    //hedge::HNextAdj,	// pointer to the next halfedge
     hedge::HEdgeData		// the previous 4 components (just more handy, you can comment this and uncomment the previous four lines)
    ,hedge::HPrevAdj	// pointer to the previous halfedge
>{};

class PFace:public vcg::Face<
     PUsedTypes
    ,vcg::face::PolyInfo // this is necessary  if you use component in vcg/simplex/face/component_polygon.h
                         // It says "this class is a polygon and the memory for its components (e.g. pointer to its vertices
                         // will be allocated dynamically")
    ,vcg::face::PFVAdj	 // Pointer to the vertices (just like FVAdj )
    ,vcg::face::PFFAdj	 // Pointer to edge-adjacent face (just like FFAdj )
    ,vcg::face::PFHAdj	 // Pointer its half -edges  ( you may need this if you use half edges)
    ,vcg::face::BitFlags // bit flags
    ,vcg::face::Normal3f // normal
    ,vcg::face::Color4b
    ,vcg::face::Qualityf
    ,vcg::face::WedgeTexCoord2f
    ,vcg::face::Mark
    ,vcg::face::EdgePlane
> {};

class PMesh: public
    vcg::tri::TriMesh<
    std::vector<PVertex>,	// the vector of vertices
    std::vector<PFace >, 						// the vector of faces
    //std::vector<PHEdge>		,						// the vector of edges
    std::vector<PEdge> 								// the vector of edges
    >{};



#endif // MESHTYPES_H
