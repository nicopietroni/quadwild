#ifndef DEFAULTMESHTYPES_H
#define DEFAULTMESHTYPES_H

#include <vcg/complex/complex.h>

/* ----- Polygon mesh ----- */

class PolyVertex;
class PolyFace;
class PolyEdge;

struct MyPolyTypes : public vcg::UsedTypes<
        vcg::Use<PolyVertex>::AsVertexType,
        vcg::Use<PolyEdge>::AsEdgeType,
        vcg::Use<PolyFace>::AsFaceType>{};

class PolyVertex : public vcg::Vertex<MyPolyTypes,
        vcg::vertex::Coord3d,
        vcg::vertex::Normal3d,
        vcg::vertex::Color4b,
        vcg::vertex::Qualityd,
        vcg::vertex::BitFlags,
        vcg::vertex::VFAdj,
        vcg::vertex::CurvatureDird>{};

class PolyFace : public vcg::Face<
        MyPolyTypes,
        vcg::face::PolyInfo,
        vcg::face::VertexRef,
        vcg::face::Normal3d,
        vcg::face::Color4b,
        vcg::face::Qualityd,
        vcg::face::BitFlags,
        vcg::face::PFVAdj,
        vcg::face::PFFAdj,
        vcg::face::PVFAdj,
        vcg::face::CurvatureDird,
        vcg::face::Mark,
        vcg::face::WedgeTexCoord2d> {};

class PolyEdge : public vcg::Edge<
        MyPolyTypes,
        vcg::edge::VertexRef,
        vcg::edge::BitFlags> {};

class PolyMesh : public vcg::tri::TriMesh<
        std::vector<PolyVertex>,
        std::vector<PolyEdge>,
        std::vector<PolyFace>> {};


/* ----- Triangle mesh ----- */

class TriangleVertex;
class TriangleFace;
struct MyTriangleTypes : public vcg::UsedTypes<
        vcg::Use<TriangleVertex>::AsVertexType,
        vcg::Use<TriangleFace>::AsFaceType>{};

class TriangleVertex : public vcg::Vertex<
        MyTriangleTypes,
        vcg::vertex::Coord3d,
        vcg::vertex::Normal3d,
        vcg::vertex::VFAdj,
        vcg::vertex::Color4b,
        vcg::vertex::Qualityd,
        vcg::vertex::BitFlags,
        vcg::vertex::CurvatureDird>{};

class TriangleFace : public vcg::Face<
        MyTriangleTypes,
        vcg::face::VertexRef,
        vcg::face::Normal3d,
        vcg::face::Color4b,
        vcg::face::Qualityd,
        vcg::face::BitFlags,
        vcg::face::FFAdj,
        vcg::face::VFAdj,
        vcg::face::CurvatureDird,
        vcg::face::Mark,
        vcg::face::WedgeTexCoord2d> {};

class TriangleMesh : public vcg::tri::TriMesh<
        std::vector<TriangleVertex>,
        std::vector<TriangleFace> > {};


#endif // DEFAULTMESHTYPES_H
