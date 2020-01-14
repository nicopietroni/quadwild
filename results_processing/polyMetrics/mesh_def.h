#ifndef MESH_DEF_H
#define MESH_DEF_H


#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/polygon_support.h>
#include <vcg/complex/algorithms/polygonal_algorithms.h>

#include <wrap/io_trimesh/import_obj.h>
#include <wrap/io_trimesh/export_obj.h>


class PolyVertex;
class PolyEdge;
class PolyFace;

struct MyUsedTypes : public vcg::UsedTypes<
        vcg::Use<PolyVertex>::AsVertexType,
        //        vcg::Use<PolyEdge>::AsEdgeType,
        vcg::Use<PolyFace>::AsFaceType
        > {};

class PolyVertex : public vcg::Vertex<
        MyUsedTypes,
        vcg::vertex::Coord3d,
        vcg::vertex::Qualityd,
        vcg::vertex::Normal3d,
        vcg::vertex::Color4b,
        vcg::vertex::Mark,
        vcg::vertex::BitFlags,
        vcg::vertex::VFAdj
        > {};

class PolyFace : public vcg::Face<
        MyUsedTypes,
        vcg::face::Color4b,
        vcg::face::Normal3d,
        vcg::face::VertexRef,
        vcg::face::BitFlags,
        vcg::face::Mark,
        vcg::face::PolyInfo,
        vcg::face::PFVAdj,
        vcg::face::PFFAdj,
        vcg::face::Qualityd
        > {};


class PolyMesh : public vcg::tri::TriMesh<
        std::vector<PolyVertex>,
        std::vector<PolyFace >
        >{};


typedef typename PolyVertex::CoordType CoordType;
typedef typename PolyVertex::ScalarType ScalarType;

typedef typename PolyMesh::VertexType VertexType;
typedef typename PolyMesh::VertexPointer VertexPointer;

typedef typename PolyMesh::FaceType FaceType;
typedef typename PolyMesh::FacePointer FacePointer;


#endif // MESH_DEF_H
