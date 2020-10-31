#pragma once
#include <memory>
#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <OpenMesh/Core/Mesh/PolyConnectivity.hh>
#include "../../container_cast.h"
#include "geodesic/geodesic_algorithm_dijkstra.h"
#include "geodesic/geodesic_algorithm_exact.h"
#include "geodesic/geodesic_algorithm_subdivision.h"

namespace kt84 {

    struct GeodesicPoint {
        OpenMesh::PolyConnectivity::VHandle vhandle;
        OpenMesh::PolyConnectivity::EHandle ehandle;
        OpenMesh::PolyConnectivity::FHandle fhandle;
        Eigen::Vector2d ecoord;         // barycentric coord for edge
        Eigen::Vector3d fcoord;         // barycentric coord for face
        GeodesicPoint(){}
        GeodesicPoint(OpenMesh::PolyConnectivity::VHandle vhandle) : vhandle(vhandle) {}
        GeodesicPoint(OpenMesh::PolyConnectivity::EHandle ehandle, const Eigen::Vector2d& ecoord) : ehandle(ehandle), ecoord(ecoord) {}
        GeodesicPoint(OpenMesh::PolyConnectivity::FHandle fhandle, const Eigen::Vector3d& fcoord) : fhandle(fhandle), fcoord(fcoord) {}
        GeodesicPoint(const geodesic::SurfacePoint& p) {
            switch (p.type()) {
            case geodesic::VERTEX:
                vhandle = OpenMesh::PolyConnectivity::VHandle(p.base_element()->id);
                break;
            case geodesic::EDGE:
                ehandle = OpenMesh::PolyConnectivity::EHandle(p.base_element()->id);
                break;
            }
        }
        operator geodesic::SurfacePoint() const {
            geodesic::SurfacePoint p;
            //...
            return p;
        }
        bool is_vertex() const { return vhandle.is_valid(); }
        bool is_edge  () const { return ehandle.is_valid(); }
        bool is_face  () const { return fhandle.is_valid(); }
    };
    
    template <class TMeshBase, class TMesh>
    struct Geodesic : public DerivedPtrHolder<TMesh, Geodesic <TMeshBase, TMesh>> {
        void geodesic_init() {
            auto mesh = get_mesh();
            std::vector<double> points;
            points.reserve(mesh->n_vertices() * 3);
            for (auto v : mesh->vertices()) {
                auto p = mesh->point(v);
                points.push_back(p[0]);
                points.push_back(p[1]);
                points.push_back(p[2]);
            }
            std::vector<unsigned> faces;
            faces.reserve(mesh->n_faces() * 3);
            for (auto f : mesh->faces()) {
                if (mesh->valence(f) != 3) {
                    std::cerr << "Error: geodesic algorithm does not work on polygonal meshes!\n";
                    assert(false);
                }
                for (auto v : mesh->fv_range(f))
                    faces.push_back(v.idx());
            }
            geodesic.mesh.initialize_mesh_data(points, faces);
            geodesic_set_algorithm_exact();      // default setting
        }
        void geodesic_set_algorithm_exact() {
            geodesic.algorithm = std::make_shared<geodesic::GeodesicAlgorithmExact>(&geodesic.mesh);
        }
        void geodesic_set_algorithm_dijkstra() {
            geodesic.algorithm = std::make_shared<geodesic::GeodesicAlgorithmDijkstra>(&geodesic.mesh);
        }
        void geodesic_set_algorithm_subdivision(int subdivision_level) {
            geodesic.algorithm = std::make_shared<geodesic::GeodesicAlgorithmSubdivision>(&geodesic.mesh, subdivision_level);
        }
        std::vector<GeodesicPoint> geodesic_compute(const GeodesicPoint& source, const GeodesicPoint& target) {
            assert(geodesic.algorithm);
            std::vector<geodesic::SurfacePoint> path;
            geodesic::SurfacePoint source_(source), target_(target);
            geodesic.algorithm->geodesic(source_, target_, path);
            return container_cast<GeodesicPoint>(path);
        }
        struct Data {
            geodesic::Mesh mesh;
            std::shared_ptr<geodesic::GeodesicAlgorithmBase> algorithm;
            Data() {}
            Data(const Data& src)
                : mesh(src.mesh)
            {
                switch (src.algorithm->type()) {
                case geodesic::GeodesicAlgorithmBase::EXACT:
                    algorithm = std::make_shared<geodesic::GeodesicAlgorithmExact>(&mesh);
                    break;
                case geodesic::GeodesicAlgorithmBase::DIJKSTRA:
                    algorithm = std::make_shared<geodesic::GeodesicAlgorithmDijkstra>(&mesh);
                    break;
                case geodesic::GeodesicAlgorithmBase::SUBDIVISION:
                {
                    auto subdivision_level = static_cast<geodesic::GeodesicAlgorithmSubdivision*>(algorithm.get())->subdivision_level();
                    algorithm = std::make_shared<geodesic::GeodesicAlgorithmSubdivision>(&mesh, subdivision_level);
                    break;
                }
                default:
                    assert(false);
                }
            }
            Data& Data::operator=(const Data& src) {
                mesh = src.mesh;
                switch (src.algorithm->type()) {
                case geodesic::GeodesicAlgorithmBase::EXACT:
                    algorithm = std::make_shared<geodesic::GeodesicAlgorithmExact>(&mesh);
                    break;
                case geodesic::GeodesicAlgorithmBase::DIJKSTRA:
                    algorithm = std::make_shared<geodesic::GeodesicAlgorithmDijkstra>(&mesh);
                    break;
                case geodesic::GeodesicAlgorithmBase::SUBDIVISION:
                {
                    auto subdivision_level = static_cast<geodesic::GeodesicAlgorithmSubdivision*>(algorithm.get())->subdivision_level();
                    algorithm = std::make_shared<geodesic::GeodesicAlgorithmSubdivision>(&mesh, subdivision_level);
                    break;
                }
                default:
                    assert(false);
                }
                return *this;
            }
        } geodesic;
    private:
        const TMesh* get_mesh() const { return DerivedPtrHolder<TMesh, Geodesic<TMeshBase, TMesh>>::derived_ptr; }
    };
}
