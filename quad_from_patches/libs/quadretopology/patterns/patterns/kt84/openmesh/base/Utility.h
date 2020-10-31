#pragma once
#include <utility>
#include <vector>
#include <OpenMesh/Core/Mesh/PolyConnectivity.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include "../../RangeAdaptor.h"
#include "../DerivedPtrHolder.h"

namespace kt84 {

template <typename TMeshBase, typename TMesh>
struct Utility : public DerivedPtrHolder<TMesh, Utility<TMeshBase, TMesh>> {
    typedef std::pair<typename TMeshBase::VHandle, typename TMeshBase::VHandle> VertexPair;
    typedef std::pair<typename TMeshBase::Point  , typename TMeshBase::Point  > PointPair;
    typedef std::pair<typename TMeshBase::HHandle, typename TMeshBase::HHandle> HalfedgePair;
    
    // NOTE: every public function is prefixed by "util_" to avoid ambiguity
    
    // circulation around vertex/face stored in std::vector
    std::vector<typename TMeshBase::VHandle> util_vv_vec (typename TMeshBase::VHandle handle) const { std::vector<typename TMeshBase::VHandle> result; result.reserve(get_mesh()->valence(handle)); for (auto x : get_mesh()->vv_range (handle)) result.push_back(x); return result; }
    std::vector<typename TMeshBase::HHandle> util_voh_vec(typename TMeshBase::VHandle handle) const { std::vector<typename TMeshBase::HHandle> result; result.reserve(get_mesh()->valence(handle)); for (auto x : get_mesh()->voh_range(handle)) result.push_back(x); return result; }
    std::vector<typename TMeshBase::HHandle> util_vih_vec(typename TMeshBase::VHandle handle) const { std::vector<typename TMeshBase::HHandle> result; result.reserve(get_mesh()->valence(handle)); for (auto x : get_mesh()->vih_range(handle)) result.push_back(x); return result; }
    std::vector<typename TMeshBase::EHandle> util_ve_vec (typename TMeshBase::VHandle handle) const { std::vector<typename TMeshBase::EHandle> result; result.reserve(get_mesh()->valence(handle)); for (auto x : get_mesh()->ve_range (handle)) result.push_back(x); return result; }
    std::vector<typename TMeshBase::FHandle> util_vf_vec (typename TMeshBase::VHandle handle) const { std::vector<typename TMeshBase::FHandle> result; result.reserve(get_mesh()->valence(handle)); for (auto x : get_mesh()->vf_range (handle)) result.push_back(x); return result; }
    std::vector<typename TMeshBase::VHandle> util_fv_vec (typename TMeshBase::FHandle handle) const { std::vector<typename TMeshBase::VHandle> result; result.reserve(get_mesh()->valence(handle)); for (auto x : get_mesh()->fv_range (handle)) result.push_back(x); return result; }
    std::vector<typename TMeshBase::HHandle> util_fh_vec (typename TMeshBase::FHandle handle) const { std::vector<typename TMeshBase::HHandle> result; result.reserve(get_mesh()->valence(handle)); for (auto x : get_mesh()->fh_range (handle)) result.push_back(x); return result; }
    std::vector<typename TMeshBase::EHandle> util_fe_vec (typename TMeshBase::FHandle handle) const { std::vector<typename TMeshBase::EHandle> result; result.reserve(get_mesh()->valence(handle)); for (auto x : get_mesh()->fe_range (handle)) result.push_back(x); return result; }
    std::vector<typename TMeshBase::FHandle> util_ff_vec (typename TMeshBase::FHandle handle) const { std::vector<typename TMeshBase::FHandle> result; result.reserve(get_mesh()->valence(handle)); for (auto x : get_mesh()->ff_range (handle)) result.push_back(x); return result; }
    
    // additional connectivity info
    typename TMeshBase::HHandle util_next_halfedge_n(typename TMeshBase::HHandle h, int n) const {
        auto mesh = get_mesh();
        for (int i = 0; i < n; ++i)
            h = mesh->next_halfedge_handle(h);
        return h;
    }
    typename TMeshBase::VHandle util_opposite_vertex(typename TMeshBase::EHandle e, typename TMeshBase::VHandle v) const {
        auto mesh = get_mesh();
        auto v0 = mesh->to_vertex_handle(mesh->halfedge_handle(e, 0));
        auto v1 = mesh->to_vertex_handle(mesh->halfedge_handle(e, 1));
        return
            v0 == v ? v1 :
            v1 == v ? v0 :
            typename TMeshBase::VHandle();
    }
    typename TMeshBase::VHandle util_opposite_vertex(typename TMeshBase::HHandle h, typename TMeshBase::VHandle v) const {
        auto mesh = get_mesh();
        return opposite_vertex(mesh, mesh->edge_handle(h), v);
    }
    typename TMeshBase::FHandle util_opposite_face(typename TMeshBase::EHandle e, typename TMeshBase::FHandle f) const {
        auto mesh = get_mesh();
        auto f0 = mesh->face_handle(mesh->halfedge_handle(e, 0));
        auto f1 = mesh->face_handle(mesh->halfedge_handle(e, 1));
        return
            f0 == f ? f1 :
            f1 == f ? f0 :
            typename TMeshBase::FHandle();
    }
    typename TMeshBase::HHandle util_halfedge_handle(typename TMeshBase::VHandle v_from, typename TMeshBase::VHandle v_to) const {
        auto mesh = get_mesh();
        for (auto voh : mesh->voh_range(v_from)) {
            if (mesh->to_vertex_handle(voh) == v_to)
                return voh;
        }
        return typename TMeshBase::HHandle();
    }
    typename TMeshBase::EHandle util_edge_handle(typename TMeshBase::VHandle v0, typename TMeshBase::VHandle v1) const {
        return get_mesh()->edge_handle(util_halfedge_handle(v0, v1));
    }
    typename TMeshBase::FHandle util_face_handle(typename TMeshBase::VHandle v0, typename TMeshBase::VHandle v1, typename TMeshBase::VHandle v2) const {
        auto mesh = get_mesh();
        for (auto voh : mesh->voh_range(v0)) {
            if (mesh->to_vertex_handle(voh) != v1)
                continue;
            auto h1 = mesh->next_halfedge_handle(voh);
            auto h2 = mesh->next_halfedge_handle(h1);
            if (mesh->to_vertex_handle(h1) != v2 || mesh->to_vertex_handle(h2) != v0)
                continue;
            return mesh->face_handle(voh);
        }
        
        return typename TMeshBase::FHandle();
    }
    
    // halfedge-related functions
    VertexPair util_halfedge_to_vertex_pair(typename TMeshBase::HHandle h) const { auto mesh = get_mesh(); return std::make_pair(mesh->from_vertex_handle(h), mesh->to_vertex_handle(h)); }
    PointPair  util_halfedge_to_point_pair (typename TMeshBase::HHandle h) const { auto mesh = get_mesh(); auto p = util_halfedge_to_vertex_pair(h); return std::make_pair(mesh->point(p.first), mesh->point(p.second)); }
    typename TMeshBase::Point util_halfedge_to_vector(typename TMeshBase::HHandle h) const {
        auto mesh = get_mesh();
        return mesh->point(mesh->  to_vertex_handle(h)) -
               mesh->point(mesh->from_vertex_handle(h));
    }
    
    // edge-related functions
    typename TMeshBase::Point::value_type util_edge_length(typename TMeshBase::EHandle e) const {
        return util_halfedge_to_vector(get_mesh()->halfedge_handle(e, 0)).length();
    }
    VertexPair   util_edge_to_vertex_pair    (typename TMeshBase::EHandle e) const { return util_halfedge_to_vertex_pair(get_mesh()->halfedge_handle(e, 0)); }
    PointPair    util_edge_to_point_pair     (typename TMeshBase::EHandle e) const { return util_halfedge_to_point_pair (get_mesh()->halfedge_handle(e, 0)); }
    HalfedgePair util_edge_to_halfedge_pair  (typename TMeshBase::EHandle e) const { auto mesh = get_mesh(); return std::make_pair(mesh->halfedge_handle(e, 0), mesh->halfedge_handle(e, 1)); }
    
    // face-related functions
    std::vector<typename TMeshBase::Point> util_face_to_points(typename TMeshBase::FHandle f) const {
        auto mesh = get_mesh();
        std::vector<typename TMeshBase::Point> result;
        result.reserve(mesh->valence(f));
        for (auto v : mesh->fv_range(f))
            result.push_back(mesh->point(v));
        return result;
    }
    typename TMeshBase::Point util_face_center(typename TMeshBase::FHandle f) const {
        auto mesh = get_mesh();
        typedef typename TMeshBase::Point Point;
        Point result =  Point::vectorized(0);
        for (auto v : mesh->fv_range(f))
            result += mesh->point(v);
        result /= mesh->valence(f);
        return result;
    }
    
    // containment query
    bool util_contains(typename TMeshBase::EHandle e, typename TMeshBase::VHandle v) const {
        auto mesh = get_mesh();
        return
            mesh->to_vertex_handle(mesh->halfedge_handle(e, 0)) == v ||
            mesh->to_vertex_handle(mesh->halfedge_handle(e, 1)) == v;
    }
    bool util_contains(typename TMeshBase::EHandle e, typename TMeshBase::HHandle h) const {
        auto mesh = get_mesh();
        return
            mesh->halfedge_handle(e, 0) == h ||
            mesh->halfedge_handle(e, 1) == h;
    }
    bool util_contains(typename TMeshBase::FHandle f, typename TMeshBase::VHandle v) const { for (auto fv : get_mesh()->fv_range(f)) if (fv == v) return true; return false; }
    bool util_contains(typename TMeshBase::FHandle f, typename TMeshBase::HHandle h) const { for (auto fh : get_mesh()->fh_range(f)) if (fh == h) return true; return false; }
    bool util_contains(typename TMeshBase::FHandle f, typename TMeshBase::EHandle e) const { for (auto fe : get_mesh()->fe_range(f)) if (fe == e) return true; return false; }

    bool util_is_triangle_mesh() const {
        auto mesh = get_mesh();
        for (auto f : mesh->faces())
            if (mesh->valence(f) != 3) return false;
        return true;
    }
    
    typename TMeshBase::VHandle util_add_vertex(const typename TMeshBase::Point& point, const typename TMeshBase::Normal& normal) {
        auto mesh = get_mesh();
        auto v = mesh->add_vertex(point);
        mesh->set_normal(v, normal);
        return v;
    }
private:
    const TMesh* get_mesh() const { return DerivedPtrHolder<TMesh, Utility<TMeshBase, TMesh>>::derived_ptr; }
          TMesh* get_mesh()       { return DerivedPtrHolder<TMesh, Utility<TMeshBase, TMesh>>::derived_ptr; }
};

}
