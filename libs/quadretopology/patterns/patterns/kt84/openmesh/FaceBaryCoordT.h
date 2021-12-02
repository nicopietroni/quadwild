#pragma once
#include <functional>
#include <OpenMesh/Core/Mesh/PolyConnectivity.hh>
#include "../BaryCoordT.h"

namespace kt84 {

template <bool SumToOne = true>
struct FaceBaryCoordT {
    OpenMesh::PolyConnectivity::FHandle f;
    BaryCoordT<SumToOne> bc;                // barycentric coordinates of this point.
    
    FaceBaryCoordT() {}
    FaceBaryCoordT(const OpenMesh::PolyConnectivity::FHandle& f_, const BaryCoordT<SumToOne>& bc_)
        : f(f_)
        , bc(bc_)
    {}
    
    template <class TMesh>
    typename TMesh::Point blend_point(const TMesh& mesh) const {
        auto v = mesh.cfv_iter(f);
        auto p0 = mesh.point(*v);
        auto p1 = mesh.point(*++v);
        auto p2 = mesh.point(*++v);
        return bc.blend(p0, p1, p2);
    }
    
    template <class TMesh>
    typename TMesh::Normal blend_normal(const TMesh& mesh) const {
        static_assert(SumToOne, "blending normals with BaryCoordZero doesn't make sense!");
        
        auto v = mesh.cfv_iter(f);
        auto n0 = mesh.normal(v);
        auto n1 = mesh.normal(++v);
        auto n2 = mesh.normal(++v);
        auto n = bc.blend(n0, n1, n2);
        n.normalize();
        return n;
    }
    
    template <class TMesh, typename TValue>
    TValue blend_value(const TMesh& mesh, const std::function<TValue(const TMesh& mesh, typename TMesh::VHandle v)>& value_getter) const {
        auto v = mesh.cfv_iter(f);
        auto x0 = value_getter(mesh, *v);
        auto x1 = value_getter(mesh, *++v);
        auto x2 = value_getter(mesh, *++v);
        return bc.blend(x0, x1, x2);
    }
};

typedef FaceBaryCoordT<true > FaceBaryCoord;
typedef FaceBaryCoordT<false> FaceBaryCoordZero;

}
