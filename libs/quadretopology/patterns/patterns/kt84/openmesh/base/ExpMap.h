#pragma once
#include <vector>
#include <queue>
#include <Eigen/Geometry>
#include <boost/optional.hpp>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyConnectivity.hh>
#include "../vector_convert.h"
#include "../FaceBaryCoordT.h"
#include "../../eigen_util.h"
#include "../../util.h"
#include "../DerivedPtrHolder.h"

namespace kt84 {

/*
    Implementation of a surface mesh parameterization algorithm descibed in the following articles:
        Ryan Schmidt, Cindy Grimm, and Brian Wyvill.
        Interactive decal compositing with discrete exponential maps.
        ACM SIGGRAPH 2006.
        
        Ryan Schmidt and Karan Singh.
        Drag, Drop, and Clone: An Interactive Interface for Surface Composition.
        Technical Report CSRG-611, Department of Computer Science, University of Toronto, 2010.
    
    NOTE: The algorithm is implemented in a slightly different way; the uv coordiante of currently processed vertex q is estimated by rotating 
    its local frame to that of its adjacent vertex p that has been processed already, which is in contrast to what is claimed in the original 
    SIGGRAPH paper (see sentences below equation (3)). This modification seems to provide better robustness especially when parameterizing highly 
    curved surface regions (probably because normals at distant locations are very different from the seed's normal).

    parameters:
        seed(s)                 seed vertex id
        seed(s)_basis_u         3D vector in object space corresponding to unit basis vector of u-coordinate.
                                (its length determines the overall scaling of the parameterization)
        seed(s)_uv              seed uv coordinate
        dist_max                termination criteria based on object-space distance
        uv_min, uv_max          termination criteria based on parameter-space bounding box
*/

struct ExpMap_VertexTraits {
    struct Data {
        bool            paramed;
        double          dist;
        double          weight;
        Eigen::Vector2d uv;
        Eigen::Vector2d basis_u;              // represents u-coordiante basis vector in object space.
                                            // (i.e., its length determines scaling of parameterization.)
                                            // living on the tangent plane defined by intrinsic bases
        
        Data()
            : paramed()
            , dist   (util::dbl_max())
            , weight ()
            , uv     (Eigen::Vector2d::Zero())
            , basis_u(Eigen::Vector2d::Zero())
        {}
    } expmap;
};

struct ExpMap_FaceTraits {
    bool expmap_paramed;
    
    ExpMap_FaceTraits()
        : expmap_paramed()
    {}
};

namespace internal {
    template <class TMeshBase>
    void expmap_compute_common(TMeshBase* mesh,
                               const std::vector<OpenMesh::PolyConnectivity::VHandle>& seeds,
                               const std::vector<Eigen::Vector3d                    >& seeds_basis_u,
                               const std::vector<Eigen::Vector2d                    >& seeds_uv,
                               double dist_max,
                               const Eigen::Vector2d& uv_min,
                               const Eigen::Vector2d& uv_max,
                               std::vector<OpenMesh::PolyConnectivity::FHandle>& expmap_paramed_faces)
    {
        // helper functions ==============================================================================================
        auto get_intr_basis = [] (const Eigen::Vector3d& normal) -> Eigen::Vector3d {
            // get unit basis vector orthogonal to input vector in intrinsical manner
            double nx = std::abs(normal.x());
            double ny = std::abs(normal.y());
            double nz = std::abs(normal.z());
            Eigen::Vector3d basis =
                nx < ny && ny < nz ? Eigen::Vector3d::UnitX() :
                ny < nz            ? Eigen::Vector3d::UnitY() :
                                     Eigen::Vector3d::UnitZ();
            eigen_util::orthonormalize(normal, basis);
            return basis;
        };
        
        auto compute_rotation = [get_intr_basis] (const Eigen::Vector3d& p_normal,
                                                  const Eigen::Vector3d& q_normal) -> Eigen::Rotation2Dd
            // compute 2D rotation which will be used to rotate basis_u of p to estimate that of q
        {
            // get intrinsic basis vectors of p, q
            auto p_intr_basis = get_intr_basis(p_normal);
            auto q_intr_basis = get_intr_basis(q_normal);
            
            // step 1: compute 3D rotation that aligns q_normal with p_normal. apply it to q_basis.
            auto rot1_axis = q_normal.cross(p_normal);
            if (!rot1_axis.isZero()) {
                rot1_axis.normalize();
                double rot1_angle = util::acos_clamped(q_normal.dot(p_normal));
                q_intr_basis = Eigen::AngleAxisd(rot1_angle, rot1_axis) * q_intr_basis;
            }
            
            // step 2: compute 2D rotation that makes (rotated) q_basis aligned with p_basis.
            double rot2_angle = util::acos_clamped(q_intr_basis.dot(p_intr_basis));
            if (q_intr_basis.cross(p_intr_basis).dot(p_normal) < 0)     // negative rotation
                rot2_angle *= -1.0f;
            
            return Eigen::Rotation2Dd(rot2_angle);
        };
        
        auto compute_uv = [get_intr_basis] (const Eigen::Vector3d& p_point ,
                                            const Eigen::Vector3d& p_normal,
                                            const Eigen::Vector2d& p_uv,
                                            const Eigen::Vector2d& p_basis_u,
                                            const Eigen::Vector3d& q_point ) -> Eigen::Vector2d
            // computes uv for point q using p's point, normal, uv, basis_u
        {
            auto p_intr_basis_x = get_intr_basis(p_normal);
            auto p_intr_basis_y = p_normal.cross(p_intr_basis_x);
            
            // project (rotate) pq into tangent plane
            Eigen::Vector3d pq = q_point - p_point;
            Eigen::Vector2d local_uv(pq.dot(p_intr_basis_x), pq.dot(p_intr_basis_y));
            local_uv *= pq.norm() / local_uv.norm();
            
            // encode 2D local_uv using propagated (non-unit) basis vectors
            auto p_basis_v = eigen_util::rotate90(p_basis_u);
            local_uv = Eigen::Vector2d(p_basis_u.dot(local_uv), p_basis_v.dot(local_uv));
            
            // account for scale difference
            local_uv /= p_basis_u.squaredNorm();
            
            Eigen::Vector2d q_uv = p_uv + local_uv;
            
            return q_uv;
        };
        // ============================================================================================== helper functions
        
        // clear all vertex data
        for (auto v : mesh->vertices())
            mesh->data(v).expmap = ExpMap_VertexTraits::Data();
        
        // priority queue for propagation front
        struct QueueElement {
            OpenMesh::PolyConnectivity::VHandle vhandle;
            double dist;
            
            QueueElement() : dist(){}
            QueueElement(OpenMesh::PolyConnectivity::VHandle vhandle_, double dist_) : vhandle(vhandle_) , dist(dist_) {}
            bool operator<(const QueueElement& rhs) const { return dist > rhs.dist; }       // note: smaller distance gets higher priority
        };
        std::priority_queue<QueueElement> candidates;
        
        // init
        for (size_t i = 0; i < seeds.size(); ++i) {
            auto seed = seeds[i];
            auto& seed_data = mesh->data(seed).expmap;
            
            Eigen::Vector3d seed_normal       = o2e(mesh->normal(seed));
            Eigen::Vector3d seed_intr_basis_x = get_intr_basis(seed_normal);
            Eigen::Vector3d seed_intr_basis_y = seed_normal.cross(seed_intr_basis_x);
            Eigen::Vector3d seed_basis_u      = seeds_basis_u[i];
            
            seed_data.weight  += 1;
            seed_data.dist     = 0;
            seed_data.uv      += seeds_uv[i];
            seed_data.basis_u += Eigen::Vector2d(seed_intr_basis_x.dot(seed_basis_u), seed_intr_basis_y.dot(seed_basis_u));
            
            candidates.push(QueueElement(seed, 0));
        }
        
        while (!candidates.empty()) {
            auto p = candidates.top().vhandle;
            candidates.pop();
            
            auto& p_data   = mesh->data(p).expmap;
            auto  p_point  = o2e(mesh->point (p));
            auto  p_normal = o2e(mesh->normal(p));
            
            if (p_data.paramed) continue;
            p_data.paramed = true;
            
            p_data.uv      /= p_data.weight;
            p_data.basis_u /= p_data.weight;
            
            if (dist_max < p_data.dist) continue;
            
            if (!Eigen::AlignedBox2d(uv_min, uv_max).contains(p_data.uv)) continue;
            
            for (auto q = mesh->vv_iter(p); q.is_valid(); ++q) {
                auto& q_data   = mesh->data(*q).expmap;
                auto  q_point  = o2e(mesh->point (*q));
                auto  q_normal = o2e(mesh->normal(*q));
                
                if (q_data.paramed) continue;
                
                // distance update
                double dist_pq = (p_point - q_point).norm();
                if (p_data.dist + dist_pq < q_data.dist)
                    q_data.dist = p_data.dist + dist_pq;
                
                // weight based on distance between p and q
                double weight = 1.0f / std::pow(dist_pq, 0.25);
                q_data.weight += weight;
                
                // estimate of data by weighted averaging
                q_data.uv      += weight * compute_uv(p_point, p_normal, p_data.uv, p_data.basis_u, q_point);
                q_data.basis_u += weight * (compute_rotation(p_normal, q_normal) * p_data.basis_u);            // be sure to use parentheses! otherwise Eigen will happily go funny:)
                
                candidates.push(QueueElement(*q, q_data.dist));
            }
        }
        
        // collect faces with its vertices all paramed
        expmap_paramed_faces.clear();
        expmap_paramed_faces.reserve(mesh->n_faces());
        for (auto f : mesh->faces()) {
            mesh->data(f).expmap_paramed = true;
            for (auto v = mesh->fv_iter(f); v.is_valid(); ++v) {
                if (!mesh->data(*v).expmap.paramed) {
                    mesh->data(f).expmap_paramed = false;
                    break;
                }
            }
            if (mesh->data(f).expmap_paramed)
                expmap_paramed_faces.push_back(f);
        }
    }
}

template <class TMeshBase, class TMesh>
struct ExpMap : public DerivedPtrHolder<TMesh, ExpMap<TMeshBase, TMesh>> {
    std::vector<typename TMeshBase::FHandle> expmap_paramed_faces;
    
    void expmap_compute(
        const std::vector<typename TMeshBase::VHandle>& seeds,
        const std::vector<Eigen::Vector3d            >& seeds_basis_u,
        const std::vector<Eigen::Vector2d            >& seeds_uv,
        double dist_max,
        const Eigen::Vector2d& uv_min = Eigen::Vector2d::Constant(-std::numeric_limits<double>::max()),
        const Eigen::Vector2d& uv_max = Eigen::Vector2d::Constant( std::numeric_limits<double>::max()))
    {
        TMesh* mesh = DerivedPtrHolder<TMesh, ExpMap<TMeshBase, TMesh>>::derived_ptr;
        internal::expmap_compute_common(mesh, seeds, seeds_basis_u, seeds_uv, dist_max, uv_min, uv_max, expmap_paramed_faces);
    }
    void expmap_compute(
        typename TMeshBase::VHandle seed,
        Eigen::Vector3d             seed_basis_u,
        Eigen::Vector2d             seed_uv,
        double dist_max,
        const Eigen::Vector2d& uv_min = Eigen::Vector2d::Constant(-std::numeric_limits<double>::max()),
        const Eigen::Vector2d& uv_max = Eigen::Vector2d::Constant( std::numeric_limits<double>::max()))
    {
        expmap_compute(
            std::vector<typename TMeshBase::VHandle>(1, seed        ),
            std::vector<Eigen::Vector3d            >(1, seed_basis_u),
            std::vector<Eigen::Vector2d            >(1, seed_uv     ),
            dist_max,
            uv_min,
            uv_max);
    }
};

// specialization for triangle meshes, with functionality to convert on-surface position (face id + barycentric coord) into uv coordinates.
template <class TTrait, class TMesh>
struct ExpMap<OpenMesh::TriMesh_ArrayKernelT<TTrait>, TMesh> : public DerivedPtrHolder<TMesh, ExpMap<OpenMesh::TriMesh_ArrayKernelT<TTrait>, TMesh>> {
    typedef OpenMesh::TriMesh_ArrayKernelT<TTrait> MeshBase;
    
    std::vector<typename MeshBase::FHandle> expmap_paramed_faces;
    
    void expmap_compute(
        const std::vector<typename MeshBase::VHandle>& seeds,
        const std::vector<Eigen::Vector3d           >& seeds_basis_u,
        const std::vector<Eigen::Vector2d           >& seeds_uv,
        double dist_max,
        const Eigen::Vector2d& uv_min = Eigen::Vector2d::Constant(-std::numeric_limits<double>::max()),
        const Eigen::Vector2d& uv_max = Eigen::Vector2d::Constant( std::numeric_limits<double>::max()))
    {
        TMesh* mesh = get_mesh();
        internal::expmap_compute_common(mesh, seeds, seeds_basis_u, seeds_uv, dist_max, uv_min, uv_max, expmap_paramed_faces);
    }
    void expmap_compute(
        typename MeshBase::VHandle seed,
        Eigen::Vector3d            seed_basis_u,
        Eigen::Vector2d            seed_uv,
        double dist_max,
        const Eigen::Vector2d& uv_min = Eigen::Vector2d::Constant(-std::numeric_limits<double>::max()),
        const Eigen::Vector2d& uv_max = Eigen::Vector2d::Constant( std::numeric_limits<double>::max()))
    {
        expmap_compute(
            std::vector<typename MeshBase::VHandle>(1, seed        ),
            std::vector<Eigen::Vector3d           >(1, seed_basis_u),
            std::vector<Eigen::Vector2d           >(1, seed_uv     ),
            dist_max,
            uv_min,
            uv_max);
    }
    void expmap_compute(
        const std::vector<FaceBaryCoord  >& seeds,
        const std::vector<Eigen::Vector3d>& seeds_basis_u,
        const std::vector<Eigen::Vector2d>& seeds_uv,
        double dist_max,
        const Eigen::Vector2d& uv_min = Eigen::Vector2d::Constant(-std::numeric_limits<double>::max()),
        const Eigen::Vector2d& uv_max = Eigen::Vector2d::Constant( std::numeric_limits<double>::max()))
    {
        TMesh* mesh = get_mesh();
        
        std::vector<typename MeshBase::VHandle> seeds2;
        std::vector<Eigen::Vector3d           > seeds2_basis_u;
        std::vector<Eigen::Vector2d           > seeds2_uv;
        
        for (size_t i = 0; i < seeds.size(); ++i) {
            auto& seed         = seeds        [i];
            auto& seed_basis_u = seeds_basis_u[i];
            auto& seed_uv      = seeds_uv     [i];
            
            Eigen::Vector3d p = o2e(seed.blend_point (*mesh));
            Eigen::Vector3d n = o2e(seed.blend_normal(*mesh));
            
            Eigen::Vector3d seed_basis_v = n.cross(seed_basis_u);
            
            double offset_scaling = 1 / seed_basis_u.squaredNorm();
            
            for (auto v = mesh->fv_iter(seed.f); v; ++v) {
                Eigen::Vector3d pq = o2e(mesh->point(v)) - p;
                Eigen::Vector2d uv_offset(pq.dot(seed_basis_u), pq.dot(seed_basis_v));
                uv_offset *= offset_scaling;
                
                seeds2        .push_back(v);
                seeds2_basis_u.push_back(seed_basis_u);
                seeds2_uv     .push_back(seed_uv + uv_offset);
            }
        }
        
        expmap_compute(seeds2, seeds2_basis_u, seeds2_uv, dist_max, uv_min, uv_max);
    }
    void expmap_compute(
        FaceBaryCoord   seed,
        Eigen::Vector3d seed_basis_u,
        Eigen::Vector2d seed_uv,
        double dist_max,
        const Eigen::Vector2d& uv_min = Eigen::Vector2d::Constant(-std::numeric_limits<double>::max()),
        const Eigen::Vector2d& uv_max = Eigen::Vector2d::Constant( std::numeric_limits<double>::max()))
    {
        expmap_compute(
            std::vector<FaceBaryCoord  >(1, seed        ),
            std::vector<Eigen::Vector3d>(1, seed_basis_u),
            std::vector<Eigen::Vector2d>(1, seed_uv     ),
            dist_max,
            uv_min,
            uv_max);
    }
    // conversion from on-surface position to uv coordinates
    boost::optional<Eigen::Vector2d> fbc_to_uv(const FaceBaryCoord& fbc) const {
        TMesh* mesh = get_mesh();
        
        if (!mesh->data(fbc.f).expmap_paramed)
            return boost::none;
        
        return fbc.blend_value<MeshBase, Eigen::Vector2d>(
            *mesh,
            [](const MeshBase& mesh, typename MeshBase::VHandle v) { return mesh.data(v).expmap.uv; }
        );
    }
    boost::optional<Eigen::Vector2d> fbc_to_uv(typename MeshBase::FHandle f, const BaryCoord& bc) const {
        return fbc_to_uv(FaceBaryCoord(f, bc));
    }
    FaceBaryCoord uv_to_fbc(Eigen::Vector2d uv) const {
        TMesh* mesh = get_mesh();
        
        for (auto f : expmap_paramed_faces) {
            // compute bounding box of this face
            Eigen::AlignedBox2d box;
            std::vector<Eigen::Vector2d> face_uv;
            for (auto v = mesh->cfv_iter(f); v.is_valid(); ++v) {
                auto expmap_uv = mesh->data(*v).expmap.uv;
                face_uv.push_back(expmap_uv);
                box.extend(expmap_uv);
            }
            
            // skip if uv is outside the bounding box
            if (!box.contains(uv)) continue;
            
            // compute barycentric coordinate
            // (1 - hit.u - hit.v) * face_uv[0] + hit.u * face_uv[1] + hit.v * face_uv[2] = uv
            // | face_uv[1] - face_uv[0], face_uv[2] - face_uv[0] | * |hit.u| = | uv - face_uv[0] |
            // |                                                  |   |hit.v|   |                 |
            Eigen::Matrix2d A;
            A << face_uv[1] - face_uv[0], face_uv[2] - face_uv[0];
            Eigen::Vector2d hit_uv = A.inverse() * (uv - face_uv[0]);
            
            // skip if barycentric coordinate is negative
            BaryCoord bc(1 - hit_uv[0] - hit_uv[1], hit_uv[0], hit_uv[1]);
            if (!bc.is_all_positive()) continue;
            
            // return result
            return FaceBaryCoord(f, bc);
        }
        
        return FaceBaryCoord();
    }
private:
    TMesh* get_mesh() const { return DerivedPtrHolder<TMesh, ExpMap<OpenMesh::TriMesh_ArrayKernelT<TTrait>, TMesh>>::derived_ptr; }
};

}
