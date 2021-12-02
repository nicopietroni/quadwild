#pragma once
#include <iostream>
#include <vector>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "../DerivedPtrHolder.h"

namespace kt84 {

template <class TMeshBase> struct DebugInfo_Vertex  ;
template <class TMeshBase> struct DebugInfo_Face    ;
template <class TMeshBase> struct DebugInfo_Edge    ;
template <class TMeshBase> struct DebugInfo_Halfedge;

template <class TMeshBase>
struct DebugInfo_Vertex {
    int idx;
    bool deleted;
    typename TMeshBase::VertexData* data;
    typename TMeshBase::Point  point ;
    typename TMeshBase::Normal normal;
    bool is_boundary;
    DebugInfo_Halfedge<TMeshBase>* halfedge;
    std::vector<DebugInfo_Vertex  <TMeshBase>*> vv;
    std::vector<DebugInfo_Face    <TMeshBase>*> vf;
    std::vector<DebugInfo_Edge    <TMeshBase>*> ve;
    std::vector<DebugInfo_Halfedge<TMeshBase>*> voh;
    std::vector<DebugInfo_Halfedge<TMeshBase>*> vih;
    
    DebugInfo_Vertex()
        : idx(-1)
        , deleted()
        , data(nullptr)
        , halfedge(nullptr)
        , is_boundary(false)
    {}
};

template <class TMeshBase>
struct DebugInfo_Face {
    int idx;
    bool deleted;
    typename TMeshBase::FaceData* data;
    typename TMeshBase::Normal normal;
    bool is_boundary;
    DebugInfo_Halfedge<TMeshBase>* halfedge;
    std::vector<DebugInfo_Vertex  <TMeshBase>*> fv;
    std::vector<DebugInfo_Edge    <TMeshBase>*> fe;
    std::vector<DebugInfo_Halfedge<TMeshBase>*> fh;
    std::vector<DebugInfo_Face    <TMeshBase>*> ff;
    
    DebugInfo_Face()
        : idx(-1)
        , deleted()
        , data(nullptr)
        , is_boundary(false)
    {}
};

template <class TTrait>
struct DebugInfo_Face<OpenMesh::TriMesh_ArrayKernelT<TTrait>> {
    typedef OpenMesh::TriMesh_ArrayKernelT<TTrait> MeshBase;
    
    int idx;
    bool deleted;
    typename MeshBase::FaceData* data;
    typename MeshBase::Normal normal;
    bool is_boundary;
    DebugInfo_Halfedge<MeshBase>* halfedge;
    DebugInfo_Vertex  <MeshBase>* fv[3];
    DebugInfo_Edge    <MeshBase>* fe[3];
    DebugInfo_Halfedge<MeshBase>* fh[3];
    DebugInfo_Face    <MeshBase>* ff[3];
    
    DebugInfo_Face()
        : idx(-1)
        , deleted()
        , data(nullptr)
        , is_boundary(false)
    {
        fv[0] = fv[1] = fv[2] = nullptr;
        fe[0] = fe[1] = fe[2] = nullptr;
        fh[0] = fh[1] = fh[2] = nullptr;
        ff[0] = ff[1] = ff[2] = nullptr;
    }
};

template <class TMeshBase>
struct DebugInfo_Edge {
    int idx;
    bool deleted;
    typename TMeshBase::EdgeData* data;
    bool is_boundary;
    DebugInfo_Vertex  <TMeshBase>* vertex  [2];
    DebugInfo_Face    <TMeshBase>* face    [2];
    DebugInfo_Halfedge<TMeshBase>* halfedge[2];
    
    DebugInfo_Edge()
        : idx(-1)
        , deleted()
        , data(nullptr)
        , is_boundary(false)
    {
        vertex  [0] = vertex  [1] = nullptr;
        face    [0] = face    [1] = nullptr;
        halfedge[0] = halfedge[1] = nullptr;
    }
};

template <class TMeshBase>
struct DebugInfo_Halfedge {
    int idx;
    bool deleted;
    typename TMeshBase::HalfedgeData* data;
    bool is_boundary;
    DebugInfo_Vertex  <TMeshBase>* to_vertex  ;
    DebugInfo_Vertex  <TMeshBase>* from_vertex;
    DebugInfo_Face    <TMeshBase>* face       ;
    DebugInfo_Halfedge<TMeshBase>* next       ;
    DebugInfo_Halfedge<TMeshBase>* prev       ;
    DebugInfo_Halfedge<TMeshBase>* opposite   ;
    DebugInfo_Edge    <TMeshBase>* edge       ;
    
    DebugInfo_Halfedge()
        : idx(-1)
        , deleted()
        , data       (nullptr)
        , is_boundary(false)
        , to_vertex  (nullptr)
        , from_vertex(nullptr)
        , face       (nullptr)
        , next       (nullptr)
        , prev       (nullptr)
        , opposite   (nullptr)
        , edge       (nullptr)
    {}
};

namespace internal {
    template <class T>
    void debugInfo_get_common(T& object, typename T::MeshBase& mesh, bool with_normals_v, bool with_normals_f, bool is_verbose) {
        // clear
        object.debugInfo = typename T::Data();
        
        // number of elements
        int nv = mesh.n_vertices ();
        int nf = mesh.n_faces    ();
        int ne = mesh.n_edges    ();
        int nh = mesh.n_halfedges();
        
        if (is_verbose) {
            std::cout << "\n[kt84/openmesh/DebugInfo] ---------------------------\n";
            std::cout << nv << " vertices\n";
            std::cout << nf << " faces\n";
            std::cout << ne << " edge\n";
            std::cout << nh << " halfedges\n";
            std::cout << "--------------------------- [kt84/openmesh/DebugInfo] \n";
        }
        
        // resize arrays
        object.debugInfo.vertices .resize(nv);
        object.debugInfo.faces    .resize(nf);
        object.debugInfo.edges    .resize(ne);
        object.debugInfo.halfedges.resize(nh);
        
        // vertex
        for (int i = 0; i < nv; ++i) {
            auto v = mesh.vertex_handle(i);
            auto& v_dbg = object.debugInfo.vertices[i];
            v_dbg.idx = v.idx();
            if (mesh.has_vertex_status()) v_dbg.deleted = mesh.status(v).deleted();
            v_dbg.data = &mesh.data(v);
            v_dbg.point = mesh.point(v);
            v_dbg.is_boundary = mesh.is_boundary(v);
            if (with_normals_v) v_dbg.normal = mesh.normal(v);
            if (mesh.halfedge_handle(v).is_valid())
                v_dbg.halfedge = &object.debugInfo.halfedges[mesh.halfedge_handle(v).idx()];
            for (auto vv  = mesh.vv_iter (v); vv .is_valid(); ++vv ) v_dbg.vv .push_back(&object.debugInfo.vertices [vv ->idx()]);
            for (auto vf  = mesh.vf_iter (v); vf .is_valid(); ++vf ) v_dbg.vf .push_back(&object.debugInfo.faces    [vf ->idx()]);
            for (auto ve  = mesh.ve_iter (v); ve .is_valid(); ++ve ) v_dbg.ve .push_back(&object.debugInfo.edges    [ve ->idx()]);
            for (auto voh = mesh.voh_iter(v); voh.is_valid(); ++voh) v_dbg.voh.push_back(&object.debugInfo.halfedges[voh->idx()]);
            for (auto vih = mesh.vih_iter(v); vih.is_valid(); ++vih) v_dbg.vih.push_back(&object.debugInfo.halfedges[vih->idx()]);
        }
        
        // edge
        for (int i = 0; i < ne; ++i) {
            auto e = mesh.edge_handle(i);
            auto& e_dbg = object.debugInfo.edges[i];
            e_dbg.idx = e.idx();
            if (mesh.has_edge_status()) e_dbg.deleted = mesh.status(e).deleted();
            e_dbg.data = &mesh.data(e);
            e_dbg.is_boundary = mesh.is_boundary(e);
            for (int j = 0; j < 2; ++j) {
                auto h = mesh.halfedge_handle(e, j);
                e_dbg.vertex  [j] = &object.debugInfo.vertices [mesh.to_vertex_handle(h).idx()];
                e_dbg.halfedge[j] = &object.debugInfo.halfedges[h.idx()];
                
                auto f = mesh.face_handle(h);
                e_dbg.face[j] = f.is_valid() ? &object.debugInfo.faces[f.idx()] : 0;
            }
        }
        
        // halfedge
        for (int i = 0; i < nh; ++i) {
            auto h = mesh.halfedge_handle(i);
            auto& h_dbg = object.debugInfo.halfedges[i];
            h_dbg.idx = h.idx();
            if (mesh.has_halfedge_status()) h_dbg.deleted = mesh.status(h).deleted();
            h_dbg.data = &mesh.data(h);
            h_dbg.is_boundary = mesh.is_boundary(h);
            h_dbg.to_vertex   = &object.debugInfo.vertices [mesh.to_vertex_handle        (h).idx()];
            h_dbg.from_vertex = &object.debugInfo.vertices [mesh.from_vertex_handle      (h).idx()];
            h_dbg.next        = &object.debugInfo.halfedges[mesh.next_halfedge_handle    (h).idx()];
            h_dbg.prev        = &object.debugInfo.halfedges[mesh.prev_halfedge_handle    (h).idx()];
            h_dbg.opposite    = &object.debugInfo.halfedges[mesh.opposite_halfedge_handle(h).idx()];
            h_dbg.edge        = &object.debugInfo.edges    [mesh.edge_handle             (h).idx()];
            
            auto f = mesh.face_handle(h);
            h_dbg.face = f.is_valid() ? &object.debugInfo.faces[f.idx()] : 0;
        }
    }
}

template <class TMeshBase, class TMesh>
struct DebugInfo : public DerivedPtrHolder<TMesh, DebugInfo<TMeshBase, TMesh>> {
    typedef TMeshBase MeshBase;
    
    static_assert(MeshBase::Point::size_ == 3, "TMeshBase::Point should be 3D vector!");
    
    struct Data {
        std::vector<DebugInfo_Vertex  <MeshBase>> vertices ;
        std::vector<DebugInfo_Face    <MeshBase>> faces    ;
        std::vector<DebugInfo_Edge    <MeshBase>> edges    ;
        std::vector<DebugInfo_Halfedge<MeshBase>> halfedges;
    } debugInfo;
    
    void debugInfo_get(bool with_normals_v = false, bool with_normals_f = false, bool is_verbose = false) {
        TMesh* mesh = DerivedPtrHolder<TMesh, DebugInfo<TMeshBase, TMesh>>::derived_ptr;
        
        internal::debugInfo_get_common(*this, *mesh, with_normals_v, with_normals_f, is_verbose);
        
        int nf = mesh->n_faces();
        
        // face
        for (int i = 0; i < nf; ++i) {
            auto f = mesh->face_handle(i);
            auto& f_dbg = debugInfo.faces[i];
            f_dbg.idx = f.idx();
            if (mesh->has_face_status()) f_dbg.deleted = mesh->status(f).deleted();
            f_dbg.data = &mesh->data(f);
            f_dbg.is_boundary = mesh->is_boundary(f);
            if (with_normals_f) f_dbg.normal = mesh->normal(f);
            f_dbg.halfedge = &debugInfo.halfedges[mesh->halfedge_handle(f).idx()];
            int n_sides = 0;
            for (auto fh = mesh->fh_iter(f); fh.is_valid(); ++fh) ++n_sides;
            f_dbg.fv.resize(n_sides, 0);
            f_dbg.ff.resize(n_sides, 0);
            f_dbg.fe.resize(n_sides, 0);
            f_dbg.fh.resize(n_sides, 0);
            auto fh = mesh->fh_iter(f);
            for (int j = 0; j < n_sides; ++j, ++fh) {
                auto v = mesh->to_vertex_handle(*fh);
                auto e = mesh->edge_handle     (*fh);
                auto f = mesh->face_handle     (mesh->opposite_halfedge_handle(*fh));
                f_dbg.fv[j] = &debugInfo.vertices [v .idx()];
                f_dbg.ff[j] = f.idx() == -1 ? nullptr : &debugInfo.faces[f.idx()];
                f_dbg.fe[j] = &debugInfo.edges    [e .idx()];
                f_dbg.fh[j] = &debugInfo.halfedges[fh->idx()];
            }
        }
    }
};

template <class TTrait, class TMesh>
struct DebugInfo<OpenMesh::TriMesh_ArrayKernelT<TTrait>, TMesh> : public DerivedPtrHolder<TMesh, DebugInfo<OpenMesh::TriMesh_ArrayKernelT<TTrait>, TMesh>> {
    typedef OpenMesh::TriMesh_ArrayKernelT<TTrait> MeshBase;
    
    struct Data {
        std::vector<DebugInfo_Vertex  <MeshBase>> vertices ;
        std::vector<DebugInfo_Face    <MeshBase>> faces    ;
        std::vector<DebugInfo_Edge    <MeshBase>> edges    ;
        std::vector<DebugInfo_Halfedge<MeshBase>> halfedges;
    } debugInfo;
    
    void debugInfo_get(bool with_normals_v = false, bool with_normals_f = false, bool is_verbose = false) {
        TMesh* mesh = DerivedPtrHolder<TMesh, DebugInfo<OpenMesh::TriMesh_ArrayKernelT<TTrait>, TMesh>>::derived_ptr;
        
        internal::debugInfo_get_common(*this, *mesh, with_normals_v, with_normals_f, is_verbose);
        
        int nf = mesh->n_faces();
        
        // face
        for (int i = 0; i < nf; ++i) {
            auto f = mesh->face_handle(i);
            auto& f_dbg = debugInfo.faces[i];
            f_dbg.idx = f.idx();
            if (mesh->has_face_status()) f_dbg.deleted = mesh->status(f).deleted();
            f_dbg.data = &mesh->data(f);
            f_dbg.is_boundary = mesh->is_boundary(f);
            if (with_normals_f) f_dbg.normal = mesh->normal(f);
            f_dbg.halfedge = &debugInfo.halfedges[mesh->halfedge_handle(f).idx()];
            auto fh = mesh->fh_iter(f);
            for (int j = 0; j < 3; ++j, ++fh) {
                auto v = mesh->to_vertex_handle(*fh);
                auto e = mesh->edge_handle     (*fh);
                auto f = mesh->face_handle     (mesh->opposite_halfedge_handle(*fh));
                f_dbg.fv[j] = &debugInfo.vertices [v .idx()];
                f_dbg.ff[j] = f.idx() == -1 ? nullptr : &debugInfo.faces[f.idx()];
                f_dbg.fe[j] = &debugInfo.edges    [e .idx()];
                f_dbg.fh[j] = &debugInfo.halfedges[fh->idx()];
            }
        }
    }
};

}
