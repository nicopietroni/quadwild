#pragma once
#include <vector>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyConnectivity.hh>
#include "../../StdUtil.h"
#include "../../EigenUtil.h"
#include "../vector_convert.h"

/*
    Implementation of the following article:
        Shi-Qing Xin, Ying He, Chi-Wing Fu.
        Efficiently Computing Exact Geodesic Loops within Finite Steps.
        IEEE Transactions on Visualization and Computer Graphics.
        vol.18, no.6, pp.879-889, June 2012.
        doi: 10.1109/TVCG.2011.119
    
    The paper's webpage: https://sites.google.com/site/xinshiqing/
}*/

namespace kt84 {

struct GeodesicLoop_VertexTraits {
    enum class StripSide {
        UNDEFINED = -1,
        LEFT,
        RIGHT
    };
    
    struct Data {
        StripSide stripSide;
        Eigen::Vector2d uv;
        Eigen::Vector2d uv_duplicate;       // for vertices on the cutting edge. corresponding to the duplicated vertex at the end of unfolded triangle strip.
        
        Data()
            : stripSide(StripSide::UNDEFINED)
            , uv(0, 0)
            , uv_duplicate(0, 0)
        {}
    } geodesicLoop;
};

struct GeodesicLoop_EdgeTraits {
    enum class VertexSide {
        UNDEFINED = -1,
        COINCIDE,
        LEFT,
        RIGHT
    } geodesicLoop_vertexSide;
    
    GeodesicLoop_EdgeTraits()
        : geodesicLoop_vertexSide(VertexSide::UNDEFINED)
    {}
};

template <class TMeshBase>
struct GeodesicLoop;

template <class TTrait>
struct GeodesicLoop<OpenMesh::TriMesh_ArrayKernelT<TTrait>> {
    typedef OpenMesh::TriMesh_ArrayKernelT<TTrait> MeshBase;
    
    struct PathPoint {
        typename MeshBase::VHandle vertex;
        typename MeshBase::EHandle edge;
        double edge_coordinate;         // supposing vi corresponds to to_vertex(halfedge(edge, i)) where i = 0 or 1, 
                                        // this represents on-edge point as: (1 - edge_coordinate) * v0 + edge_coordinate * v1
        enum class Type {
            UNDEFINED,
            VERTEX,
            EDGE
        };
        Type type() const {
            return
                 vertex.is_valid() && !edge.is_valid() ? Type::VERTEX :
                !vertex.is_valid() &&  edge.is_valid() ? Type::EDGE   : Type::UNDEFINED;
        }
        PathPoint() : edge_coordinate(-1) {}
        PathPoint(typename MeshBase::VHandle vertex_) : vertex(vertex_), edge_coordinate(-1) {}
        PathPoint(typename MeshBase::EHandle edge_, double edge_coordinate_) : edge(edge_), edge_coordinate(edge_coordinate_) {}
        bool operator==(const PathPoint& rhs) const {
            return
                (type() == Type::VERTEX && vertex == rhs.vertex) ||
                (type() == Type::EDGE   && edge   == rhs.edge   && edge_coordinate == rhs.edge_coordinate);
        }
    };
    typedef std::vector<typename MeshBase::EHandle> EdgeSequence;
    typedef std::vector<PathPoint> Path;
    
    virtual ~GeodesicLoop() {}
    Path geodesicLoop_compute(const EdgeSequence& closed_edge_sequence) {
        if (closed_edge_sequence.front() != closed_edge_sequence.back())
            // semantics of closed sequence/loop: front and back are the same
            return Path();
        
        auto mesh = dynamic_cast<MeshBase*>(this);
        assert(mesh);
        
        auto E = closed_edge_sequence;
        
        // helper functions=====================================================================================================================================
        auto determine_side_around_vertex = [mesh] (typename MeshBase::VHandle v, const PathPoint& point_prev, const PathPoint& point_next) -> void {
            // init
            for (auto e = mesh->ve_iter(v); e; ++e)
                mesh->data(e).geodesicLoop_vertexSide = GeodesicLoop_EdgeTraits::VertexSide::UNDEFINED;
            
            // set side information immediately available from point_prev and point_next
            for (auto h = mesh->voh_iter(v); h; ++h) {
                auto v_to            = mesh->to_vertex_handle(h);
                auto e               = mesh->edge_handle(h);
                auto e_prev          = mesh->edge_handle(mesh->prev_halfedge_handle(h));
                auto e_next          = mesh->edge_handle(mesh->prev_halfedge_handle(h));
                auto e_opposite_next = mesh->edge_handle(mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(h)));
                
                bool prev_is_left = v_to == point_next.vertex;
                
                if (v_to == point_prev.vertex || v_to == point_next.vertex) {
                    mesh->data(e              ).geodesicLoop_vertexSide = GeodesicLoop_EdgeTraits::VertexSide::COINCIDE;
                    mesh->data(e_prev         ).geodesicLoop_vertexSide = prev_is_left ? GeodesicLoop_EdgeTraits::VertexSide::LEFT  : GeodesicLoop_EdgeTraits::VertexSide::RIGHT;
                    mesh->data(e_opposite_next).geodesicLoop_vertexSide = prev_is_left ? GeodesicLoop_EdgeTraits::VertexSide::RIGHT : GeodesicLoop_EdgeTraits::VertexSide::LEFT ;
                    
                } else if (e_next == point_prev.edge || e_next == point_next.edge) {
                    mesh->data(e_prev).geodesicLoop_vertexSide = prev_is_left ? GeodesicLoop_EdgeTraits::VertexSide::LEFT  : GeodesicLoop_EdgeTraits::VertexSide::RIGHT;
                    mesh->data(e     ).geodesicLoop_vertexSide = prev_is_left ? GeodesicLoop_EdgeTraits::VertexSide::RIGHT : GeodesicLoop_EdgeTraits::VertexSide::LEFT ;
                }
            }
            
            // find halfedge outgoing from v which is not yet assigned side info with its previous one side info assigned already
            auto h_start = mesh->halfedge_handle(v);
            for (; ; ) {
                auto e = mesh->edge_handle(h_start);
                if (mesh->data(e).geodesicLoop_vertexSide != GeodesicLoop_EdgeTraits::VertexSide::UNDEFINED)
                    continue;
                
                auto e_prev = mesh->edge_handle(mesh->prev_halfedge_handle(h_start));
                if (mesh->data(e_prev).geodesicLoop_vertexSide != GeodesicLoop_EdgeTraits::VertexSide::UNDEFINED)
                    break;
                
                h_start = mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(h_start));
            }
            
            // set side flag for other halfedges by propagation
            for (auto h = h_start; ; ) {
                auto e = mesh->edge_handle(h);
                auto& vertexSide = mesh->data(e).geodesicLoop_vertexSide;
                
                if (vertexSide == GeodesicLoop_EdgeTraits::VertexSide::UNDEFINED) {
                    // assign the same side as the previous edge
                    auto e_prev = mesh->edge_handle(mesh->prev_halfedge_handle(h_start));
                    vertexSide = mesh->data(e_prev).geodesicLoop_vertexSide;
                }
                
                h = mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(h));
                if (h == h_start)
                    break;
            }
        };
        auto get_edge_point = [mesh] (const PathPoint& point) -> Eigen::Vector3d {
            assert(point.type() == PathPoint::Type::EDGE);
            auto v0v1 = o2e(OpenMeshUtil::edge_to_point_pair(*mesh, point.edge));
            return (1 - point.edge_coordinate) * v0v1.first + point.edge_coordinate * v0v1.second;
        };
        auto compute_incident_angles = [mesh, get_edge_point] (typename MeshBase::VHandle v, const PathPoint& point_prev, const PathPoint& point_next)-> std::pair<double, double> {
            auto result = std::make_pair<double, double>(0, 0);
            
            Eigen::Vector3d edge_point_prev = get_edge_point(point_prev);
            Eigen::Vector3d edge_point_next = get_edge_point(point_next);
            
            for (auto h = mesh->voh_iter(v); h; ++h) {
                auto h_next = mesh->next_halfedge_handle(h);
                auto h_prev = mesh->prev_halfedge_handle(h);
            
                Eigen::Vector3d d0 = -o2e(OpenMeshUtil::halfedge_to_vector(*mesh, h_prev)).normalized();
                Eigen::Vector3d d1 =  o2e(OpenMeshUtil::halfedge_to_vector(*mesh, h     )).normalized();
                
                auto side0 = mesh->data(mesh->edge_handle(h_prev)).geodesicLoop_vertexSide;
                auto side1 = mesh->data(mesh->edge_handle(h     )).geodesicLoop_vertexSide;
                
                if (side0 == GeodesicLoop_EdgeTraits::VertexSide::LEFT && side1 == GeodesicLoop_EdgeTraits::VertexSide::RIGHT) {
                    // this face is at transition from left side to right side
                    assert(point_next.edge == mesh->edge_handle(h_next));
                    
                    Eigen::Vector3d dm = (edge_point_next - o2e(mesh->point(v))).normalized();
                    
                    double angle_left  = StdUtil::acos_clamped(d0.dot(dm));
                    double angle_right = StdUtil::acos_clamped(d1.dot(dm));
                    
                    result.first  += angle_left ;
                    result.second += angle_right;
                    
                } else if (side0 == GeodesicLoop_EdgeTraits::VertexSide::RIGHT && side1 == GeodesicLoop_EdgeTraits::VertexSide::LEFT) {
                    // this face is at transition from right side to left side
                    assert(point_prev.edge == mesh->edge_handle(h_next));
                    
                    Eigen::Vector3d dm = (edge_point_prev - o2e(mesh->point(v))).normalized();
                    
                    double angle_left  = StdUtil::acos_clamped(d1.dot(dm));
                    double angle_right = StdUtil::acos_clamped(d0.dot(dm));
                    
                    result.first  += angle_left ;
                    result.second += angle_right;
                    
                
                } else {
                    // simply add angle subtending the face to appropriate side
                    double angle = StdUtil::acos_clamped(d0.dot(d1));
                    
                    auto side = side0 == GeodesicLoop_EdgeTraits::VertexSide::COINCIDE ? side1 : side0;
                    
                    if (side == GeodesicLoop_EdgeTraits::VertexSide::LEFT)
                        result.first  += angle;
                    else
                        result.second += angle;
                }
            }
            
            return result;
        };
        auto is_inward_convex         = [mesh, determine_side_around_vertex, compute_incident_angles] (const EdgeSequence& E, typename MeshBase::VHandle v, const PathPoint& point_prev, const PathPoint& point_next) -> bool {
            determine_side_around_vertex(v, point_prev, point_next);
            auto angles = compute_incident_angles(v, point_prev, point_next);
            
            auto vertexSide = GeodesicLoop_EdgeTraits::VertexSide::UNDEFINED;
            for (auto e = mesh->ve_iter(v); e; ++e) {
                if (StdUtil::find(E, e) != E.end()) {
                    vertexSide = mesh->data(e).geodesicLoop_vertexSide;
                    break;
                }
            }
            
            return
                vertexSide == GeodesicLoop_EdgeTraits::VertexSide::LEFT  && angles.first  >= StdUtil::pi() ||
                vertexSide == GeodesicLoop_EdgeTraits::VertexSide::RIGHT && angles.second >= StdUtil::pi();
        };
        auto is_fully_convex          = [mesh, determine_side_around_vertex, compute_incident_angles] (typename MeshBase::VHandle v, const PathPoint& point_prev, const PathPoint& point_next) -> bool {
            determine_side_around_vertex(v, point_prev, point_next);
            auto angles = compute_incident_angles(v, point_prev, point_next);
            
            return angles.first >= StdUtil::pi() && angles.second >= StdUtil::pi();
        };
        auto unfold_edge_sequence     = [mesh] (const EdgeSequence& E) -> void {
            auto e = E.begin();
            auto h = mesh->halfedge_handle(*e, 0);
            
            // ensure h is oriented from right to left
            if (mesh->data(mesh->to_vertex_handle(h)).geodesicLoop.stripSide != GeodesicLoop_VertexTraits::StripSide::LEFT)
                h = mesh->halfedge_handle(*e, 1);
            
            // fix first edge arbitrarily
            mesh->data(mesh->from_vertex_handle(h)).geodesicLoop.uv << 0, 1;
            mesh->data(mesh->to_vertex_handle  (h)).geodesicLoop.uv << 0, 0;
            
            MeshBase::VHandle endpoint[2] = { mesh->from_vertex_handle(h), mesh->to_vertex_handle(h) };
            
            while (true) {
                auto h_next = mesh->next_halfedge_handle(h);
                
                auto v0 = mesh->from_vertex_handle(h);
                auto v1 = mesh->to_vertex_handle  (h);
                auto v2 = mesh->to_vertex_handle  (h_next);
                // assuming that uv for v0 and v1 is known, compute uv for v2

                Eigen::Vector3d p0 = o2e(mesh->point(v0));
                Eigen::Vector3d p1 = o2e(mesh->point(v1));
                Eigen::Vector3d p2 = o2e(mesh->point(v2));
                Eigen::Vector3d n  = o2e(mesh->normal(mesh->face_handle(h)));
                
                Eigen::Vector3d p0p1 = p1 - p0;
                Eigen::Vector3d p0p2 = p2 - p0;
                Eigen::Vector3d p0pH = n.cross(p0p1);
                
                double a = p0p2.dot(p0p1) / p0p1.squaredNorm();
                double b = p0p2.dot(p0pH) / p0pH.squaredNorm();
                
                Eigen::Vector2d q0 = mesh->data(v0).geodesicLoop.uv;
                Eigen::Vector2d q1 = mesh->data(v1).geodesicLoop.uv;
                
                Eigen::Vector2d q0q1 = q1 - q0;
                Eigen::Vector2d q0qH = EigenUtil::rotate90(q0q1);
                
                Eigen::Vector2d q2 = q0 + a * q0q1 + b * q0qH;
                
                // store result to approrpiate place by checking if v2 is endpoint or not
                if (v2 == endpoint[0] || v2 == endpoint[1])
                    mesh->data(v2).geodesicLoop.uv_duplicate = q2;
                else
                    mesh->data(v2).geodesicLoop.uv           = q2;
                
                // terminate if reaching the end
                if (e == --E.end())
                    break;
                
                // go to next edge
                ++e;
                h = h_next;
                if (mesh->edge_handle(h) != *e)
                    h = mesh->next_halfedge_handle(h);      // ensure h corresponds to e
                h = mesh->opposite_halfedge_handle(h);
            }
        };
        auto compute_relaxed_geodesic = [mesh] (const EdgeSequence& E, typename MeshBase::VHandle v_relaxed) -> Path {
            // reference: http://digestingduck.blogspot.ch/2010/03/simple-stupid-funnel-algorithm.html
            
            // contract uv of the other endpoint (opposite side of v_relaxed) of E's first edge to uv of v_relaxed
            {
                auto v_other = OpenMeshUtil::opposite_vertex(*mesh, E[0], v_relaxed);
                auto& v_relaxed_data = mesh->data(v_relaxed).geodesicLoop;
                auto& v_other_data   = mesh->data(v_other  ).geodesicLoop;
            
                v_other_data.uv           = v_relaxed_data.uv;
                v_other_data.uv_duplicate = v_relaxed_data.uv_duplicate;
            }
            
            std::vector<MeshBase::VHandle> apices;
            apices.push_back(v_relaxed);
            
            MeshBase::VHandle apex;
            MeshBase::VHandle funnel_left;
            MeshBase::VHandle funnel_right;
            apex = funnel_left = funnel_right = v_relaxed;
            
            for (auto e = E.begin() + 2; e != --E.end(); ++e) {
                auto portal_left  = mesh->to_vertex_handle(mesh->halfedge_handle(*e, 0));
                auto portal_right = mesh->to_vertex_handle(mesh->halfedge_handle(*e, 1));
                
                // ensure correct side
                bool is_portal_flipped = false;
                if (mesh->data(portal_left).geodesicLoop.stripSide == GeodesicLoop_VertexTraits::StripSide::RIGHT) {
                    std::swap(portal_left, portal_right);
                    is_portal_flipped = true;
                }
                
                Eigen::Vector2d apex_pos         = mesh->data(apex        ).geodesicLoop.uv;
                Eigen::Vector2d funnel_left_pos  = mesh->data(funnel_left ).geodesicLoop.uv;
                Eigen::Vector2d funnel_right_pos = mesh->data(funnel_right).geodesicLoop.uv;
                Eigen::Vector2d portal_left_pos  = mesh->data(portal_left ).geodesicLoop.uv;
                Eigen::Vector2d portal_right_pos = mesh->data(portal_right).geodesicLoop.uv;
                
                if (e == --E.end()) {
                    portal_left_pos  = mesh->data(portal_left ).geodesicLoop.uv_duplicate;
                    portal_right_pos = mesh->data(portal_right).geodesicLoop.uv_duplicate;
                }
                
                if (EigenUtil::triangle_area(apex_pos, funnel_right_pos, portal_right_pos) <= 0) {
                    // update right funnel
                    if (apex == funnel_right || EigenUtil::triangle_area(apex_pos, funnel_left_pos, portal_right_pos) > 0) {
                        // tighten right funnel
                        funnel_right = portal_right;
                        
                    } else {
                        // right portal occludes left funnel -> reset apex as left funnel
                        apex = funnel_right = funnel_left;
                        
                        // add new apex
                        apices.push_back(apex);
                        
                        // TODO: restart scan!
                        
                        continue;
                    }
                }
                
                if (EigenUtil::triangle_area(apex_pos, funnel_left_pos, portal_left_pos) >= 0) {
                    // update left funnel
                    if (apex == funnel_left || EigenUtil::triangle_area(apex_pos, funnel_right_pos, portal_left_pos) < 0) {
                        // tighten left funnel
                        funnel_left = portal_left;
                        
                    } else {
                        // left portal occludes right funnel -> reset apex as right funnel
                        apex = funnel_left = funnel_right;
                        
                        // add new apex
                        apices.push_back(apex);
                        
                        // TODO: restart scan!
                        
                        continue;
                    }
                }
            }
            
            // back of apices
            apices.push_back(v_relaxed);
            
            // go through E, adding apices and edge points that apices intersect.
            Path result;
            auto apex0 = apices.begin();
            auto apex1 = apex0 + 1;
            for (auto e : E) {
                if (OpenMeshUtil::contains(*mesh, e, *apex0)) {
                    result.push_back(*apex0);
                
                } else if (OpenMeshUtil::contains(*mesh, e, *apex1)) {
                    result.push_back(*apex1);
                    ++apex0;
                    ++apex1;
                    if (apex1 == --apices.end())
                        break;
                
                } else {
                    // get intersection between (apex0, apex1) and (v0, v1) where vi is to_vertex(halfedge(e, i))
                    Eigen::Vector2d apex0_pos = mesh->data(*apex0).geodesicLoop.uv;
                    Eigen::Vector2d apex1_pos = mesh->data(*apex1).geodesicLoop.uv;
                    if (apex1 == --apices.end())
                        apex1_pos = mesh->data(*apex1).geodesicLoop.uv_duplicate;
                    
                    auto v0v1 = OpenMeshUtil::edge_to_vertex_pair(*mesh, e);
                    Eigen::Vector2d v0_pos = mesh->data(v0v1.first ).geodesicLoop.uv;
                    Eigen::Vector2d v1_pos = mesh->data(v0v1.second).geodesicLoop.uv;
                    
                    double s, t;
                    EigenUtil::intersection(apex0_pos, apex1_pos, v0_pos, v1_pos, s, t);
                    
                    // insert edge point
                    result.push_back(PathPoint(e, t));
                }
            }
            result.erase(std::unique(result.begin(), result.end()), result.end());

            return result;
        };
        auto next_edge_seuqnece       = [mesh] (const EdgeSequence& E, Path path) -> EdgeSequence {
            EdgeSequence result;
            
            for (auto point = path.begin(); point != --path.end(); ++point) {
                if (point->type() == PathPoint::Type::EDGE) {
                    // very simple:)
                    result.push_back(point->edge);
                    continue;
                }
                
                auto v = point->vertex;
                auto point_prev = point == path.begin() ? --path.end() : (point - 1);
                auto point_next = point == --path.end() ? path.begin() : (point + 1);
                
                // find an edge in E which doesn't contain v
                auto search_pos = E.begin();
                for (; search_pos != --E.end(); ++search_pos) {
                    if (!OpenMeshUtil::contains(*mesh, *search_pos, v))
                        break;
                }
                if (search_pos == --E.end())
                    // degenerate case (probably?)
                    return EdgeSequence();
                
                // advance search_pos until it finds the first edge that contains v
                while (!OpenMeshUtil::contains(*mesh, *search_pos, v)) {
                    ++search_pos;
                    if (search_pos == --E.end())
                        search_pos = E.begin();
                }
                
                // add sequence of edges that contain v
                EdgeSequence v_edges_old;
                while (OpenMeshUtil::contains(*mesh, *search_pos, v)) {
                    v_edges_old.push_back(*search_pos);
                    ++search_pos;
                    if (search_pos == --E.end())
                        search_pos = E.begin();
                }
                
                // collect ordered array of edges around v except for those coinciding with prev/next path vertex (if any)
                EdgeSequence v_edges;
                for (auto h = mesh->voh_iter(v); h; ++h) {
                    auto v2 = mesh->to_vertex_handle(h);
                    if (v2 == point_prev->vertex || v2 == point_next->vertex)
                        continue;
                    
                    v_edges.push_back(mesh->edge_handle(h));
                }
                
                // reorder the array s.t. edge sequence contained in the old one comes first
                auto rotate_pos = v_edges.begin();
                for (; rotate_pos != v_edges.end(); ++rotate_pos) {
                    auto pos_prev = rotate_pos == v_edges.begin() ? --v_edges.end() : (rotate_pos - 1);
                    if (std::find(v_edges_old.begin(), v_edges_old.end(), *pos_prev  ) == v_edges_old.end() &&
                        std::find(v_edges_old.begin(), v_edges_old.end(), *rotate_pos) != v_edges_old.end())
                    {
                        break;
                    }
                }
                std::rotate(v_edges_old.begin(), rotate_pos, v_edges_old.end());
                
                // get new edge sequence around v
                auto insert_pos = v_edges.begin();
                for (; insert_pos != v_edges.end(); ++insert_pos) {
                    if (std::find(v_edges_old.begin(), v_edges_old.end(), *insert_pos) == v_edges_old.end())
                        break;
                }
                EdgeSequence v_edges_new;
                v_edges_new.insert(v_edges_new.end(), insert_pos, v_edges.end());
                
                // reverse the order if needed
                if (v_edges.front() == v_edges_old.front())
                    std::reverse(v_edges_new.begin(), v_edges_new.end());
                
                if (v_edges_new.empty())
                    // v is on mesh boundary, cannot improve edge sequence here...
                    result.insert(result.end(), v_edges_old.begin(), v_edges_old.end());
                else
                    result.insert(result.end(), v_edges_new.begin(), v_edges_new.end());
            }
            
            // to comply with the semantics of closed sequence
            result.push_back(result.front());
            
            return result;
        };
        // =====================================================================================================================================helper functions
        
        while (true) {
            // clear side flag for all mesh vertices
            for (auto v : mesh->vertices())
                mesh->data(v).geodesicLoop.stripSide = GeodesicLoop_VertexTraits::StripSide::UNDEFINED;
            
            // set side flag for endpoints of the first edge in E
            for (int i = 0; i < 2; ++i) {
                auto h = mesh->halfedge_handle(E[0], i);
                auto h_next  = mesh->next_halfedge_handle(h);
                auto h_next2 = mesh->next_halfedge_handle(h_next);
                
                if (mesh->edge_handle(h_next) == E[1] || mesh->edge_handle(h_next2) == E[1]) {
                    // this halfedge is on the side facing E[1] -> set its to/from vertex to left/right side
                    mesh->data(mesh->to_vertex_handle  (h)).geodesicLoop.stripSide = GeodesicLoop_VertexTraits::StripSide::LEFT;
                    mesh->data(mesh->from_vertex_handle(h)).geodesicLoop.stripSide = GeodesicLoop_VertexTraits::StripSide::RIGHT;
                    break;
                }
            }
            
            // set side flag for other vertices in E
            for (auto e = ++E.begin(); e != --E.end(); ++e) {
                auto v0v1 = OpenMeshUtil::edge_to_vertex_pair(*mesh, *e);
                auto& v0_side = mesh->data(v0v1.first ).geodesicLoop.stripSide;
                auto& v1_side = mesh->data(v0v1.second).geodesicLoop.stripSide;
                
                // either v0 or v1 should be assigned side flag already
                if (v0_side == GeodesicLoop_VertexTraits::StripSide::UNDEFINED)
                    v0_side = static_cast<GeodesicLoop_VertexTraits::StripSide>((static_cast<int>(v1_side) + 1) % 2);
                else
                    v1_side = static_cast<GeodesicLoop_VertexTraits::StripSide>((static_cast<int>(v0_side) + 1) % 2);
            }
            
            // initialize cutting edge in E as its shortest
            typename MeshBase::EHandle e_cut;
            {
                double length_min = StdUtil::dbl_max();
                for (auto e : E) {
                    auto v0v1 = OpenMeshUtil::edge_to_point_pair(*mesh, e);
                    double length_v0v1 = (v0v1.first - v0v1.second).length();
                    if (length_v0v1 < length_min) {
                        length_min = length_v0v1;
                        e_cut = e;
                    }
                }
            }
            
            // initialize v_relaxed as an endpoint of e_cut
            typename MeshBase::VHandle v_relaxed = mesh->to_vertex_handle(mesh->halfedge_handle(e_cut, 0));
            
            // iterate until finding E-restricted geodesic loop
            Path path;
            while (true) {
                // bring e_cut to front
                E.pop_back();
                StdUtil::bring_front(E, e_cut);
                E.push_back(E.front());
                
                // compute 2D unfolding of E by cutting at front
                unfold_edge_sequence(E);
                
                // compute v-relaxed geodesic loop by running funnel algorithm
                path = compute_relaxed_geodesic(E, v_relaxed);
                
                // check if v_relaxed is inward convex (i.e. incident angle is greater than or equal to pi)
                if (is_inward_convex(E, v_relaxed, path[path.size() - 2], path[1]))
                    // algorithm will eventually reach here
                    break;
                
                // update v_relaxed
                auto side_v_relaxed = mesh->data(v_relaxed).geodesicLoop.stripSide;
                typename MeshBase::VHandle v_touch_opposite_side;
                typename MeshBase::VHandle v_touch_same_side;
                for (auto point = ++path.begin(); point != --path.end(); ++point) {
                    if (point->type() != PathPoint::Type::VERTEX)
                        continue;
                    
                    auto side = mesh->data(point->vertex).geodesicLoop.stripSide;
                    if (side != side_v_relaxed) {
                        v_touch_opposite_side = point->vertex;
                        break;
                        
                    } else {
                        v_touch_same_side = point->vertex;
                    }
                }
                
                if (v_touch_opposite_side.is_valid() || v_touch_same_side.is_valid()) {
                    // v-relaxed geodesic passed another vertex in E -> switch to it
                    v_relaxed = v_touch_opposite_side.is_valid() ? v_touch_opposite_side : v_touch_same_side;
                    
                    // update e_cut
                    for (auto e : E) {
                        if (OpenMeshUtil::contains(*mesh, e, v_relaxed)) {
                            e_cut = e;
                            break;
                        }
                    }
                    
                } else {
                    // otherwise, try the other endpoint of e_cut
                    v_relaxed = OpenMeshUtil::opposite_vertex(*mesh, e_cut, v_relaxed);
                }
            }
            
            // check if every vertex the path passes through is fully convex (i.e. locally shortest)
            bool is_path_fully_convex = true;
            for (auto point = ++path.begin(); point != --path.end(); ++point) {
                if (point->type() != PathPoint::Type::VERTEX)
                    continue;
                
                auto point_prev = *(point - 1);
                auto point_next = *(point + 1);
                if (!is_fully_convex(point->vertex, point_prev, point_next)) {
                    is_path_fully_convex = false;
                    break;
                }
            }
            if (is_path_fully_convex)
                return path;

            // update E
            E = next_edge_seuqnece(E, path);
            
            if (E.empty())
                // result in degenerate point
                return Path();
        }
        
        // should never reach here. algorithm is supposed to terminate within finite steps
        return Path();
    }
};

}
