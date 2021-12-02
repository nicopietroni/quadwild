#pragma once
// Remember to add TRILIBRARY to the compiler's preprocessor definitions

extern "C" {
#define ANSI_DECLARATORS
#define REAL double
#define VOID void
#include <triangle.h>
#undef ANSI_DECLARATORS
#undef REAL
#undef VOID
}

#include <vector>

namespace kt84 {
    namespace triangle_util {
        inline void init(triangulateio& data) {
            data.numberofpoints             = 0;
            data.numberofpointattributes    = 0;
            data.numberoftriangles          = 0;
            data.numberofcorners            = 0;
            data.numberoftriangleattributes = 0;
            data.numberofsegments           = 0;
            data.numberofholes              = 0;
            data.numberofregions            = 0;
            data.numberofedges              = 0;
            data.pointlist                  = nullptr;
            data.pointattributelist         = nullptr;
            data.pointmarkerlist            = nullptr;
            data.trianglelist               = nullptr;
            data.triangleattributelist      = nullptr;
            data.trianglearealist           = nullptr;
            data.neighborlist               = nullptr;
            data.segmentlist                = nullptr;
            data.segmentmarkerlist          = nullptr;
            data.holelist                   = nullptr;
            data.regionlist                 = nullptr;
            data.edgelist                   = nullptr;
            data.edgemarkerlist             = nullptr;
            data.normlist                   = nullptr;
        }
        inline void free(triangulateio& data) {
            trifree(data.pointlist            );
            trifree(data.pointattributelist   );
            trifree(data.pointmarkerlist      );
            trifree(data.trianglelist         );
            trifree(data.triangleattributelist);
            trifree(data.trianglearealist     );
            trifree(data.neighborlist         );
            trifree(data.segmentlist          );
            trifree(data.segmentmarkerlist    );
            //trifree(data.holelist             );          // Triangle never allocates memory for this field
            //trifree(data.regionlist           );          // Triangle never allocates memory for this field
            trifree(data.edgelist             );
            trifree(data.edgemarkerlist       );
            trifree(data.normlist             );
        }
        struct Format_node {
            std::vector<double> pointlist;
            std::vector<double> pointattributelist;
            std::vector<int   > pointmarkerlist;
        };
        struct Format_ele {
            int numberofcorners;
            std::vector<int   > trianglelist;
            std::vector<double> triangleattributelist;
            Format_ele() : numberofcorners() {}
        };
        struct Format_poly {
            Format_node node;
            std::vector<int> segmentlist;
            std::vector<int> segmentmarkerlist;
            std::vector<double> holelist;
            // TODO: region list
        };
        struct Format_edge {
            std::vector<int> edgelist;
            std::vector<int> edgemarkerlist;
        };
        inline void triangulate(const Format_poly& in_poly, Format_poly& out_poly, Format_ele& out_ele, Format_edge& out_edge, char* switches = "pzQ") {
            triangulateio in, out;
            init(in);
            init(out);
            
            //----+
            // IN |
            //----+
            // points
            in.numberofpoints = in_poly.node.pointlist.size() / 2;
            in.pointlist      = const_cast<double*>(&in_poly.node.pointlist[0]);
            // point attributes
            in.numberofpointattributes = in_poly.node.pointattributelist.size() / in.numberofpoints;
            if (!in_poly.node.pointattributelist.empty())
                in.pointattributelist = const_cast<double*>(&in_poly.node.pointattributelist[0]);
            // point markers
            if (!in_poly.node.pointmarkerlist.empty())
                in.pointmarkerlist = const_cast<int*>(&in_poly.node.pointmarkerlist[0]);
            // segments
            in.numberofsegments = in_poly.segmentlist.size() / 2;
            if (!in_poly.segmentlist.empty())
                in.segmentlist = const_cast<int*>(&in_poly.segmentlist[0]);
            // segment markers
            if (!in_poly.segmentmarkerlist.empty())
                in.segmentmarkerlist = const_cast<int*>(&in_poly.segmentmarkerlist[0]);
            // holes
            in.numberofholes = in_poly.holelist.size() / 2;
            if (!in_poly.holelist.empty())
                in.holelist = const_cast<double*>(&in_poly.holelist[0]);
            
            ::triangulate(switches, &in, &out, 0);
            
            //-----+
            // OUT |
            //-----+
            // points
            out_poly.node.pointlist.resize(2 * out.numberofpoints);
            for (size_t i = 0; i < out_poly.node.pointlist.size(); ++i)
                out_poly.node.pointlist[i] = out.pointlist[i];
            // point attributes
            out_poly.node.pointattributelist.resize(out.numberofpointattributes * out.numberofpoints);
            for (size_t i = 0; i < out_poly.node.pointattributelist.size(); ++i)
                out_poly.node.pointattributelist[i] = out.pointattributelist[i];
            // point markers
            if (out.pointmarkerlist != nullptr) {
                out_poly.node.pointmarkerlist.resize(out.numberofpoints);
                for (size_t i = 0; i < out_poly.node.pointmarkerlist.size(); ++i)
                    out_poly.node.pointmarkerlist[i] = out.pointmarkerlist[i];
            } else {
                out_poly.node.pointmarkerlist.clear();
            }
            // segments
            if (out.segmentmarkerlist != nullptr) {
                out_poly.segmentlist.resize(2 * out.numberofsegments);
                for (size_t i = 0; i < out_poly.segmentlist.size(); ++i)
                    out_poly.segmentlist[i] = out.segmentlist[i];
            } else {
                out_poly.segmentlist.clear();
            }
            // segment markers
            if (out.segmentmarkerlist != nullptr) {
                out_poly.segmentmarkerlist.resize(out.numberofsegments);
                for (size_t i = 0; i < out_poly.segmentmarkerlist.size(); ++i)
                    out_poly.segmentmarkerlist[i];
            } else {
                out_poly.segmentmarkerlist.clear();
            }
            // triangles
            out_ele.numberofcorners = out.numberofcorners;
            out_ele.trianglelist.resize(out.numberofcorners * out.numberoftriangles);
            for (size_t i = 0; i < out_ele.trianglelist.size(); ++i)
                out_ele.trianglelist[i] = out.trianglelist[i];
            // triangle attributes
            out_ele.triangleattributelist.resize(out.numberoftriangleattributes * out.numberoftriangles);
            for (size_t i = 0; i < out_ele.triangleattributelist.size(); ++i)
                out_ele.triangleattributelist[i] = out.triangleattributelist[i];
            // edges
            if (out.edgelist != nullptr) {
                out_edge.edgelist.resize(2 * out.numberofedges);
                for (size_t i = 0; i < out_edge.edgelist.size(); ++i)
                    out_edge.edgelist[i] = out.edgelist[i];
            } else {
                out_edge.edgelist.clear();
            }
            // edge markers
            if (out.edgemarkerlist != nullptr) {
                out_edge.edgemarkerlist.resize(out.numberofedges);
                for (size_t i = 0; i < out_edge.edgemarkerlist.size(); ++i)
                    out_edge.edgemarkerlist[i] = out.edgemarkerlist[i];
            } else {
                out_edge.edgemarkerlist.clear();
            }
            
            free(out);
        }
        inline void triangulate(const Format_node& in_node, Format_ele& out_ele, char* switches = "zQ") {
            Format_poly in_poly, out_poly;
            Format_edge out_edge;
            in_poly.node = in_node;
            triangulate(in_poly, out_poly, out_ele, out_edge, switches);
        }
        inline void triangulate(const Format_poly& in_poly, Format_poly& out_poly, Format_ele& out_ele, char* switches = "zQ") {
            Format_edge out_edge;
            triangulate(in_poly, out_poly, out_ele, out_edge, switches);
        }
    }
}
