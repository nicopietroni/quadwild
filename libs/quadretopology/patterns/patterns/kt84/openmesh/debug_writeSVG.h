#pragma once
#include <cstdio>
#include <sstream>
#include "simple_svg_1.0.0.hpp"

namespace {
#ifdef _WIN32
    inline FILE* popen(const char* command, const char* mode) { return _popen(command, mode); }
    inline int pclose(FILE* file) { return _pclose(file); }
#endif
}
namespace kt84 {

template <class TMesh>
inline void debug_writeSVG(
    const TMesh& mesh,
    const char* filename          = "debug.svg",
    bool export_pdf               = false,
    double      dimension         = 500.0,
    double      stroke_width      = 1.0,
    bool        draw_vidx         = true,
    bool        draw_fidx         = true,
    bool        draw_eidx         = true,
    const char* font_family       = "Courier",
    double      font_size         = 10.0,
    double      font_stroke_width = 1.0)
{
    auto p0 = mesh.point(mesh.vertex_handle(0));
    double xmin = p0[0], xmax = xmin;
    double ymin = p0[1], ymax = ymin;
    for (auto v : mesh.vertices()) {
        if (mesh.has_vertex_status() && mesh.status(v).deleted()) continue;
        auto p = mesh.point(v);
        xmin = std::min<double>(xmin, p[0]);
        xmax = std::max<double>(xmax, p[0]);
        ymin = std::min<double>(ymin, p[1]);
        ymax = std::max<double>(ymax, p[1]);
    }
    double margin = (xmax - xmin + ymax - ymin) * 0.1;
    xmin -= margin;     xmax += margin;
    ymin -= margin;     ymax += margin;
    double scaling = dimension / std::min<double>(xmax - xmin, ymax - ymin);
        
    svg::Document doc(filename, svg::Layout(svg::Dimensions(scaling * (xmax - xmin), scaling * (ymax - ymin))));

    auto make_svgPoint = [&] (const typename TMesh::Point& p) {
        return svg::Point(scaling * (p[0] - xmin),
                            scaling * (p[1] - ymin));
    };
    auto int2str = [](int i) {
        std::stringstream ss;
        ss << i;
        return ss.str();
    };
        
    for (auto f : mesh.faces()) {
        if (mesh.has_face_status() && mesh.status(f).deleted()) continue;
        svg::Polygon polygon(svg::Fill(svg::Color::Silver), svg::Stroke(font_stroke_width, svg::Color::Black));
        for (auto fv = mesh.cfv_iter(f); fv.is_valid(); ++fv)
            polygon << make_svgPoint(mesh.point(*fv));
        doc << polygon;
    }
        
    svg::Font font(font_size, font_family);
    if (draw_vidx) {
        for (auto v : mesh.vertices()) {
            if (mesh.has_vertex_status() && mesh.status(v).deleted()) continue;
            doc << svg::Text(make_svgPoint(mesh.point(v)), int2str(v.idx()), svg::Fill(), font, svg::Stroke(font_stroke_width, svg::Color::Red));
        }
    }
    if (draw_fidx) {
        for (auto f : mesh.faces()) {
            if (mesh.has_face_status() && mesh.status(f).deleted()) continue;
            typename TMesh::Point p;
            p.vectorize(0);
            int cnt = 0;
            for (auto fv = mesh.cfv_iter(f); fv.is_valid(); ++fv) {
                p += mesh.point(*fv);
                ++cnt;
            }
            p *= 1.0 / cnt;
            doc << svg::Text(make_svgPoint(p), int2str(f.idx()), svg::Fill(), font, svg::Stroke(font_stroke_width, svg::Color::Blue));
        }
    }
    if (draw_eidx) {
        for (auto e : mesh.edges()) {
            if (mesh.has_edge_status() && mesh.status(e).deleted()) continue;
            auto p0 = mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(e, 0)));
            auto p1 = mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(e, 1)));
            doc << svg::Text(make_svgPoint((p0 + p1) * 0.5), int2str(e.idx()), svg::Fill(), font, svg::Stroke(1, svg::Color::Brown));
        }
    }
    doc.save();
    if (export_pdf) {
        std::stringstream command;
        command << "inkscape --export-pdf=" << filename << ".pdf " << filename;
        FILE* p = popen(command.str().c_str(), "r");
        const int BUF_SIZE = 4096;
        char buf[BUF_SIZE];
        buf[0] = 0;
        while (fgets(buf, BUF_SIZE, p) != nullptr) {}
        pclose(p);
    }
}

}
