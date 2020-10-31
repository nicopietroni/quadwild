#pragma once
#include "decl.h"
#include <kt84/openmesh/append_quad_strip.h>

template <typename PatchT>
void patchgen::add_padding(PatchT& patch, const PatchParam& param) {
    const int num_sides = param.get_num_sides();
    // IMPORTANT NOTE: boundary halfedges are oriented CLOCKWISE!
    typename PatchT::HHandle h_boundary_front;
    for (auto v : patch.vertices()) {
        if (patch.data(v).patchgen.corner_index == 0) {
            h_boundary_front = patch.halfedge_handle(v);
            break;
        }
    }

    for (int i = num_sides - 1; i >= 0; --i) {
        auto h_boundary_back = h_boundary_front;
        while (patch.data(patch.to_vertex_handle(h_boundary_back)).patchgen.corner_index == -1)
            h_boundary_back = patch.next_halfedge_handle(h_boundary_back);

        for (int j = 0; j < param.p[i]; ++j) {
            // clear corner flag for current corner vertices
            patch.data(patch.from_vertex_handle(h_boundary_front)).patchgen.corner_index = -1;
            patch.data(patch.to_vertex_handle  (h_boundary_back )).patchgen.corner_index = -1;

            kt84::append_quad_strip(patch, h_boundary_front, h_boundary_back);

            // set corner flag for new corner vertices
            patch.data(patch.from_vertex_handle(h_boundary_front)).patchgen.corner_index = (i + 1) % num_sides;
            patch.data(patch.to_vertex_handle  (h_boundary_back )).patchgen.corner_index = i;
        }

        // go to next (clockwise) side
        h_boundary_front = patch.next_halfedge_handle(h_boundary_back);
    }
    assert(patch.data(patch.from_vertex_handle(h_boundary_front)).patchgen.corner_index == 0);
}
