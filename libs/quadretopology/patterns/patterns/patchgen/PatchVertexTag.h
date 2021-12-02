#pragma once
#include <kt84/util.h>

namespace patchgen {

enum struct PatchVertexTag {
    None = -1,
    C0,
    C1,
    C2,
    C3,
    C4,
    C5,
    V0,
    V1,
    V2,
    V3,
    V4,
    V5,
    V6,
    V7,
    V8,
};

template <typename PatchT>
typename PatchT::VHandle add_tagged_vertex (PatchT& patch, int index, bool is_corner) {
    auto v = patch.add_vertex(typename PatchT::Point());
    auto& vdata = patch.data(v);
    vdata.patchgen.corner_index = is_corner ? index : -1;
    vdata.patchgen.tag          = kt84::util::add_enum(is_corner ? PatchVertexTag::C0 : PatchVertexTag::V0, index);
    return v;
}

template <typename PatchT>
typename PatchT::VHandle find_tagged_vertex(const PatchT& patch, const PatchVertexTag& tag) {
    for (auto v : patch.vertices())
        if (patch.data(v).patchgen.tag == tag) return v;
    return typename PatchT::VHandle();
}

}
