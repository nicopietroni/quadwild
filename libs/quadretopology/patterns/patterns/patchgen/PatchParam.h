#pragma once
#include "Permutation.h"

namespace patchgen {
    struct PatchParam {
        int pattern_id;
        Eigen::VectorXi l;              // original input without permutation
        Permutation permutation;
        Eigen::VectorXi p;              // padding
        Eigen::VectorXi q;              // padding inserted inside the pattern (p[i] and q[i] have the same effect on the boundary. i.e. same basis vector)
        int x;                          // number of edgeflows
        int y;
        int z;                          // parameters describing auxiliary edgeflow counts (only relevant for 6-sided) TODO: better explanation...
        int w;
        PatchParam()
            : pattern_id(-1)
            , x(-1)
            , y(-1)
            , z(-1)
            , w(-1)
        {}
        int get_num_sides() const { return l.size(); }
        Eigen::VectorXi get_l_permuted() const { return permutation(l); }
    };
}