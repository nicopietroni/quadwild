#include "decl.h"
#include "Pattern_all.h"

namespace {
    template <int NumSides, int PatternID>
    bool get_default_parameter_sub(patchgen::PatchParam& param) {
        for (param.permutation.init(NumSides); param.permutation.is_valid(); param.permutation.next())
            if (patchgen::Pattern<NumSides, PatternID>::get_default_parameter(param.get_l_permuted(), param))
                return true;
        return false;
    }
}

patchgen::PatchParam patchgen::get_default_parameter(const Eigen::VectorXi& l) {
    int num_sides = l.size();

    PatchParam param;
    param.l = l;
    param.p = param.q = Eigen::VectorXi::Constant(num_sides, -1);
    if (num_sides == 2) {
        if (get_default_parameter_sub<2, 0>(param)) return param;
        if (get_default_parameter_sub<2, 1>(param)) return param;
    } else if (num_sides == 3) {
        if (get_default_parameter_sub<3, 0>(param)) return param;
        if (get_default_parameter_sub<3, 1>(param)) return param;
        if (get_default_parameter_sub<3, 2>(param)) return param;
    } else if (num_sides == 4) {
        if (get_default_parameter_sub<4, 0>(param)) return param;
        if (get_default_parameter_sub<4, 1>(param)) return param;
        if (get_default_parameter_sub<4, 2>(param)) return param;
        if (get_default_parameter_sub<4, 3>(param)) return param;
        if (get_default_parameter_sub<4, 4>(param)) return param;
    } else if (num_sides == 5) {
        if (get_default_parameter_sub<5, 0>(param)) return param;
        if (get_default_parameter_sub<5, 1>(param)) return param;
        if (get_default_parameter_sub<5, 2>(param)) return param;
        if (get_default_parameter_sub<5, 3>(param)) return param;
    } else if (num_sides == 6) {
        if (get_default_parameter_sub<6, 0>(param)) return param;
        if (get_default_parameter_sub<6, 1>(param)) return param;
        if (get_default_parameter_sub<6, 2>(param)) return param;
        if (get_default_parameter_sub<6, 3>(param)) return param;
    }
    
    assert(false);
    return {};
}
