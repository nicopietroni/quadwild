#include "decl.h"
#include "Pattern_all.h"
#include <kt84/util.h>

bool patchgen::switch_pattern(PatchParam& param, bool is_forward) {
    const PatchParam param_old = param;
    const int num_sides = param.get_num_sides();

    for (param.permutation.next(is_forward);; param.permutation.next(is_forward)) {
        if (!param.permutation.is_valid()) {
            // permutation reached to end -> switch to next pattern, reset permutation
            (is_forward ? kt84::util::incr_mod<int> : kt84::util::decr_mod<int>)(param.pattern_id, num_patterns[num_sides]);
            param.permutation.init(num_sides, is_forward ? 0 : (2 * num_sides - 1));
        }
        if (param.pattern_id == param_old.pattern_id && param.permutation == param_old.permutation) {
            // came back to the original param -> there's no other pattern to switch
            param = param_old;
            return false;
        }
        auto l_permuted = param.get_l_permuted();
        if (num_sides == 2) {
            if (param.pattern_id == 0 && Pattern<2, 0>::get_default_parameter(l_permuted, param)) return true;
            if (param.pattern_id == 1 && Pattern<2, 1>::get_default_parameter(l_permuted, param)) return true;
        } else if (num_sides == 3) {
            if (param.pattern_id == 0 && Pattern<3, 0>::get_default_parameter(l_permuted, param)) return true;
            if (param.pattern_id == 1 && Pattern<3, 1>::get_default_parameter(l_permuted, param)) return true;
            if (param.pattern_id == 2 && Pattern<3, 2>::get_default_parameter(l_permuted, param)) return true;
            if (param.pattern_id == 3 && Pattern<3, 3>::get_default_parameter(l_permuted, param)) return true;
        } else if (num_sides == 4) {
            if (param.pattern_id == 0 && Pattern<4, 0>::get_default_parameter(l_permuted, param)) return true;
            if (param.pattern_id == 1 && Pattern<4, 1>::get_default_parameter(l_permuted, param)) return true;
            if (param.pattern_id == 2 && Pattern<4, 2>::get_default_parameter(l_permuted, param)) return true;
            if (param.pattern_id == 3 && Pattern<4, 3>::get_default_parameter(l_permuted, param)) return true;
            if (param.pattern_id == 4 && Pattern<4, 4>::get_default_parameter(l_permuted, param)) return true;
        } else if (num_sides == 5) {
            if (param.pattern_id == 0 && Pattern<5, 0>::get_default_parameter(l_permuted, param)) return true;
            if (param.pattern_id == 1 && Pattern<5, 1>::get_default_parameter(l_permuted, param)) return true;
            if (param.pattern_id == 2 && Pattern<5, 2>::get_default_parameter(l_permuted, param)) return true;
            if (param.pattern_id == 3 && Pattern<5, 3>::get_default_parameter(l_permuted, param)) return true;
            if (param.pattern_id == 4 && Pattern<5, 4>::get_default_parameter(l_permuted, param)) return true;
        } else if (num_sides == 6) {
            if (param.pattern_id == 0 && Pattern<6, 0>::get_default_parameter(l_permuted, param)) return true;
            if (param.pattern_id == 1 && Pattern<6, 1>::get_default_parameter(l_permuted, param)) return true;
            if (param.pattern_id == 2 && Pattern<6, 2>::get_default_parameter(l_permuted, param)) return true;
            if (param.pattern_id == 3 && Pattern<6, 3>::get_default_parameter(l_permuted, param)) return true;
        }
    }
    assert(false);
    return false;
}
