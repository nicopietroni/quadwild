#pragma once
#include "decl.h"
#include "Pattern_all.h"

template <typename PatchT>
void patchgen::generate_subtopology(int num_sides, int pattern_id, const PatchParam& param, PatchT& patch) {
    typedef void (*func_type)(const PatchParam& param, PatchT& patch);
    static const func_type func_table_2[2] = {
        Pattern<2, 0>::generate_subtopology,
        Pattern<2, 1>::generate_subtopology
    };
    static const func_type func_table_3[4] = {
        Pattern<3, 0>::generate_subtopology,
        Pattern<3, 1>::generate_subtopology,
        Pattern<3, 2>::generate_subtopology,
        Pattern<3, 3>::generate_subtopology,
    };
    static const func_type func_table_4[5] = {
        Pattern<4, 0>::generate_subtopology,
        Pattern<4, 1>::generate_subtopology,
        Pattern<4, 2>::generate_subtopology,
        Pattern<4, 3>::generate_subtopology,
        Pattern<4, 4>::generate_subtopology
    };
    static const func_type func_table_5[5] = {
        Pattern<5, 0>::generate_subtopology,
        Pattern<5, 1>::generate_subtopology,
        Pattern<5, 2>::generate_subtopology,
        Pattern<5, 3>::generate_subtopology,
        Pattern<5, 4>::generate_subtopology,
    };
    static const func_type func_table_6[4] = {
        Pattern<6, 0>::generate_subtopology,
        Pattern<6, 1>::generate_subtopology,
        Pattern<6, 2>::generate_subtopology,
        Pattern<6, 3>::generate_subtopology,
    };
    static const func_type* func_table[7] = {
        nullptr,
        nullptr,
        func_table_2,
        func_table_3,
        func_table_4,
        func_table_5,
        func_table_6,
    };
    func_table[num_sides][pattern_id](param, patch);
}
