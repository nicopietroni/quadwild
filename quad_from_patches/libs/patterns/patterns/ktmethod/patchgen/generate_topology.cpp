#include "Pattern.h"
#include "decl.h"
//#include "add_padding.h"
#include "generate_subtopology.h"
//#include <kt84/util.h>
//#include <kt84/openmesh/flip_faces.h>

template <typename PatchT>
void patchgen::generate_topology(const PatchParam& param, PatchT& patch) {
    int num_sides = param.get_num_sides();
    std::cout << "pattern=" << param.pattern_id;
    std::string param_str = get_param_str(num_sides, param.pattern_id, param);
    if (!param_str.empty())
        std::cout << ", param=" << param_str;
    std::cout << "\n";

    // generate topology
    generate_subtopology(num_sides, param.pattern_id, param, patch);
        
    //add_padding(patch, param);

    // reordering corner index, flip faces
    /*for (auto v : patch.vertices()) {
        int& corner_index = patch.data(v).patchgen.corner_index;
        if (corner_index == -1) continue;
        corner_index = param.permutation[corner_index];
        if (param.permutation.is_flipped()) corner_index = (corner_index + 1) % num_sides;      // Be careful!
    }
    if (param.permutation.is_flipped())
        kt84::flip_faces(patch);*/
}

template <typename PatchT>
void patchgen::generate_topology(const Eigen::VectorXi& l, PatchParam& param, PatchT& patch) {
    int num_sides = l.size();
    std::cout << "input=[";
    for (int i = 0; i < num_sides; ++i)
        std::cout << " " << l[i];
    std::cout << " ]\n";

    param = get_default_parameter(l);
    std::cout << "\n";

    generate_topology(param, patch);
}
