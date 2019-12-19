#pragma once
#include "Pattern.h"

namespace patchgen {
    // generate topology according to the input pattern parameter
    template <typename PatchT>
    void generate_topology(const PatchParam& param, PatchT& patch);
    
    // compute a default pattern parameter from the input l, then generate topology accordingly
    template <typename PatchT>
    void generate_topology(const Eigen::VectorXi& l, PatchParam& param, PatchT& patch);
    
    // generate an intermediate mesh (after edge loops insertion, before padding)
    template <typename PatchT>
    void generate_subtopology(int num_sides, int pattern_id, const PatchParam& param, PatchT& patch);
    
    template <typename PatchT>
    void add_padding(PatchT& patch, const PatchParam& param);
    
    Eigen::MatrixXd& get_constraint_matrix(int num_sides, int pattern_id);
    inline int get_num_variables(int num_sides, int pattern_id) {
        return get_constraint_matrix(num_sides, pattern_id).cols();
    }
    
    Eigen::VectorXd get_constraint_rhs(int num_sides, int pattern_id, const Eigen::VectorXi& l);
    
    int& get_variable(int num_sides, int pattern_id, PatchParam& param, int index);
    
    PatchParam get_default_parameter(const Eigen::VectorXi& l);
    
    //VariableIndicators& get_variable_indicators(int num_sides, int pattern_id);
    
    std::string get_param_str(int num_sides, int pattern_id, const PatchParam& param);
    
    bool adjust_parameter(PatchParam& param, int variable_index, bool is_increase);
    
    bool switch_pattern(PatchParam& param, bool is_forward);
    
    Eigen::Vector2d get_boundary_geometry(int num_sides, double t);
}
