#pragma once
#include "PatchParam.h"
#include "PatchVertexTag.h"
#include <vector>
#include <string>
#include <utility>

namespace patchgen {
    typedef std::vector<std::vector<std::pair<PatchVertexTag, PatchVertexTag>>> VariableIndicators;
    
    template <int NumSides, int PatternID>
    struct Pattern {
        static Eigen::MatrixXd& get_constraint_matrix();
        
        static Eigen::VectorXd get_constraint_rhs(const Eigen::VectorXi& l);
        
        static int& get_variable(PatchParam& param, int index);
        
        static bool get_default_parameter(const Eigen::VectorXi& l, PatchParam& param);
        
        template <typename PatchT>
        static void generate_subtopology(const PatchParam& param, PatchT& patch);
        
        static VariableIndicators& get_variable_indicators();
        
        static std::string get_param_str(const PatchParam& param);
    };
}
