#pragma once
#include "Pattern.h"
#include "ILP.h"
#include <kt84/eigen_def.h>
#include <kt84/openmesh/edgeloop.h>
#include <sstream>

/*
equation for pattern 2:
  |0|     |1|     |0|     |1|   |3|    |1|    |1|   |l0|
p0|1| + p1|0| + p2|1| + p3|0| + |1| + x|1| + y|0| = |l1|
  |0|     |1|     |0|     |1|   |1|    |0|    |0|   |l2|
  |1|     |0|     |1|     |0|   |1|    |0|    |1|   |l3|
*/
namespace patchgen {
    template <>
    struct Pattern<4, 2> {
        static Eigen::MatrixXd& get_constraint_matrix() {
            static Eigen::MatrixXd constraint_matrix;
            if (constraint_matrix.size() == 0) {
                constraint_matrix.resize(4, 6);
                constraint_matrix << 0, 1, 0, 1, 1, 1,
                                     1, 0, 1, 0, 1, 0,
                                     0, 1, 0, 1, 0, 0,
                                     1, 0, 1, 0, 0, 1;
            }
            return constraint_matrix;
        }
    
        static Eigen::VectorXd get_constraint_rhs(const Eigen::VectorXi& l) {
            return Eigen::Vector4d(l[0] - 3,
                                   l[1] - 1,
                                   l[2] - 1,
                                   l[3] - 1);
        }
    
        static int& get_variable(PatchParam& param, int index) {
            if (index < 4) return param.p[index];
            if (index == 4) return param.x;
            return param.y;
        }

        static bool get_default_parameter(const Eigen::VectorXi& l, PatchParam& param) {
            auto& constraint_matrix = get_constraint_matrix();
            const int num_variables = constraint_matrix.cols();
            ILP ilp(num_variables);
            ilp.add_constraint(constraint_matrix, EQ, get_constraint_rhs(l));
        
            // arbitrary constraints and objective
            // p0 <= p2
            ilp.add_constraint(kt84::make_Vector6d(1, 0, -1, 0, 0, 0), LE, 0);
            // p1 <= p3
            ilp.add_constraint(kt84::make_Vector6d(0, 1, 0, -1, 0, 0), LE, 0);
            // maximize p0+p1
            ilp.set_objective(kt84::make_Vector6d(1, 1, 0, 0, 0, 0), true);
            if (!ilp.solve()) return false;
        
            auto variables = ilp.get_variables();
            for (int i = 0; i < num_variables; ++i)
                get_variable(param, i) = variables[i];
        
            param.pattern_id = 2;
            return true;
        }
        template <typename PatchT>
        static void generate_subtopology(const PatchParam& param, PatchT& patch) {
            /*
        |        C3-------------C2
        |        |            _/|
        |        |          _/  |
        |   y--> |        _/    | <--x
        |        |      _/      |
        |        |     /        |
        |        C0---V0---V1---C1
        |             x--^    ^--y
            */
            patch.clear();
            typename PatchT::VHandle C[4];
            typename PatchT::VHandle V[2];
            for (int i = 0; i < 4; ++i) C[i] = add_tagged_vertex(patch, i, true );
            for (int i = 0; i < 2; ++i) V[i] = add_tagged_vertex(patch, i, false);
            patch.add_face(C[0], V[0], C[2], C[3]);
            patch.add_face(V[0], V[1], C[1], C[2]);
        
            auto h_insert_x = patch.halfedge_handle(V[1]);  // corresponds to V1-V0
            for (int i = 0; i < param.x; ++i)
                insert_edgeloop(patch, h_insert_x);
            auto h_insert_y = patch.halfedge_handle(C[1]);  // corresponds to C1-V1
            for (int i = 0; i < param.y; ++i)
                insert_edgeloop(patch, h_insert_y);
        }
        static VariableIndicators& get_variable_indicators() {
            static VariableIndicators variable_indicators;
            if (variable_indicators.empty()) {
                variable_indicators.resize(2);
                //variable_indicators[0].push_back(std::make_pair(PatchVertexTag::V0, PatchVertexTag::V1));       // for x (not editable anyway)
                //variable_indicators[0].push_back(std::make_pair(PatchVertexTag::C1, PatchVertexTag::C2));
                //variable_indicators[1].push_back(std::make_pair(PatchVertexTag::C0, PatchVertexTag::C3));       // for y (not editable anyway)
                //variable_indicators[1].push_back(std::make_pair(PatchVertexTag::C1, PatchVertexTag::V1));
            }
            return variable_indicators;
        }
        static std::string get_param_str(const PatchParam& param) {
            std::stringstream ss;
            ss << "p0=" << param.p[0] 
               << "_p1=" << param.p[1];
            return ss.str();
        }
    };
}
