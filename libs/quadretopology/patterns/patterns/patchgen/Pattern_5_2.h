#pragma once
#include "Pattern.h"
#include "ILP.h"
#include <kt84/eigen_def.h>
#include <kt84/openmesh/edgeloop.h>
#include <sstream>

/*
equation for pattern 2:
  |0|     |1|     |0|     |0|     |1|     |0|     |1|     |1|   |4|    |2|   |l0|
  |1|     |0|     |1|     |0|     |0|     |1|     |0|     |0|   |1|    |0|   |l1|
p0|0| + p1|1| + p2|0| + p3|1| + p4|0| + q0|0| + q1|1| + q4|0| + |1| + x|0| = |l2|
  |0|     |0|     |1|     |0|     |1|     |0|     |0|     |1|   |1|    |0|   |l3|
  |1|     |0|     |0|     |1|     |0|     |1|     |0|     |0|   |1|    |0|   |l4|
*/
namespace patchgen {
    template <>
    struct Pattern<5, 2> {
        static Eigen::MatrixXd& get_constraint_matrix() {
            static Eigen::MatrixXd constraint_matrix;
            if (constraint_matrix.size() == 0) {
                constraint_matrix.resize(5, 9);
                constraint_matrix << 0, 1, 0, 0, 1, 0, 1, 1, 2,
                                     1, 0, 1, 0, 0, 1, 0, 0, 0,
                                     0, 1, 0, 1, 0, 0, 1, 0, 0,
                                     0, 0, 1, 0, 1, 0, 0, 1, 0,
                                     1, 0, 0, 1, 0, 1, 0, 0, 0;
            }
            return constraint_matrix;
        }
    
        static Eigen::VectorXd get_constraint_rhs(const Eigen::VectorXi& l) {
            return kt84::make_Vector5d(l[0] - 4,
                                       l[1] - 1,
                                       l[2] - 1,
                                       l[3] - 1,
                                       l[4] - 1);
        }
    
        static int& get_variable(PatchParam& param, int index) {
            if (index < 5) return param.p[index];
            if (index == 5) return param.q[0];
            if (index == 6) return param.q[1];
            if (index == 7) return param.q[4];
            return param.x;
        }

        static bool get_default_parameter(const Eigen::VectorXi& l, PatchParam& param) {
            auto& constraint_matrix = get_constraint_matrix();
            const int num_variables = constraint_matrix.cols();
            ILP ilp(num_variables);
            ilp.add_constraint(constraint_matrix, EQ, get_constraint_rhs(l));
        
            // arbitrary constraints and objective
            // xmin
            ilp.set_objective(kt84::make_Vector9d(0, 0, 0, 0, 0, 0, 0, 0, 1), false);
            if (!ilp.solve()) return false;
            int xmin = ilp.get_variables()[8];
            // xmax
            ilp.refresh();
            ilp.set_objective(kt84::make_Vector9d(0, 0, 0, 0, 0, 0, 0, 0, 1), true);
            if (!ilp.solve()) return false;
            int xmax = ilp.get_variables()[8];
            // xmid-1<=x<=xmid+1
            ilp.refresh();
            int xmid = (xmin + xmax) / 2;
            ilp.add_constraint(kt84::make_Vector9d(0, 0, 0, 0, 0, 0, 0, 0, 1), LE, xmid + 1);
            ilp.add_constraint(kt84::make_Vector9d(0, 0, 0, 0, 0, 0, 0, 0, 1), GE, xmid - 1);
            // p0<=q0
            ilp.add_constraint(kt84::make_Vector9d(1, 0, 0, 0, 0, -1, 0, 0, 0), LE, 0);
            // p1<=q1
            ilp.add_constraint(kt84::make_Vector9d(0, 1, 0, 0, 0, 0, -1, 0, 0), LE, 0);
            // p4<=q4
            ilp.add_constraint(kt84::make_Vector9d(0, 0, 0, 0, 1, 0, 0, -1, 0), LE, 0);
            // maximize p0+p1+p2+p3+p4
            ilp.set_objective(kt84::make_Vector9d(1, 1, 1, 1, 1, 0, 0, 0, 0), true);
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
        |                         _C3
        |                       _/ | \_
        |                     _/   |   \_
        |                   _/     |     \_
        |                 _/       |       \_
        |        q4-->  _/         |         \_   <--q1
        |             _/           |           \_
        |           _/             |             \_
        |         _/               V4              \_
        |       _/                /|\                \_
        |     _/                 / | \                 \_
        |   C4                  /  |  \                  C2
        |    \                 /   |   \                 /
        |     \               /    |    \               /
        |      \             /    V3     \             /
        |       \           /     /\      \           /
        |        \         /     /  \      \         /
        |   q0--> \       /     /    \      \       /  <--q0
        |          \     /     /      \      \     /
        |           \   /     /        \      \   /
        |            \ /     /          \      \ /
        |            C0----V0-----V1---V2-------C1
        |            x--^  q1--^     ^--q4   ^--x
            */
            patch.clear();
            typename PatchT::VHandle C[5];
            typename PatchT::VHandle V[5];
            for (int i = 0; i < 5; ++i) C[i] = add_tagged_vertex(patch, i, true );
            for (int i = 0; i < 5; ++i) V[i] = add_tagged_vertex(patch, i, false);
        
            patch.add_face(C[0], V[0], V[3], V[4]);
            patch.add_face(C[1], V[4], V[3], V[2]);
            patch.add_face(V[0], V[1], V[2], V[3]);
            patch.add_face(V[4], C[3], C[4], C[0]);
            patch.add_face(V[4], C[1], C[2], C[3]);
        
            auto h_insert_x = patch.halfedge_handle(V[0]);  // corresponds to V0-C0
            for (int i = 0; i < param.x; ++i)
                insert_edgeloop(patch, h_insert_x);
            auto h_insert_q0 = patch.halfedge_handle(C[0]);  // corresponds to C0-C4
            for (int i = 0; i < param.q[0]; ++i)
                insert_edgeloop(patch, h_insert_q0);
            auto h_insert_q1 = patch.halfedge_handle(V[1]);  // corresponds to V1-V0
            for (int i = 0; i < param.q[1]; ++i)
                insert_edgeloop(patch, h_insert_q1);
            auto h_insert_q4 = patch.halfedge_handle(V[2]);  // corresponds to V2-V1
            for (int i = 0; i < param.q[4]; ++i)
                insert_edgeloop(patch, h_insert_q4);
        }
        static VariableIndicators& get_variable_indicators() {
            static VariableIndicators variable_indicators;
            if (variable_indicators.empty()) {
                variable_indicators.resize(4);
                variable_indicators[0].push_back(std::make_pair(PatchVertexTag::C0, PatchVertexTag::C4));       // for q0
                variable_indicators[0].push_back(std::make_pair(PatchVertexTag::C1, PatchVertexTag::C2));
                variable_indicators[0].push_back(std::make_pair(PatchVertexTag::C3, PatchVertexTag::V4));
                variable_indicators[1].push_back(std::make_pair(PatchVertexTag::V0, PatchVertexTag::V1));       // for q1
                variable_indicators[1].push_back(std::make_pair(PatchVertexTag::V2, PatchVertexTag::V3));
                variable_indicators[1].push_back(std::make_pair(PatchVertexTag::C1, PatchVertexTag::V4));
                variable_indicators[1].push_back(std::make_pair(PatchVertexTag::C2, PatchVertexTag::C3));
                variable_indicators[2].push_back(std::make_pair(PatchVertexTag::V2, PatchVertexTag::V1));       // for q4
                variable_indicators[2].push_back(std::make_pair(PatchVertexTag::V0, PatchVertexTag::V3));
                variable_indicators[2].push_back(std::make_pair(PatchVertexTag::C0, PatchVertexTag::V4));
                variable_indicators[2].push_back(std::make_pair(PatchVertexTag::C4, PatchVertexTag::C3));
                variable_indicators[3].push_back(std::make_pair(PatchVertexTag::C0, PatchVertexTag::V0));       // for x
                variable_indicators[3].push_back(std::make_pair(PatchVertexTag::C1, PatchVertexTag::V2));
                variable_indicators[3].push_back(std::make_pair(PatchVertexTag::V3, PatchVertexTag::V4));
            }
            return variable_indicators;
        }
        static std::string get_param_str(const PatchParam& param) {
            std::stringstream ss;
            ss << "q0=" << param.q[0] 
               << "_q1=" << param.q[1]
               << "_q4=" << param.q[4]
               << "_x=" << param.x;
            return ss.str();
        }
    };
}
