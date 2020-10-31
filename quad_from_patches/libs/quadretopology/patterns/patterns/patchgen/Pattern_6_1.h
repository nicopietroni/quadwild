#pragma once
#include "Pattern.h"
#include "ILP.h"
#include <kt84/eigen_def.h>
#include <kt84/openmesh/edgeloop.h>
#include <sstream>

/*
  |0|     |1|     |0|     |0|     |0|     |1|   |2|    |1|    |1|    |0|    |0|   |l0|
  |1|     |0|     |1|     |0|     |0|     |0|   |2|    |1|    |0|    |1|    |0|   |l1|
p0|0| + p1|1| + p2|0| + p3|1| + p4|0| + p5|0| + |1| + x|0| + y|0| + z|0| + w|1| = |l2|
  |0|     |0|     |1|     |0|     |1|     |0|   |1|    |0|    |1|    |0|    |0|   |l3|
  |0|     |0|     |0|     |1|     |0|     |1|   |1|    |0|    |0|    |1|    |0|   |l4|
  |1|     |0|     |0|     |0|     |1|     |0|   |1|    |0|    |0|    |0|    |1|   |l5|
*/
namespace patchgen {
    template <>
    struct Pattern<6, 1> {
        static Eigen::MatrixXd& get_constraint_matrix() {
            static Eigen::MatrixXd constraint_matrix;
            if (constraint_matrix.size() == 0) {
                constraint_matrix.resize(6, 10);
                constraint_matrix << 0, 1, 0, 0, 0, 1, 1, 1, 0, 0,
                                     1, 0, 1, 0, 0, 0, 1, 0, 1, 0,
                                     0, 1, 0, 1, 0, 0, 0, 0, 0, 1,
                                     0, 0, 1, 0, 1, 0, 0, 1, 0, 0,
                                     0, 0, 0, 1, 0, 1, 0, 0, 1, 0,
                                     1, 0, 0, 0, 1, 0, 0, 0, 0, 1;
            }
            return constraint_matrix;
        }
    
        static Eigen::VectorXd get_constraint_rhs(const Eigen::VectorXi& l) {
            return kt84::make_Vector6d(l[0] - 2,
                                       l[1] - 2,
                                       l[2] - 1,
                                       l[3] - 1,
                                       l[4] - 1,
                                       l[5] - 1);
        }
    
        static int& get_variable(PatchParam& param, int index) {
            if (index < 6) return param.p[index];
            if (index == 6) return param.x;
            if (index == 7) return param.y;
            if (index == 8) return param.z;
            return param.w;
        }

        static bool get_default_parameter(const Eigen::VectorXi& l, PatchParam& param) {
            auto& constraint_matrix = get_constraint_matrix();
            const int num_variables = constraint_matrix.cols();
            ILP ilp(num_variables);
            ilp.add_constraint(constraint_matrix, EQ, get_constraint_rhs(l));
        
            // arbitrary constraints and objective
            // xmin
            ilp.set_objective(kt84::make_Vector10d(0, 0, 0, 0, 0, 0, 1, 0, 0, 0), false);
            if (!ilp.solve()) return false;
            int xmin = ilp.get_variables()[6];
            // xmax
            ilp.refresh();
            ilp.set_objective(kt84::make_Vector10d(0, 0, 0, 0, 0, 0, 1, 0, 0, 0), true);
            if (!ilp.solve()) return false;
            int xmax = ilp.get_variables()[6];
            // xmid-1<=x<=xmid+1
            ilp.refresh();
            int xmid = (xmin + xmax) / 2;
            ilp.add_constraint(kt84::make_Vector10d(0, 0, 0, 0, 0, 0, 1, 0, 0, 0), LE, xmid + 1);
            ilp.add_constraint(kt84::make_Vector10d(0, 0, 0, 0, 0, 0, 1, 0, 0, 0), GE, xmid - 1);
            // maximize p0+p1+p2+p3+p4+p5
            ilp.set_objective(kt84::make_Vector10d(1, 1, 1, 1, 1, 1, 0, 0, 0, 0), true);
            if (!ilp.solve()) return false;
        
            auto variables = ilp.get_variables();
            for (int i = 0; i < num_variables; ++i)
                get_variable(param, i) = variables[i];
        
            param.pattern_id = 1;
            return true;
        }
        template <typename PatchT>
        static void generate_subtopology(const PatchParam& param, PatchT& patch) {
            /*
        |                v--y
        |           C4--------C3
        |          /  \        \
        |   z-->  /    \        \  <--w
        |        /      \        \
        |       C5      V3-------C2
        |        \     /  \      /   <--x
        |   w-->  \   /   V2----V1
        |          \ /    /    /   <--z
        |           C0---V0---C1
        |           x--^    ^--y
            */
            patch.clear();
            typename PatchT::VHandle C[6];
            typename PatchT::VHandle V[4];
            for (int i = 0; i < 6; ++i) C[i] = add_tagged_vertex(patch, i, true );
            for (int i = 0; i < 4; ++i) V[i] = add_tagged_vertex(patch, i, false);
        
            patch.add_face(C[0], V[0], V[2], V[3]);
            patch.add_face(C[1], V[1], V[2], V[0]);
            patch.add_face(V[1], C[2], V[3], V[2]);
            patch.add_face(C[2], C[3], C[4], V[3]);
            patch.add_face(C[4], C[5], C[0], V[3]);
        
            auto h_insert_x = patch.halfedge_handle(V[0]);  // corresponds to V0-C0
            for (int i = 0; i < param.x; ++i)
                insert_edgeloop(patch, h_insert_x);
            auto h_insert_y = patch.halfedge_handle(C[1]);  // corresponds to C1-V0
            for (int i = 0; i < param.y; ++i)
                insert_edgeloop(patch, h_insert_y);
            auto h_insert_z = patch.halfedge_handle(V[1]);  // corresponds to V1-C1
            for (int i = 0; i < param.z; ++i)
                insert_edgeloop(patch, h_insert_z);
            auto h_insert_w = patch.halfedge_handle(C[3]);  // corresponds to C3-C2
            for (int i = 0; i < param.w; ++i)
                insert_edgeloop(patch, h_insert_w);
        }
        static VariableIndicators& get_variable_indicators() {
            static VariableIndicators variable_indicators;
            if (variable_indicators.empty()) {
                variable_indicators.resize(4);
                variable_indicators[0].push_back(std::make_pair(PatchVertexTag::C0, PatchVertexTag::V0));       // for x
                variable_indicators[0].push_back(std::make_pair(PatchVertexTag::V3, PatchVertexTag::V2));
                variable_indicators[0].push_back(std::make_pair(PatchVertexTag::C2, PatchVertexTag::V1));
                variable_indicators[1].push_back(std::make_pair(PatchVertexTag::C1, PatchVertexTag::V0));       // for y
                variable_indicators[1].push_back(std::make_pair(PatchVertexTag::V1, PatchVertexTag::V2));
                variable_indicators[1].push_back(std::make_pair(PatchVertexTag::C2, PatchVertexTag::V3));
                variable_indicators[1].push_back(std::make_pair(PatchVertexTag::C3, PatchVertexTag::C4));
                variable_indicators[2].push_back(std::make_pair(PatchVertexTag::C1, PatchVertexTag::V1));       // for z
                variable_indicators[2].push_back(std::make_pair(PatchVertexTag::V0, PatchVertexTag::V2));
                variable_indicators[2].push_back(std::make_pair(PatchVertexTag::C0, PatchVertexTag::V3));
                variable_indicators[2].push_back(std::make_pair(PatchVertexTag::C5, PatchVertexTag::C4));
                variable_indicators[3].push_back(std::make_pair(PatchVertexTag::C0, PatchVertexTag::C5));       // for w
                variable_indicators[3].push_back(std::make_pair(PatchVertexTag::V3, PatchVertexTag::C4));
                variable_indicators[3].push_back(std::make_pair(PatchVertexTag::C2, PatchVertexTag::C3));
            }
            return variable_indicators;
        }
        static std::string get_param_str(const PatchParam& param) {
            std::stringstream ss;
            ss << "x=" << param.x
               << "_y=" << param.y
               << "_z=" << param.z
               << "_w=" << param.w;
            return ss.str();
        }
    };
}
