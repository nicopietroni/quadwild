#pragma once
#include "Pattern.h"
#include "ILP.h"
#include <sstream>
#include <vcg/complex/algorithms/create/platonic.h>
#include "edgeloop.h"
/*
equation for pattern 1:
p0|0| + p1|2| + |2| + x|1| + y|1| = |l0|
  |2|     |0|   |2|    |1|    |1|   |l1|
*/
namespace patchgen {
    template <>
    struct Pattern<2, 1> {
        static Eigen::MatrixXd& get_constraint_matrix() {
            static Eigen::MatrixXd constraint_matrix;
            if (constraint_matrix.size() == 0) {
                constraint_matrix.resize(2, 4);
                constraint_matrix << 0, 2, 1, 1,
                                     2, 0, 1, 1;
            }
            return constraint_matrix;
        }
        static Eigen::VectorXd get_constraint_rhs(const Eigen::VectorXi& l) {
            return Eigen::Vector2d(l[0] - 2, l[1] - 2);
        }
        static int& get_variable(PatchParam& param, int index) {
            if (index < 2) return param.p[index];
            if (index == 2) return param.x;
            return param.y;
        }
        static bool get_default_parameter(const Eigen::VectorXi& l, PatchParam& param) {
            auto& constraint_matrix = get_constraint_matrix();
            const int num_variables = constraint_matrix.cols();
            ILP ilp(num_variables);
            ilp.add_constraint(constraint_matrix, EQ, get_constraint_rhs(l));
    
            // arbitrary constraints and objective
            // xmin
            ilp.set_objective(Eigen::Vector4d(0, 0, 1, 0), false);
            if (!ilp.solve()) return false;
            int xmin = ilp.get_variables()[2];
            // xmax
            ilp.refresh();
            ilp.set_objective(Eigen::Vector4d(0, 0, 1, 0), true);
            if (!ilp.solve()) return false;
            int xmax = ilp.get_variables()[2];
            // xmid-1<=x<=xmid+1
            ilp.refresh();
            int xmid = (xmin + xmax) / 2;
            ilp.add_constraint(Eigen::Vector4d(0, 0, 1, 0), LE, xmid + 1);
            ilp.add_constraint(Eigen::Vector4d(0, 0, 1, 0), GE, xmid - 1);
            // maximize p0+p1
            ilp.set_objective(Eigen::Vector4d(1, 1, 0, 0), true);
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
        |   y--v __V1__  v--x
        |     __/      \__
        |    /            \
        |   C0            C1
        |    \__        __/
        |       \__V0__/
        |   x--^       ^--y
            */
         /*   patch.clear();
            typename PatchT::VHandle C[2];
            typename PatchT::VHandle V[2];
            for (int i = 0; i < 2; ++i) C[i] = add_tagged_vertex(patch, i, true );
            for (int i = 0; i < 2; ++i) V[i] = add_tagged_vertex(patch, i, false);
        
            patch.add_face(C[0], V[0], C[1], V[1]);
        
            auto h_insert_x = patch.halfedge_handle(V[0]);  // corresponds to V0-C0
            for (int i = 0; i < param.x; ++i)
                insert_edgeloop(patch, h_insert_x);
            auto h_insert_y = patch.halfedge_handle(C[1]);  // corresponds to C1-V0
            for (int i = 0; i < param.y; ++i)
                insert_edgeloop(patch, h_insert_y);
          */
            patch.Clear();
            bool hasAttributeL = vcg::tri::HasPerVertexAttribute(patch,"LeftSide");
            bool hasAttributeR = vcg::tri::HasPerVertexAttribute(patch,"RightSide");
            if(hasAttributeL)
                vcg::tri::Allocator<PatchT>::DeletePerVertexAttribute(patch,"LeftSide");
            if(hasAttributeR)
                vcg::tri::Allocator<PatchT>::DeletePerVertexAttribute(patch,"RightSide");

            auto side_indexL=vcg::tri::Allocator<PatchT>::template GetPerVertexAttribute<int>(patch,std::string("LeftSide"));
            auto side_indexR=vcg::tri::Allocator<PatchT>::template GetPerVertexAttribute<int>(patch,std::string("RightSide"));
            vcg::tri::Allocator<PatchT>::AddVertices(patch,4);
            patch.vert[0].SetS(); //C0
            patch.vert[1].SetS(); //C1
            // Populate index sides
            side_indexL[0]=1;
            side_indexR[0]=0;
            side_indexL[1]=0;
            side_indexR[1]=1;
            side_indexL[2]=0;  //V0
            side_indexR[2]=0;
            side_indexL[3]=1;  //V1
            side_indexR[3]=1;

            typename PatchT::FaceIterator pfi=vcg::tri::Allocator<PatchT>::AddFaces(patch,1);
            (*pfi).Alloc(4); //C[0], V[0], C[1], V[1]
            (*pfi).V(0)=&patch.vert[0];(*pfi).V(1)=&patch.vert[2];(*pfi).V(2)=&patch.vert[1];(*pfi).V(3)=&patch.vert[3];

            vcg::face::Pos<typename PatchT::FaceType> startPos=getPosFromIndex(patch,0,2);
            for (int i = 0; i < param.x; ++i)
               insert_edgeloop<PatchT>(patch,startPos);

            startPos=getPosFromIndex(patch,2,1);
            for (int i = 0; i < param.y; ++i)
               insert_edgeloop<PatchT>(patch,startPos);

//            cout<<"patch 2 -- 1 "<<endl;
        }
        /*static VariableIndicators& get_variable_indicators() {
            static VariableIndicators variable_indicators;
            if (variable_indicators.empty()) {
                variable_indicators.resize(2);
                variable_indicators[0].push_back(std::make_pair(PatchVertexTag::C0, PatchVertexTag::V0));       // for x
                variable_indicators[0].push_back(std::make_pair(PatchVertexTag::C1, PatchVertexTag::V1));
                variable_indicators[1].push_back(std::make_pair(PatchVertexTag::C0, PatchVertexTag::V1));       // for y
                variable_indicators[1].push_back(std::make_pair(PatchVertexTag::C1, PatchVertexTag::V0));
            }
            return variable_indicators;

        }*/
        static std::string get_param_str(const PatchParam& param) {
            std::stringstream ss;
            ss << "x=" << param.x
               << "_y=" << param.y;
            return ss.str();
        }
    };
}
