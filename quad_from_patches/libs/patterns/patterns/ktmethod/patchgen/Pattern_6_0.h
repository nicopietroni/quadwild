#pragma once
#include "Pattern.h"
#include "ILP.h"
#include <sstream>
#include <vcg/complex/algorithms/create/platonic.h>
#include "edgeloop.h"
/*
equation for pattern 0:
  |0|     |1|     |0|     |0|     |0|     |1|   |1|    |1|   |l0|
  |1|     |0|     |1|     |0|     |0|     |0|   |1|    |0|   |l1|
p0|0| + p1|1| + p2|0| + p3|1| + p4|0| + p5|0| + |1| + x|0| = |l2|
  |0|     |0|     |1|     |0|     |1|     |0|   |1|    |1|   |l3|
  |0|     |0|     |0|     |1|     |0|     |1|   |1|    |0|   |l4|
  |1|     |0|     |0|     |0|     |1|     |0|   |1|    |0|   |l5|
*/
namespace patchgen {
    template <>
    struct Pattern<6, 0> {
        static Eigen::MatrixXd& get_constraint_matrix() {
            static Eigen::MatrixXd constraint_matrix;
            if (constraint_matrix.size() == 0) {
                constraint_matrix.resize(6, 7);
                constraint_matrix << 0, 1, 0, 0, 0, 1, 1,
                                     1, 0, 1, 0, 0, 0, 0,
                                     0, 1, 0, 1, 0, 0, 0,
                                     0, 0, 1, 0, 1, 0, 1,
                                     0, 0, 0, 1, 0, 1, 0,
                                     1, 0, 0, 0, 1, 0, 0;
            }
            return constraint_matrix;
        }
    
        static Eigen::VectorXd get_constraint_rhs(const Eigen::VectorXi& l) {
            return kt84::make_Vector6d(l[0] - 1,
                                       l[1] - 1,
                                       l[2] - 1,
                                       l[3] - 1,
                                       l[4] - 1,
                                       l[5] - 1);
        }
    
        static int& get_variable(PatchParam& param, int index) {
            if (index < 6) return param.p[index];
            return param.x;
        }

        static bool get_default_parameter(const Eigen::VectorXi& l, PatchParam& param) {
            auto& constraint_matrix = get_constraint_matrix();
            const int num_variables = constraint_matrix.cols();
            ILP ilp(num_variables);
            ilp.add_constraint(constraint_matrix, EQ, get_constraint_rhs(l));
        
            // arbitrary constraints and objective
            // xmin
            ilp.set_objective(kt84::make_Vector7d(0, 0, 0, 0, 0, 0, 1), false);
            if (!ilp.solve()) return false;
            int xmin = ilp.get_variables()[6];
            // xmax
            ilp.refresh();
            ilp.set_objective(kt84::make_Vector7d(0, 0, 0, 0, 0, 0, 1), true);
            if (!ilp.solve()) return false;
            int xmax = ilp.get_variables()[6];
            // xmid-1<=x<=xmid+1
            ilp.refresh();
            int xmid = (xmin + xmax) / 2;
            ilp.add_constraint(kt84::make_Vector7d(0, 0, 0, 0, 0, 0, 1), LE, xmid + 1);
            ilp.add_constraint(kt84::make_Vector7d(0, 0, 0, 0, 0, 0, 1), GE, xmid - 1);
            // maximize p0+p1+p2+p3+p4+p5
            ilp.set_objective(kt84::make_Vector7d(1, 1, 1, 1, 1, 1, 0), true);
            if (!ilp.solve()) return false;
        
            auto variables = ilp.get_variables();
            for (int i = 0; i < num_variables; ++i)
                get_variable(param, i) = variables[i];
        
            param.pattern_id = 0;
            return true;
        }
        template <typename PatchT>
        static void generate_subtopology(const PatchParam& param, PatchT& patch) {
            /*
        |      v--x
        |   C4----C3
        |  /       \
        | /         \
        |C5----------C2
        | \         /
        |  \       /
        |   C0----C1
        |      ^--x
            */
         /*   patch.clear();
            typename PatchT::VHandle C[6];
            for (int i = 0; i < 6; ++i) C[i] = add_tagged_vertex(patch, i, true );
        
            patch.add_face(C[0], C[1], C[2], C[5]);
            patch.add_face(C[2], C[3], C[4], C[5]);
        
            auto h_insert_x = patch.halfedge_handle(C[1]);  // corresponds to C1-C0
            for (int i = 0; i < param.x; ++i)
                insert_edgeloop(patch, h_insert_x);
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
            vcg::tri::Allocator<PatchT>::AddVertices(patch,6);
            patch.vert[0].SetS(); //C0
            patch.vert[1].SetS(); //C1
            patch.vert[2].SetS(); //C2
            patch.vert[3].SetS(); //C3
            patch.vert[4].SetS(); //C4
            patch.vert[5].SetS(); //C5

            // Populate index sides
            side_indexL[0]=1;
            side_indexR[0]=0;
            side_indexL[1]=0;
            side_indexR[1]=5;
            side_indexL[2]=5;
            side_indexR[2]=4;
            side_indexL[3]=4;
            side_indexR[3]=3;
            side_indexL[4]=3;
            side_indexR[4]=2;
            side_indexL[5]=2;
            side_indexR[5]=1;

            typename PatchT::FaceIterator pfi=vcg::tri::Allocator<PatchT>::AddFaces(patch,2);
            (*pfi).Alloc(4); //C[0], C[1], C[2], C[5]
            (*pfi).V(0)=&patch.vert[0];(*pfi).V(1)=&patch.vert[1];(*pfi).V(2)=&patch.vert[2];(*pfi).V(3)=&patch.vert[5];
            pfi++;
            (*pfi).Alloc(4); //C[2], C[3], C[4], C[5]
            (*pfi).V(0)=&patch.vert[2];(*pfi).V(1)=&patch.vert[3];(*pfi).V(2)=&patch.vert[4];(*pfi).V(3)=&patch.vert[5];

            vcg::face::Pos<typename PatchT::FaceType> startPos=getPosFromIndex(patch,0,1);
            for (int i = 0; i < param.x; ++i)
               insert_edgeloop<PatchT>(patch,startPos);

//            cout<<"patch 6 -0 "<<endl;
        }
        /*static VariableIndicators& get_variable_indicators() {
            static VariableIndicators variable_indicators;
            if (variable_indicators.empty()) {
                variable_indicators.resize(1);
                variable_indicators[0].push_back(std::make_pair(PatchVertexTag::C0, PatchVertexTag::C1));       // for x
                variable_indicators[0].push_back(std::make_pair(PatchVertexTag::C5, PatchVertexTag::C2));
                variable_indicators[0].push_back(std::make_pair(PatchVertexTag::C4, PatchVertexTag::C3));
            }
            return variable_indicators;
        }*/
        static std::string get_param_str(const PatchParam& param) {
            std::stringstream ss;
            ss << "x=" << param.x;
            return ss.str();
        }
    };
}
