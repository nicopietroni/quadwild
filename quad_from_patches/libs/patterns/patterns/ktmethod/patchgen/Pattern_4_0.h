#pragma once
#include "Pattern.h"
#include "ILP.h"
#include <vcg/complex/algorithms/create/platonic.h>
/*
equation for pattern 0:
  |0|     |1|   |1|   |l0|
p0|1| + p1|0| + |1| = |l1|
  |0|     |1|   |1|   |l2|
  |1|     |0|   |1|   |l3|
*/
namespace patchgen {
    template <>
    struct Pattern<4, 0> {
        static Eigen::MatrixXd& get_constraint_matrix() {
            static Eigen::MatrixXd constraint_matrix;
            if (constraint_matrix.size() == 0) {
                constraint_matrix.resize(4, 2);
                constraint_matrix << 0, 1,
                                     1, 0,
                                     0, 1,
                                     1, 0;
            }
            return constraint_matrix;
        }
    
        static Eigen::VectorXd get_constraint_rhs(const Eigen::VectorXi& l) {
            return Eigen::Vector4d(l[0] - 1,
                                   l[1] - 1,
                                   l[2] - 1,
                                   l[3] - 1);
        }
    
        static int& get_variable(PatchParam& param, int index) {
            return param.p[index];
        }

        static bool get_default_parameter(const Eigen::VectorXi& l, PatchParam& param) {
            auto& constraint_matrix = get_constraint_matrix();
            const int num_variables = constraint_matrix.cols();
            ILP ilp(num_variables);
            ilp.add_constraint(constraint_matrix, EQ, get_constraint_rhs(l));
        
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
        |   C3----C2
        |   |      |
        |   |      |
        |   C0----C1
            */
            /*patch.Clear();
            typename PatchT::VHandle C[4];
            for (int i = 0; i < 4; ++i) C[i] = add_tagged_vertex(patch, i, true );
            patch.add_face(C[0], C[1], C[2], C[3]);*/
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
            patch.vert[2].SetS(); //C2
            patch.vert[3].SetS(); //C3

            // Populate index sides
            side_indexL[0]=1;
            side_indexR[0]=0;
            side_indexL[1]=0;
            side_indexR[1]=3;
            side_indexL[2]=3;
            side_indexR[2]=2;
            side_indexL[3]=2;
            side_indexR[3]=1;

            typename PatchT::FaceIterator pfi=vcg::tri::Allocator<PatchT>::AddFaces(patch,1);
            (*pfi).Alloc(4);
            for(int i=0;i<4;i++)
                (*pfi).V(i)=&patch.vert[i];

//            cout<<"patch 4 - 0 "<<endl;
        }
        /*static VariableIndicators& get_variable_indicators() {
            static VariableIndicators variable_indicators;
            return variable_indicators;
        }*/
        static std::string get_param_str(const PatchParam& param) { return ""; }
    };
}
