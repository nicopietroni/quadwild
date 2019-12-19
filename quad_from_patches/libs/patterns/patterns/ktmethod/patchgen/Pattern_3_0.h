#pragma once
#include "Pattern.h"
#include "ILP.h"
#include <vcg/complex/algorithms/create/platonic.h>
/*
equation for pattern 0:
  |0|     |1|     |1|   |2|   |l0|
p0|1| + p1|0| + p2|1| + |1| = |l1|
  |1|     |1|     |0|   |1|   |l2|
*/
namespace patchgen {
    template <>
    struct Pattern<3, 0> {
        static Eigen::MatrixXd& get_constraint_matrix() {
            static Eigen::MatrixXd constraint_matrix;
            if (constraint_matrix.size() == 0) {
                constraint_matrix.resize(3, 3);
                constraint_matrix << 0, 1, 1,
                                     1, 0, 1,
                                     1, 1, 0;
            }
            return constraint_matrix;
        }
        static Eigen::VectorXd get_constraint_rhs(const Eigen::VectorXi& l) {
            return Eigen::Vector3d(l[0] - 2,
                                   l[1] - 1,
                                   l[2] - 1);
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
        //selected vertices are corners
        template <typename PatchT>
        static void generate_subtopology(const PatchParam& param, PatchT& patch) {
            /*
        |        C2
        |       /  \
        |      /    \ 
        |     /      \
        |   C0---V0---C1
            */
            patch.Clear();
            /*typename PatchT::VHandle C[3];
            for (int i = 0; i < 3; ++i) C[i] = add_tagged_vertex(patch, i, true );
            auto V0 = add_tagged_vertex(patch, 0, false);
        
            patch.add_face(C[0], V0, C[1], C[2]);*/
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
            // Populate index sides
            side_indexL[0]=1;
            side_indexR[0]=0;
            side_indexL[1]=0;
            side_indexR[1]=2;
            side_indexL[2]=2;
            side_indexR[2]=1;
            side_indexL[3]=0; //V0
            side_indexR[3]=0;

            typename PatchT::FaceIterator pfi=vcg::tri::Allocator<PatchT>::AddFaces(patch,1);
            (*pfi).Alloc(4);
            (*pfi).V(0)=&patch.vert[0];(*pfi).V(1)=&patch.vert[3];(*pfi).V(2)=&patch.vert[1];(*pfi).V(3)=&patch.vert[2];
//            cout<<"patch 3 -0 "<<endl;
        }
        /*static VariableIndicators& get_variable_indicators() {
            static VariableIndicators variable_indicators;
            return variable_indicators;
        }*/
        static std::string get_param_str(const PatchParam& param) { return ""; }
    };
}
