#pragma once
#include "Pattern.h"
#include "ILP.h"
#include <vcg/complex/algorithms/create/platonic.h>
#include "edgeloop.h"
/*
equation for pattern 0:
  |0|     |1|     |0|     |0|     |1|   |2|   |l0|
  |1|     |0|     |1|     |0|     |0|   |1|   |l1|
p0|0| + p1|1| + p2|0| + p3|1| + p4|0| + |1| = |l2|
  |0|     |0|     |1|     |0|     |1|   |1|   |l3|
  |1|     |0|     |0|     |1|     |0|   |1|   |l4|
*/
namespace patchgen {
    template <>
    struct Pattern<5, 0> {
        static Eigen::MatrixXd& get_constraint_matrix() {
            static Eigen::MatrixXd constraint_matrix;
            if (constraint_matrix.size() == 0) {
                constraint_matrix.resize(5, 5);
                constraint_matrix << 0, 1, 0, 0, 1,
                                     1, 0, 1, 0, 0,
                                     0, 1, 0, 1, 0,
                                     0, 0, 1, 0, 1,
                                     1, 0, 0, 1, 0;
            }
            return constraint_matrix;
        }
    
        static Eigen::VectorXd get_constraint_rhs(const Eigen::VectorXi& l) {
            return kt84::make_Vector5d(l[0] - 2,
                                       l[1] - 1,
                                       l[2] - 1,
                                       l[3] - 1,
                                       l[4] - 1);
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
        |        ___C3__
        |      _/   |   \_   
        |     /     |     \
        |   C4      |      C2
        |     \     |     /
        |      \    |    /
        |      C0---V0--C1
            */
        /*    patch.clear();
            typename PatchT::VHandle C[5];
            for (int i = 0; i < 5; ++i) C[i] = add_tagged_vertex(patch, i, true );
            auto V0 = add_tagged_vertex(patch, 0, false);
        
            patch.add_face(V0, C[3], C[4], C[0]);
            patch.add_face(V0, C[1], C[2], C[3]);*/
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

            // Populate index sides
            side_indexL[0]=1;
            side_indexR[0]=0;
            side_indexL[1]=0;
            side_indexR[1]=4;
            side_indexL[2]=4;
            side_indexR[2]=3;
            side_indexL[3]=3;
            side_indexR[3]=2;
            side_indexL[4]=2;
            side_indexR[4]=1;
            side_indexL[5]=0;   //V0
            side_indexR[5]=0;

            typename PatchT::FaceIterator pfi=vcg::tri::Allocator<PatchT>::AddFaces(patch,2);
            (*pfi).Alloc(4); //V0, C[3], C[4], C[0]
            (*pfi).V(0)=&patch.vert[5];(*pfi).V(1)=&patch.vert[3];(*pfi).V(2)=&patch.vert[4];(*pfi).V(3)=&patch.vert[0];
            pfi++;
            (*pfi).Alloc(4); //V0, C[1], C[2], C[3]
            (*pfi).V(0)=&patch.vert[5];(*pfi).V(1)=&patch.vert[1];(*pfi).V(2)=&patch.vert[2];(*pfi).V(3)=&patch.vert[3];

//            cout<<"patch 5 -0 "<<endl;
        }
        /*static VariableIndicators& get_variable_indicators() {
            static VariableIndicators variable_indicators;
            return variable_indicators;
        }*/
        static std::string get_param_str(const PatchParam& param) { return ""; }
    };
}
