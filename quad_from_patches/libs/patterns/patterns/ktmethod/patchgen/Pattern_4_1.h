#pragma once
#include "Pattern.h"
#include "ILP.h"
#include <sstream>
#include <vcg/complex/algorithms/create/platonic.h>
#include "edgeloop.h"
/*
equation for pattern 1:
  |0|     |1|     |0|     |1|   |2|    |1|   |l0|
p0|1| + p1|0| + p2|1| + p3|0| + |2| + x|1| = |l1|
  |0|     |1|     |0|     |1|   |1|    |0|   |l2|
  |1|     |0|     |1|     |0|   |1|    |0|   |l3|
*/
namespace patchgen {
    template <>
    struct Pattern<4, 1> {
        static Eigen::MatrixXd& get_constraint_matrix() {
            static Eigen::MatrixXd constraint_matrix;
            if (constraint_matrix.size() == 0) {
                constraint_matrix.resize(4, 5);
                constraint_matrix << 0, 1, 0, 1, 1,
                                     1, 0, 1, 0, 1,
                                     0, 1, 0, 1, 0,
                                     1, 0, 1, 0, 0;
            }
            return constraint_matrix;
        }
    
        static Eigen::VectorXd get_constraint_rhs(const Eigen::VectorXi& l) {
            return Eigen::Vector4d(l[0] - 2,
                                   l[1] - 2,
                                   l[2] - 1,
                                   l[3] - 1);
        }
    
        static int& get_variable(PatchParam& param, int index) {
            if (index < 4) return param.p[index];
            return param.x;
        }

        static bool get_default_parameter(const Eigen::VectorXi& l, PatchParam& param) {
            auto& constraint_matrix = get_constraint_matrix();
            const int num_variables = constraint_matrix.cols();
            ILP ilp(num_variables);
            ilp.add_constraint(constraint_matrix, EQ, get_constraint_rhs(l));
        
            // arbitrary constraints and objective
            // p0 <= p2
            ilp.add_constraint(kt84::make_Vector5d(1, 0, -1, 0, 0), LE, 0);
            // p1 <= p3
            ilp.add_constraint(kt84::make_Vector5d(0, 1, 0, -1, 0), LE, 0);
            // maximize p0+p1
            ilp.set_objective(kt84::make_Vector5d(1, 1, 0, 0, 0), true);
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
        |   C3--------C2
        |   |\        |
        |   | \       | <--x
        |   |  \      |
        |   |   \     |
        |   |    V2---V1
        |   |    |    |
        |   |    |    |
        |   |    |    |
        |   |    |    |
        |   C0---V0---C1
        |      ^--x
            */
        /*    patch.clear();
            typename PatchT::VHandle C[4];
            typename PatchT::VHandle V[3];
            for (int i = 0; i < 4; ++i) C[i] = add_tagged_vertex(patch, i, true );
            for (int i = 0; i < 3; ++i) V[i] = add_tagged_vertex(patch, i, false);
            patch.add_face(C[0], V[0], V[2], C[3]);
            patch.add_face(V[0], C[1], V[1], V[2]);
            patch.add_face(V[1], C[2], C[3], V[2]);
        
            auto h_insert_x = patch.halfedge_handle(V[0]);  // corresponds to V0-C0
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
            vcg::tri::Allocator<PatchT>::AddVertices(patch,7);
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
            side_indexL[4]=0;   //V0
            side_indexR[4]=0;
            side_indexL[5]=3;   //V1
            side_indexR[5]=3;
            side_indexL[6]=-1;   //V2
            side_indexR[6]=-1;

            typename PatchT::FaceIterator pfi=vcg::tri::Allocator<PatchT>::AddFaces(patch,3);
            (*pfi).Alloc(4);
            (*pfi).V(0)=&patch.vert[0];(*pfi).V(1)=&patch.vert[4];(*pfi).V(2)=&patch.vert[6];(*pfi).V(3)=&patch.vert[3];
            pfi++;
            (*pfi).Alloc(4);
            (*pfi).V(0)=&patch.vert[4];(*pfi).V(1)=&patch.vert[1];(*pfi).V(2)=&patch.vert[5];(*pfi).V(3)=&patch.vert[6];
            pfi++;
            (*pfi).Alloc(4);
            (*pfi).V(0)=&patch.vert[5];(*pfi).V(1)=&patch.vert[2];(*pfi).V(2)=&patch.vert[3];(*pfi).V(3)=&patch.vert[6];

            vcg::face::Pos<typename PatchT::FaceType> startPos=getPosFromIndex(patch,0,4);
            for (int i = 0; i < param.x; ++i)
               insert_edgeloop<PatchT>(patch,startPos);

//            cout<<"patch 4 -- 1 "<<endl;

        }
        /*static VariableIndicators& get_variable_indicators() {
            static VariableIndicators variable_indicators;
            if (variable_indicators.empty()) {
                variable_indicators.resize(1);
                //variable_indicators[0].push_back(std::make_pair(PatchVertexTag::C0, PatchVertexTag::V0));       // for x (not editable anyway)
                //variable_indicators[0].push_back(std::make_pair(PatchVertexTag::C2, PatchVertexTag::V1));
                //variable_indicators[0].push_back(std::make_pair(PatchVertexTag::C3, PatchVertexTag::V2));
            }
            return variable_indicators;
        }*/
        static std::string get_param_str(const PatchParam& param) {
            std::stringstream ss;
            ss << "p0=" << param.p[0] 
               << "_p1=" << param.p[1];
            return ss.str();
        }
    };
}
