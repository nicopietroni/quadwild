#pragma once
#include "Pattern.h"
#include "ILP.h"
#include <sstream>
#include <vcg/complex/algorithms/create/platonic.h>
#include "edgeloop.h"
/*
equation for pattern 4:
  |0|     |1|     |0|     |1|     |1|   |4|    |2|    |1|   |l0|
p0|1| + p1|0| + p2|1| + p3|0| + q1|0| + |2| + x|0| + y|1| = |l1|
  |0|     |1|     |0|     |1|     |1|   |1|    |0|    |0|   |l2|
  |1|     |0|     |1|     |0|     |0|   |1|    |0|    |0|   |l3|
*/
namespace patchgen {
    template <>
    struct Pattern<4, 4> {
        static Eigen::MatrixXd& get_constraint_matrix() {
            static Eigen::MatrixXd constraint_matrix;
            if (constraint_matrix.size() == 0) {
                constraint_matrix.resize(4, 7);
                constraint_matrix << 0, 1, 0, 1, 1, 2, 1,
                                     1, 0, 1, 0, 0, 0, 1,
                                     0, 1, 0, 1, 1, 0, 0,
                                     1, 0, 1, 0, 0, 0, 0;
            }
            return constraint_matrix;
        }
    
        static Eigen::VectorXd get_constraint_rhs(const Eigen::VectorXi& l) {
            return Eigen::Vector4d(l[0] - 4,
                                   l[1] - 2,
                                   l[2] - 1,
                                   l[3] - 1);
        }
    
        static int& get_variable(PatchParam& param, int index) {
            if (index < 4) return param.p[index];
            if (index == 4) return param.q[1];
            if (index == 5) return param.x;
            return param.y;
        }

        static bool get_default_parameter(const Eigen::VectorXi& l, PatchParam& param) {
            auto& constraint_matrix = get_constraint_matrix();
            const int num_variables = constraint_matrix.cols();
            ILP ilp(num_variables);
            ilp.add_constraint(constraint_matrix, EQ, get_constraint_rhs(l));
        
            // arbitrary constraints and objective
            // p0 <= p2
            ilp.add_constraint(kt84::make_Vector7d(1, 0, -1, 0, 0, 0, 0), LE, 0);
            // p1 <= p3
            ilp.add_constraint(kt84::make_Vector7d(0, 1, 0, -1, 0, 0, 0), LE, 0);
            // p3 <= q1
            ilp.add_constraint(kt84::make_Vector7d(0, 0, 0, 1, -1, 0, 0), LE, 0);
            // maximize p0+p1
            ilp.set_objective(kt84::make_Vector7d(1, 1, 0, 0, 0, 0, 0), true);
            if (!ilp.solve()) return false;
        
            auto variables = ilp.get_variables();
            for (int i = 0; i < num_variables; ++i)
                get_variable(param, i) = variables[i];
        
            param.pattern_id = 4;
            return true;
        }
        template <typename PatchT>
        static void generate_subtopology(const PatchParam& param, PatchT& patch) {
            /*
        |            v--q1
        |   C3-----------------C2
        |   | \              _/ |
        |   |  \           _/   |
        |   |   \        _/     | <--y
        |   |    \      /       |
        |   |     V6--V5        V3
        |   |     |    |\_    _/|
        |   |     |    |  \  /  |
        |   |     |    |   V4   |
        |   |     |    |   |    |
        |   C0---V0---V1---V2---C1
        |  x--^  q1--^ y--^    ^--x
            */
          /*  patch.clear();
            typename PatchT::VHandle C[4];
            typename PatchT::VHandle V[7];
            for (int i = 0; i < 4; ++i) C[i] = add_tagged_vertex(patch, i, true );
            for (int i = 0; i < 7; ++i) V[i] = add_tagged_vertex(patch, i, false);
            patch.add_face(C[0], V[0], V[6], C[3]);
            patch.add_face(V[0], V[1], V[5], V[6]);
            patch.add_face(V[1], V[2], V[4], V[5]);
            patch.add_face(V[2], C[1], V[3], V[4]);
            patch.add_face(V[3], C[2], V[5], V[4]);
            patch.add_face(C[2], C[3], V[6], V[5]);
        
            auto h_insert_x = patch.halfedge_handle(V[0]);  // corresponds to V0-C0
            for (int i = 0; i < param.x; ++i)
                insert_edgeloop(patch, h_insert_x);
            auto h_insert_y = patch.halfedge_handle(V[2]);  // corresponds to V2-V1
            for (int i = 0; i < param.y; ++i)
                insert_edgeloop(patch, h_insert_y);
            auto h_insert_q1 = patch.halfedge_handle(V[1]);  // corresponds to V1-V0
            for (int i = 0; i < param.q[1]; ++i)
                insert_edgeloop(patch, h_insert_q1);
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
            vcg::tri::Allocator<PatchT>::AddVertices(patch,11);
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
            side_indexL[5]=0;   //V1
            side_indexR[5]=0;
            side_indexL[6]=0;   //V2
            side_indexR[6]=0;
            side_indexL[7]=3;   //V3
            side_indexR[7]=3;
            side_indexL[8]=-1;   //V4
            side_indexR[8]=-1;
            side_indexL[9]=-1;   //V5
            side_indexR[9]=-1;
            side_indexL[10]=-1;   //V6
            side_indexR[10]=-1;

            typename PatchT::FaceIterator pfi=vcg::tri::Allocator<PatchT>::AddFaces(patch,6);
            (*pfi).Alloc(4); //C[0], V[0], V[6], C[3]
            (*pfi).V(0)=&patch.vert[0];(*pfi).V(1)=&patch.vert[4];(*pfi).V(2)=&patch.vert[10];(*pfi).V(3)=&patch.vert[3];
            pfi++;
            (*pfi).Alloc(4);  //V[0], V[1], V[5], V[6]
            (*pfi).V(0)=&patch.vert[4];(*pfi).V(1)=&patch.vert[5];(*pfi).V(2)=&patch.vert[9];(*pfi).V(3)=&patch.vert[10];
            pfi++;
            (*pfi).Alloc(4); //V[1], V[2], V[4], V[5]
            (*pfi).V(0)=&patch.vert[5];(*pfi).V(1)=&patch.vert[6];(*pfi).V(2)=&patch.vert[8];(*pfi).V(3)=&patch.vert[9];
            pfi++;
            (*pfi).Alloc(4); //V[2], C[1], V[3], V[4]
            (*pfi).V(0)=&patch.vert[6];(*pfi).V(1)=&patch.vert[1];(*pfi).V(2)=&patch.vert[7];(*pfi).V(3)=&patch.vert[8];
            pfi++;
            (*pfi).Alloc(4); //V[3], C[2], V[5], V[4]
            (*pfi).V(0)=&patch.vert[7];(*pfi).V(1)=&patch.vert[2];(*pfi).V(2)=&patch.vert[9];(*pfi).V(3)=&patch.vert[8];
            pfi++;
            (*pfi).Alloc(4); //C[2], C[3], V[6], V[5]
            (*pfi).V(0)=&patch.vert[2];(*pfi).V(1)=&patch.vert[3];(*pfi).V(2)=&patch.vert[10];(*pfi).V(3)=&patch.vert[9];

            vcg::face::Pos<typename PatchT::FaceType> startPos=getPosFromIndex(patch,0,4);
            for (int i = 0; i < param.x; ++i)
               insert_edgeloop<PatchT>(patch,startPos);

            startPos=getPosFromIndex(patch,5,6);
            for (int i = 0; i < param.y; ++i)
               insert_edgeloop<PatchT>(patch,startPos);

            startPos=getPosFromIndex(patch,2,3);
            for (int i = 0; i < param.q[1]; ++i)
               insert_edgeloop<PatchT>(patch,startPos);

//            cout<<"patch 4 -- 4"<<endl;
        }
        /*static VariableIndicators& get_variable_indicators() {
            static VariableIndicators variable_indicators;
            if (variable_indicators.empty()) {
                variable_indicators.resize(3);
                variable_indicators[0].push_back(std::make_pair(PatchVertexTag::V0, PatchVertexTag::V1));       // for q1
                variable_indicators[0].push_back(std::make_pair(PatchVertexTag::V5, PatchVertexTag::V6));
                variable_indicators[0].push_back(std::make_pair(PatchVertexTag::C2, PatchVertexTag::C3));
                //variable_indicators[1].push_back(std::make_pair(PatchVertexTag::C0, PatchVertexTag::V0));       // for x (not editable anyway)
                //variable_indicators[1].push_back(std::make_pair(PatchVertexTag::C1, PatchVertexTag::V2));
                //variable_indicators[1].push_back(std::make_pair(PatchVertexTag::C2, PatchVertexTag::V5));
                //variable_indicators[1].push_back(std::make_pair(PatchVertexTag::C3, PatchVertexTag::V6));
                //variable_indicators[1].push_back(std::make_pair(PatchVertexTag::V3, PatchVertexTag::V4));
                //variable_indicators[2].push_back(std::make_pair(PatchVertexTag::C2, PatchVertexTag::V3));       // for y (not editable anyway)
                //variable_indicators[2].push_back(std::make_pair(PatchVertexTag::V1, PatchVertexTag::V2));
                //variable_indicators[2].push_back(std::make_pair(PatchVertexTag::V4, PatchVertexTag::V5));
            }
            return variable_indicators;
        }*/
        static std::string get_param_str(const PatchParam& param) {
            std::stringstream ss;
            ss << "p0=" << param.p[0] 
               << "_p1=" << param.p[1]
               << "_q1=" << param.q[1];
            return ss.str();
        }
    };
}
