#pragma once
#include "Pattern.h"
#include "ILP.h"
#include <sstream>
#include <vcg/complex/algorithms/create/platonic.h>
#include "edgeloop.h"
/*
equation for pattern 3:
  |0|     |1|     |0|     |0|     |0|     |1|     |0|   |4|    |2|    |1|    |1|   |l0|
  |1|     |0|     |1|     |0|     |0|     |0|     |0|   |2|    |0|    |1|    |0|   |l1|
p0|0| + p1|1| + p2|0| + p3|1| + p4|0| + p5|0| + q3|1| + |1| + x|0| + y|0| + z|0| = |l2|
  |0|     |0|     |1|     |0|     |1|     |0|     |0|   |1|    |0|    |0|    |1|   |l3|
  |0|     |0|     |0|     |1|     |0|     |1|     |1|   |1|    |0|    |0|    |0|   |l4|
  |1|     |0|     |0|     |0|     |1|     |0|     |0|   |1|    |0|    |0|    |0|   |l5|
*/
namespace patchgen {
    template <>
    struct Pattern<6, 3> {
        static Eigen::MatrixXd& get_constraint_matrix() {
            static Eigen::MatrixXd constraint_matrix;
            if (constraint_matrix.size() == 0) {
                constraint_matrix.resize(6, 10);
                constraint_matrix << 0, 1, 0, 0, 0, 1, 0, 2, 1, 1,
                                     1, 0, 1, 0, 0, 0, 0, 0, 1, 0,
                                     0, 1, 0, 1, 0, 0, 1, 0, 0, 0,
                                     0, 0, 1, 0, 1, 0, 0, 0, 0, 1,
                                     0, 0, 0, 1, 0, 1, 1, 0, 0, 0,
                                     1, 0, 0, 0, 1, 0, 0, 0, 0, 0;
            }
            return constraint_matrix;
        }
    
        static Eigen::VectorXd get_constraint_rhs(const Eigen::VectorXi& l) {
            return kt84::make_Vector6d(l[0] - 4,
                                       l[1] - 2,
                                       l[2] - 1,
                                       l[3] - 1,
                                       l[4] - 1,
                                       l[5] - 1);
        }
    
        static int& get_variable(PatchParam& param, int index) {
            if (index < 6) return param.p[index];
            if (index == 6) return param.q[3];
            if (index == 7) return param.x;
            if (index == 8) return param.y;
            return param.z;
        }

        static bool get_default_parameter(const Eigen::VectorXi& l, PatchParam& param) {
            auto& constraint_matrix = get_constraint_matrix();
            const int num_variables = constraint_matrix.cols();
            ILP ilp(num_variables);
            ilp.add_constraint(constraint_matrix, EQ, get_constraint_rhs(l));
        
            // arbitrary constraints
            // xmin
            ilp.set_objective(kt84::make_Vector10d(0, 0, 0, 0, 0, 0, 0, 1, 0, 0), false);
            if (!ilp.solve()) return false;
            int xmin = ilp.get_variables()[7];
            // xmax
            ilp.refresh();
            ilp.set_objective(kt84::make_Vector10d(0, 0, 0, 0, 0, 0, 0, 1, 0, 0), true);
            if (!ilp.solve()) return false;
            int xmax = ilp.get_variables()[7];
            // xmid-1<=x<=xmid+1
            ilp.refresh();
            int xmid = (xmin + xmax) / 2;
            ilp.add_constraint(kt84::make_Vector10d(0, 0, 0, 0, 0, 0, 0, 1, 0, 0), LE, xmid + 1);
            ilp.add_constraint(kt84::make_Vector10d(0, 0, 0, 0, 0, 0, 0, 1, 0, 0), GE, xmid - 1);
            // maximize p0+p1+p2+p3+p4+p5
            ilp.set_objective(kt84::make_Vector10d(1, 1, 1, 1, 1, 1, 0, 0, 0, 0), true);
            if (!ilp.solve()) return false;
        
            auto variables = ilp.get_variables();
            for (int i = 0; i < num_variables; ++i)
                get_variable(param, i) = variables[i];
        
            param.pattern_id = 3;
            return true;
        }
        template <typename PatchT>
        static void generate_subtopology(const PatchParam& param, PatchT& patch) {
            /*
        |                    v--z                  
        |            C4------------C3              
        |           /  \          /  \             
        |          /    \        /    \            
        |  q3-->  /      \      /      \    <--q3
        |        /        V7--V8        \          
        |       /        /     |         \         
        |      /        /      |          \        
        |     /        /       |           \       
        |   C5--------V6-------V5----------C2      
        |    \         |      /  \        /        
        |     \        |     /    \      /  <--y
        |      \       |     |     V4---V3         
        |       \      |     /    /    /           
        |        \     |    /    /    /            
        |        C0---V0---V1---V2---C1            
        |        x--^    A   ^--y  ^--x
        |                |                         
        |                z                         
            */
         /*   patch.clear();
            typename PatchT::VHandle C[6];
            typename PatchT::VHandle V[9];
            for (int i = 0; i < 6; ++i) C[i] = add_tagged_vertex(patch, i, true );
            for (int i = 0; i < 9; ++i) V[i] = add_tagged_vertex(patch, i, false);
        
            patch.add_face(C[0], V[0], V[6], C[5]);
            patch.add_face(V[0], V[1], V[5], V[6]);
            patch.add_face(V[1], V[2], V[4], V[5]);
            patch.add_face(V[2], C[1], V[3], V[4]);
            patch.add_face(V[3], C[2], V[5], V[4]);
            patch.add_face(C[2], C[3], V[8], V[5]);
            patch.add_face(C[3], C[4], V[7], V[8]);
            patch.add_face(C[4], C[5], V[6], V[7]);
            patch.add_face(V[5], V[8], V[7], V[6]);
        
            auto h_insert_x = patch.halfedge_handle(V[0]);  // corresponds to V0-C0
            for (int i = 0; i < param.x; ++i)
                insert_edgeloop(patch, h_insert_x);
            auto h_insert_y = patch.halfedge_handle(V[2]);  // corresponds to V2-V1
            for (int i = 0; i < param.y; ++i)
                insert_edgeloop(patch, h_insert_y);
            auto h_insert_z = patch.halfedge_handle(V[1]);  // corresponds to V1-V0
            for (int i = 0; i < param.z; ++i)
                insert_edgeloop(patch, h_insert_z);
            auto h_insert_q3 = patch.halfedge_handle(C[3]);  // corresponds to C3-C2
            for (int i = 0; i < param.q[3]; ++i)
                insert_edgeloop(patch, h_insert_q3);
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
            vcg::tri::Allocator<PatchT>::AddVertices(patch,15);
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
            side_indexL[6]=0;  //V0
            side_indexR[6]=0;
            side_indexL[7]=0;  //V1
            side_indexR[7]=0;
            side_indexL[8]=0;  //V2
            side_indexR[8]=0;
            side_indexL[9]=5;  //V3
            side_indexR[9]=5;
            side_indexL[10]=-1;  //V4
            side_indexR[10]=-1;
            side_indexL[11]=-1;  //V5
            side_indexR[11]=-1;
            side_indexL[12]=-1;  //V6
            side_indexR[12]=-1;
            side_indexL[13]=-1;  //V7
            side_indexR[13]=-1;
            side_indexL[14]=-1;  //V8
            side_indexR[14]=-1;

            typename PatchT::FaceIterator pfi=vcg::tri::Allocator<PatchT>::AddFaces(patch,9);
            (*pfi).Alloc(4); //C[0], V[0], V[6], C[5]
            (*pfi).V(0)=&patch.vert[0];(*pfi).V(1)=&patch.vert[6];(*pfi).V(2)=&patch.vert[12];(*pfi).V(3)=&patch.vert[5];
            pfi++;
            (*pfi).Alloc(4); //V[0], V[1], V[5], V[6]
            (*pfi).V(0)=&patch.vert[6];(*pfi).V(1)=&patch.vert[7];(*pfi).V(2)=&patch.vert[11];(*pfi).V(3)=&patch.vert[12];
            pfi++;
            (*pfi).Alloc(4); //V[1], V[2], V[4], V[5]
            (*pfi).V(0)=&patch.vert[7];(*pfi).V(1)=&patch.vert[8];(*pfi).V(2)=&patch.vert[10];(*pfi).V(3)=&patch.vert[11];
            pfi++;
            (*pfi).Alloc(4); //V[2], C[1], V[3], V[4]
            (*pfi).V(0)=&patch.vert[8];(*pfi).V(1)=&patch.vert[1];(*pfi).V(2)=&patch.vert[9];(*pfi).V(3)=&patch.vert[10];
            pfi++;
            (*pfi).Alloc(4); //V[3], C[2], V[5], V[4]
            (*pfi).V(0)=&patch.vert[9];(*pfi).V(1)=&patch.vert[2];(*pfi).V(2)=&patch.vert[11];(*pfi).V(3)=&patch.vert[10];
            pfi++;
            (*pfi).Alloc(4); //C[2], C[3], V[8], V[5]
            (*pfi).V(0)=&patch.vert[2];(*pfi).V(1)=&patch.vert[3];(*pfi).V(2)=&patch.vert[14];(*pfi).V(3)=&patch.vert[11];
            pfi++;
            (*pfi).Alloc(4); //C[3], C[4], V[7], V[8]
            (*pfi).V(0)=&patch.vert[3];(*pfi).V(1)=&patch.vert[4];(*pfi).V(2)=&patch.vert[13];(*pfi).V(3)=&patch.vert[14];
            pfi++;
            (*pfi).Alloc(4); //C[4], C[5], V[6], V[7]
            (*pfi).V(0)=&patch.vert[4];(*pfi).V(1)=&patch.vert[5];(*pfi).V(2)=&patch.vert[12];(*pfi).V(3)=&patch.vert[13];
            pfi++;
            (*pfi).Alloc(4); //V[5], V[8], V[7], V[6]
            (*pfi).V(0)=&patch.vert[11];(*pfi).V(1)=&patch.vert[14];(*pfi).V(2)=&patch.vert[13];(*pfi).V(3)=&patch.vert[12];

            vcg::face::Pos<typename PatchT::FaceType> startPos=getPosFromIndex(patch,0,6);
            for (int i = 0; i < param.x; ++i)
               insert_edgeloop<PatchT>(patch,startPos);

            startPos=getPosFromIndex(patch,9,2);
            for (int i = 0; i < param.y; ++i)
               insert_edgeloop<PatchT>(patch,startPos);

            startPos=getPosFromIndex(patch,3,4);
            for (int i = 0; i < param.z; ++i)
               insert_edgeloop<PatchT>(patch,startPos);

            startPos=getPosFromIndex(patch,4,5);
            for (int i = 0; i < param.q[3]; ++i)
               insert_edgeloop<PatchT>(patch,startPos);

//            cout<<"patch 6 -- 3 "<<endl;
        }
        /*static VariableIndicators& get_variable_indicators() {
            static VariableIndicators variable_indicators;
            if (variable_indicators.empty()) {
                variable_indicators.resize(4);
                variable_indicators[0].push_back(std::make_pair(PatchVertexTag::C4, PatchVertexTag::C5));       // for q3
                variable_indicators[0].push_back(std::make_pair(PatchVertexTag::V7, PatchVertexTag::V6));
                variable_indicators[0].push_back(std::make_pair(PatchVertexTag::V8, PatchVertexTag::V5));
                variable_indicators[0].push_back(std::make_pair(PatchVertexTag::C3, PatchVertexTag::C2));
                variable_indicators[1].push_back(std::make_pair(PatchVertexTag::C0, PatchVertexTag::V0));       // for x
                variable_indicators[1].push_back(std::make_pair(PatchVertexTag::C5, PatchVertexTag::V6));
                variable_indicators[1].push_back(std::make_pair(PatchVertexTag::C4, PatchVertexTag::V7));
                variable_indicators[1].push_back(std::make_pair(PatchVertexTag::C3, PatchVertexTag::V8));
                variable_indicators[1].push_back(std::make_pair(PatchVertexTag::C2, PatchVertexTag::V5));
                variable_indicators[1].push_back(std::make_pair(PatchVertexTag::V3, PatchVertexTag::V4));
                variable_indicators[1].push_back(std::make_pair(PatchVertexTag::C1, PatchVertexTag::V2));
                variable_indicators[2].push_back(std::make_pair(PatchVertexTag::V1, PatchVertexTag::V2));       // for y
                variable_indicators[2].push_back(std::make_pair(PatchVertexTag::V5, PatchVertexTag::V4));
                variable_indicators[2].push_back(std::make_pair(PatchVertexTag::C2, PatchVertexTag::V3));
                variable_indicators[3].push_back(std::make_pair(PatchVertexTag::V0, PatchVertexTag::V1));       // for z
                variable_indicators[3].push_back(std::make_pair(PatchVertexTag::V6, PatchVertexTag::V5));
                variable_indicators[3].push_back(std::make_pair(PatchVertexTag::V7, PatchVertexTag::V8));
                variable_indicators[3].push_back(std::make_pair(PatchVertexTag::C4, PatchVertexTag::C3));
            }
            return variable_indicators;
        }*/
        static std::string get_param_str(const PatchParam& param) {
            std::stringstream ss;
            ss << "q3=" << param.q[3]
               << "_x=" << param.x
               << "_y=" << param.y
               << "_z=" << param.z;
            return ss.str();
        }
    };
}
