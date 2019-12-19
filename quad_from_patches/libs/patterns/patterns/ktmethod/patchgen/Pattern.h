#pragma once
#include "PatchParam.h"
//#include "PatchVertexTag.h"
#include <vector>
#include <string>
#include <utility>
#include <functional>
namespace patchgen {
    //typedef std::vector<std::vector<std::pair<PatchVertexTag, PatchVertexTag>>> VariableIndicators;
    
    template <int NumSides, int PatternID>
    struct Pattern {
        static Eigen::MatrixXd& get_constraint_matrix();
        
        static Eigen::VectorXd get_constraint_rhs(const Eigen::VectorXi& l);
        
        static int& get_variable(PatchParam& param, int index);
        
        static bool get_default_parameter(const Eigen::VectorXi& l, PatchParam& param);
        
        template <typename PatchT>
        static void generate_subtopology(const PatchParam& param, PatchT& patch);
        
        //static VariableIndicators& get_variable_indicators();
        
        static std::string get_param_str(const PatchParam& param);
    };
}
namespace  kt84{
    Eigen::VectorXd make_Vector6d( double a1, double a2,double a3,double a4,double a5,double a6);
    Eigen::VectorXd make_Vector7d( double a0, double a1, double a2,double a3,double a4,double a5,double a6);
    Eigen::VectorXd make_Vector5d( double a1, double a2,double a3,double a4,double a5);
    Eigen::VectorXd make_Vector8d( double a1, double a2,double a3,double a4,double a5, double a6,double a7,double a8);
    Eigen::VectorXd make_Vector9d( double a1, double a2,double a3,double a4,double a5, double a6,double a7,double a8,double a9);
    Eigen::VectorXd make_Vector10d( double a1, double a2,double a3,double a4,double a5, double a6, double a7,double a8,double a9,double a10);

    //template <class MeshType,class HalfEdgeType>
    //void insert_edgeloop(MeshType& mesh,HalfEdgeType& h_start);
}
