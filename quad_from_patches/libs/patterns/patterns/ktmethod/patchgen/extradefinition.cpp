#include "patterns/ktmethod/patchgen/Pattern.h"

namespace  kt84{
    Eigen::VectorXd make_Vector6d( double a1, double a2,double a3,double a4,double a5,double a6){
        Eigen::VectorXd result(6);
        result<< a1,a2,a3,a4,a5,a6;
        return result;
    }
    Eigen::VectorXd make_Vector7d( double a0, double a1, double a2,double a3,double a4,double a5,double a6){
        Eigen::VectorXd result(7);
        result<< a0,a1,a2,a3,a4,a5,a6;
        return result;
    }
    Eigen::VectorXd make_Vector5d( double a1, double a2,double a3,double a4,double a5){
        Eigen::VectorXd result(5);
        result<< a1,a2,a3,a4,a5;
        return result;
    }
    Eigen::VectorXd make_Vector8d( double a1, double a2,double a3,double a4,double a5, double a6,double a7,double a8){
        Eigen::VectorXd result(8);
        result<< a1,a2,a3,a4,a5,a6,a7,a8;
        return result;
    }
    Eigen::VectorXd make_Vector9d( double a1, double a2,double a3,double a4,double a5, double a6,double a7,double a8,double a9){
        Eigen::VectorXd result(9);
        result<< a1,a2,a3,a4,a5,a6,a7,a8,a9;
        return result;
    }
    Eigen::VectorXd make_Vector10d( double a1, double a2,double a3,double a4,double a5, double a6, double a7,double a8,double a9,double a10){
        Eigen::VectorXd result(10);
        result<< a1,a2,a3,a4,a5,a6,a7,a8,a9,a10;
        return result;
    }

}
