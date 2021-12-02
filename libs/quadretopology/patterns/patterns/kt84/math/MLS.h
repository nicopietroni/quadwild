#pragma once

#include "PolynomialBasisGen.h"
#include "RBFKernel.h"
#include "FiniteDifferentiator.h"
#include <Eigen/Dense>
#include <utility>
#include <vector>

namespace kt84 {

template <int _DimIn, int _DimOut, class _RBFKernel_Core, int _DegreePolynomial>
struct MLS
    : public FiniteDifferentiator<MLS<_DimIn, _DimOut, _RBFKernel_Core, _DegreePolynomial>, _DimIn, _DimOut>
{
    enum {
        DimIn  = _DimIn,
        DimOut = _DimOut,
        DegreePolynomial = _DegreePolynomial,
    };
    
    typedef Eigen::Matrix<double, DimIn , 1> Point;
    typedef Eigen::Matrix<double, DimOut, 1> Value;
    typedef Eigen::Matrix<double, DimOut, DimIn> Gradient;
    typedef PolynomialBasisGen<DimIn, DegreePolynomial> PolynomialBasisGen;
    typedef RBFKernel_Bivariate<DimIn, _RBFKernel_Core> Kernel;
    typedef std::pair<Point, Value> Constraint;
    
    static_assert(Kernel::Decaying, "MLS requires its kernel to be decaying!");
    
    std::vector<Constraint> constraints;
    Kernel kernel;
    
    void clear_constraints() {
        constraints.clear();
    }
    void add_constraint(const Point& point, const Value& value) {
        constraints.push_back(Constraint(point, value));
    }
    Value operator()(const Point& point) const {
        const int P = PolynomialBasisGen::DimOut;
        Eigen::Matrix<double, P, P> A;
        Eigen::Matrix<double, P, DimOut> b;
        A.setZero();
        b.setZero();
        for (size_t i = 0; i < constraints.size(); ++i) {
            const Point& point_i = constraints[i].first;
            const Value& value_i = constraints[i].second;
            double w_i = kernel(point, point_i);
            if (w_i == std::numeric_limits<double>::infinity())
                return value_i;
            w_i = w_i * w_i;
            PolynomialBasisGen::Basis basis_i = PolynomialBasisGen::basis(point_i);
            A += w_i * basis_i * basis_i.transpose();
            b += w_i * basis_i * value_i.transpose();
        }
        Eigen::Matrix<double, P, DimOut> c = A.colPivHouseholderQr().solve(b);
        return c.transpose() * PolynomialBasisGen::basis(point);
    }
};

// template specialization for 0-degree polynomial due to Eigen's inability to handle 1x1 matrix (which is scalar) coherently
template <int _DimIn, int _DimOut, class _RBFKernel_Core>
struct MLS<_DimIn, _DimOut, _RBFKernel_Core, 0>
    : public FiniteDifferentiator<MLS<_DimIn, _DimOut, _RBFKernel_Core, 0>, _DimIn, _DimOut>
{
    enum {
        DimIn  = _DimIn,
        DimOut = _DimOut,
        DegreePolynomial = 0,
    };
    
    typedef Eigen::Matrix<double, DimIn , 1> Point;
    typedef Eigen::Matrix<double, DimOut, 1> Value;
    typedef Eigen::Matrix<double, DimOut, DimIn> Gradient;
    typedef PolynomialBasisGen<DimIn, 0> PolynomialBasisGen;
    typedef RBFKernel_Bivariate<DimIn, _RBFKernel_Core> Kernel;
    typedef std::pair<Point, Value> Constraint;
    
    static_assert(Kernel::Decaying, "MLS requires its kernel to be decaying!");
    
    std::vector<Constraint> constraints;
    Kernel kernel;
    
    void clear_constraints() {
        constraints.clear();
    }
    void add_constraint(const Point& point, const Value& value) {
        constraints.push_back(Constraint(point, value));
    }
    Value operator()(const Point& point) const {
        Eigen::Matrix<double, 1, DimOut> b;
        double A = 0;
        b.setZero();
        for (size_t i = 0; i < constraints.size(); ++i) {
            const Point& point_i = constraints[i].first;
            const Value& value_i = constraints[i].second;
            double w_i = kernel(point, point_i);
            if (w_i == std::numeric_limits<double>::infinity())
                return value_i;
            w_i = w_i * w_i;
            PolynomialBasisGen::Basis basis_i = PolynomialBasisGen::basis(point_i);
            A += w_i * basis_i * basis_i.transpose();
            b += w_i * basis_i * value_i.transpose();
        }
        Eigen::Matrix<double, 1, DimOut> c = b / A;
        return c.transpose() * PolynomialBasisGen::basis(point);
    }
};


}

