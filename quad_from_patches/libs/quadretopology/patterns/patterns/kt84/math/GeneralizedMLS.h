#pragma once

#include "../PolynomialBasisGen.h"
#include "RBFKernel.h"
#include "FiniteDifferentiator.h"
#include <Eigen/Dense>
#include <boost/tuple/tuple.hpp>
#include <boost/static_assert.hpp>
#include <vector>

namespace kt84 {

template <int _DimIn, int _DimOut, class _RBFKernel_Core, int _DegreePolynomial>
struct GeneralizedMLS
    : public FiniteDifferentiator<GeneralizedMLS<_DimIn, _DimOut, _RBFKernel_Core, _DegreePolynomial>, _DimIn, _DimOut>
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
    typedef boost::tuple<Point, Value, Gradient> Constraint;
    
    BOOST_STATIC_ASSERT(Kernel::Decaying);
    
    std::vector<Constraint> constraints;
    Kernel kernel;
    
    void clear_constraints() {
        constraints.clear();
    }
    void add_constraint(const Point& point, const Value& value, const Gradient& gradient) {
        constraints.push_back(Constraint(point, value, gradient));
    }
    Value operator()(const Point& point) const {
        const int P = PolynomialBasisGen::DimOut;
        Eigen::Matrix<double, P, P> A;
        Eigen::Matrix<double, P, DimOut> b;
        A.setZero();
        b.setZero();
        for (size_t i = 0; i < constraints.size(); ++i) {
            const Point& point_i = constraints[i].get<0>();
            const Value& value_i = constraints[i].get<1>();
            const Gradient& gradient_i = constraints[i].get<2>();
            double w_i = kernel(point, point_i);
            if (w_i == std::numeric_limits<double>::infinity())
                return value_i;
            w_i = w_i * w_i;
            // value constraint
            PolynomialBasisGen::Basis basis_i = PolynomialBasisGen::basis(point_i);
            A += w_i * basis_i * basis_i.transpose();
            b += w_i * basis_i * value_i.transpose();
            // gradient constraint
            PolynomialBasisGen::Gradient basis_i_gradient = PolynomialBasisGen::gradient(point_i);
            A += w_i * basis_i_gradient * basis_i_gradient.transpose();
            b += w_i * basis_i_gradient * gradient_i.transpose();
        }
        Eigen::Matrix<double, P, DimOut> c = A.colPivHouseholderQr().solve(b);
        return c.transpose() * PolynomialBasisGen::basis(point);
    }
};

}

