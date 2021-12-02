#pragma once

#include "PolynomialBasisGen.h"
#include "RBFKernel.h"
#include "FiniteDifferentiator.h"
#include <Eigen/Dense>
#include <utility>
#include <vector>

namespace kt84 {

template <int _DimIn, int _DimOut, class _RBFKernel_Core, int _DegreePolynomial>
struct RBF
    : public FiniteDifferentiator<RBF<_DimIn, _DimOut, _RBFKernel_Core, _DegreePolynomial>, _DimIn, _DimOut>
{
    enum {
        DimIn  = _DimIn,
        DimOut = _DimOut,
        DegreePolynomial = _DegreePolynomial,
    };
    
    typedef Eigen::Matrix<double, DimIn , 1> Point;                 // TODO: treat 1x1 matrix as scalar using Matrix11ToScalar
    typedef Eigen::Matrix<double, DimOut, 1> Value;
    typedef Eigen::Matrix<double, DimOut, DimIn> Gradient;
    typedef PolynomialBasisGenT<DimIn, DegreePolynomial> PolynomialBasisGen;
    typedef RBFKernel_Bivariate<DimIn, _RBFKernel_Core> Kernel;
    typedef std::pair<Point, Value> Constraint;
    
    std::vector<Constraint> constraints;
    Kernel kernel;
    Eigen::Matrix<double, -1, DimOut> weights;
    Eigen::MatrixXd A_matrix;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> A_factorized;
    
    void clear_constraints() {
        constraints.clear();
    }
    void add_constraint(const Point& point, const Value& value) {
        constraints.push_back(Constraint(point, value));
    }
    void factorize() {
        const int P = PolynomialBasisGen::DimOut;
        const size_t n = constraints.size();
        const int m = n + P;
        A_matrix = Eigen::MatrixXd::Zero(m, m);
        for (size_t i = 0; i < n; ++i) {
            const Point& point_i = constraints[i].first;
            // rbf part
            A_matrix(i, i) = kernel.univariate(0);
            for (size_t j = i + 1; j < n; ++j) {
                const Point& point_j = constraints[j].first;
                A_matrix(i, j) = A_matrix(j, i) = kernel(point_i, point_j);
            }
            // polynomial part
            A_matrix.block<P, 1>(n, i) << PolynomialBasisGen::basis(point_i);
            A_matrix.block<1, P>(i, n) = A_matrix.block<P, 1>(n, i).transpose();
        }
        A_factorized.compute(A_matrix);         // factorize
    }
    void solve() {
        const int P = PolynomialBasisGen::DimOut;
        const size_t n = constraints.size();
        const int m = n + P;
        Eigen::Matrix<double, -1, DimOut> b;
        b.setZero(m, DimOut);
        // constraint part
        for (size_t i = 0; i < n; ++i) {
            const Value& value_i = constraints[i].second;
            b.row(i).transpose() << value_i;
        }
        // polynomial part is just 0
        weights = A_factorized.solve(b);        // solve
    }
    void factorize_and_solve() {
        factorize();
        solve();
    }
    Value operator()(const Point& point) const {
        const int P = PolynomialBasisGen::DimOut;
        Value result = Value::Zero();
        // rbf part
        for (size_t i = 0; i < constraints.size(); ++i) {
            const Point& point_i = constraints[i].first;
            result += kernel(point, point_i) * weights.row(i).transpose();
        }
        // polynomial part
        auto basis = PolynomialBasisGen::basis(point);
        result += (basis.transpose() * weights.bottomRows(P)).transpose();
        return result;
    }
    Gradient gradient(const Point& point) const {
        const int P = PolynomialBasisGen::DimOut;
        Gradient result = Gradient::Zero();
        // rbf part
        for (size_t i = 0; i < constraints.size(); ++i) {
            const Point& point_i = constraints[i].first;
            result += weights.row(i).transpose() * kernel.gradient(point, point_i);
        }
        // polynomial part
        auto b_gradient = PolynomialBasisGen::gradient(point);
        result += weights.bottomRows(P).transpose() * b_gradient;
        return result;
    }
};

}

