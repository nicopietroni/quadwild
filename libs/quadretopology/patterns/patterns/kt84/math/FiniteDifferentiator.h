#pragma once

#include <Eigen/Core>
#include <boost/static_assert.hpp>

namespace kt84 {

template <class _Func, int _DimIn, int _DimOut>
struct FiniteDifferentiator {
    typedef Eigen::Matrix<double, _DimIn , 1> Point;
    typedef Eigen::Matrix<double, _DimOut, _DimIn> Gradient;
    typedef Eigen::Matrix<double, _DimIn, _DimIn> Hessian;
    
    Gradient gradient_fd(const Point& point, double epsilon) const {
        Gradient result;
        const _Func& f = *reinterpret_cast<const _Func*>(this);
        for (int i = 0; i < _DimIn; ++i) {
            Point delta = Point::Unit(_DimIn, i) * epsilon;
            result.col(i) = (f(point + delta) - f(point - delta)) / (2 * epsilon);
        }
        return result;
    }
    Hessian hessian_fd(const Point& point, double epsilon) const {
        BOOST_STATIC_ASSERT(_DimOut == 1);
        Hessian result;
        for (int i = 0; i < _DimIn; ++i) {
            Point delta = Point::Unit(_DimIn, i) * epsilon;
            result.row(i) = (gradient_fd(point + delta, epsilon) - gradient_fd(point - delta, epsilon)) / (2 * epsilon);
        }
        return result;
    }
};

}

