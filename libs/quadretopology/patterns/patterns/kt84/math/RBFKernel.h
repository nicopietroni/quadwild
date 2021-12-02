#pragma once

#include <cmath>
#include <limits>
#include <Eigen/Core>

namespace kt84 {

// the definition of functions alpha(r) and beta(r) implemented in each kernel |
//-----------------------------------------------------------------------------+
/*

inner product of vectors x and y:
    x * y := x^T y (scalar)
outer product of vectors x and y: 
    x % y := x y^T (matrix)

notations:
    x := (x_1, ..., x_N)
    r := |x| = sqrt(x * x)

univariate RBF kernel:
    phi(r): R+ --> R;

first & second derivatives of phi:
    phi'(r) : R+ --> R
    phi"(r) : R+ --> R

examples of phi:
    phi(r) := |r|^3                             // cubic
    phi(r) := exp(-r^2/(2 * sigma^2))           // Gaussian

multivariate RBF kernel:
    f(x): R^N --> R
    f(x) := phi(r)
          = phi(|x|)

gradient:
    g(x): R^N --> R^N
    g(x) := df/dx(x)
          = phi'(r) / r * x

definition of alpha:
    alpha(r) := phi'(r) / r
    g(x) = alpha(r) * x

Hessian:
    H(x): R^N --> R^(NxN)
    H(x) := d^2f/dx^2
          = dg/dx
          = {(r * phi"(r) - phi'(r)) / r^3} * (x % x) + {phi'(r) / r} * I

definition of beta:
    beta(r) := (r * phi"(r) - phi'(r)) / r^3
    H(x) = alpha(r) * I + beta(r) * (x % x);

*/

// Gradient is always defined even for singular (r=0) cases --> can be used with Hermite RBF |
//-------------------------------------------------------------------------------------------+
struct RBFKernel_Gaussian {
    struct Param {
        double sigma;
        Param(double sigma = 1) : sigma(sigma) {}
    } param;
    double operator()       (double r) const { return std::exp(-r * r / (2 * param.sigma * param.sigma)); }
    double derivative_first (double r) const { return operator()(r) * r / (-param.sigma * param.sigma); }
    double derivative_second(double r) const { return operator()(r) * (r * r - param.sigma * param.sigma) / (param.sigma * param.sigma * param.sigma * param.sigma); }
    double alpha_singular() const { return -1 / (param.sigma * param.sigma); }
    static const bool AlphaAlwaysDefined = true;
    static const bool Decaying = true;
};

struct RBFKernel_SquaredInverse {
    struct Param {
        double epsilon;
        Param(double epsilon = 0) : epsilon(epsilon) {}
    } param;
    double operator()(double r) const {
        if (r == 0 && param.epsilon == 0)
            return std::numeric_limits<double>::infinity();         // indication of hard constraint (for MLS)
        return 1 / (r * r + param.epsilon * param.epsilon);
    }
    double derivative_first(double r) const {
        double t = operator()(r);
        return -2 * r * t * t;
    }
    double derivative_second(double r) const {
        double t = operator()(r);
        return (6 * r * r - 2 * param.epsilon * param.epsilon) * t * t * t;
    }
    double alpha_singular() const {
        return -2 / (param.epsilon * param.epsilon * param.epsilon * param.epsilon);
    }
    static const bool AlphaAlwaysDefined = true;
    static const bool Decaying = true;
};

struct RBFKernel_Wendland {
    // Wendland, H. 
    // Piecewise polynomial, positive definite and compactly supported radial basis functions of minimal degree. 
    // Advances in Computational Mathematics 4, 389--396, 1995.
    struct Param {
        double support;
        Param(double support = 1) : support(support) {}
    } param;
    double operator()(double r) const {
        double t = r / param.support;
        if (t > 1)
            return 0;
        return (1 - t) * (1 - t) * (1 - t) * (1 - t) * (4 * t + 1);
    }
    double derivative_first(double r) const {
        double t = r / param.support;
        if (t > 1)
            return 0;
        return -20 * t * (1 - t) * (1 - t) * (1 - t) / param.support;
    }
    double derivative_second(double r) const {
        double t = r / param.support;
        if (t > 1)
            return 0;
        return 20 * (1 - t) * (1 - t) * (4 * t - 1) / (param.support * param.support);
    }
    double alpha_singular() const { return -20 / (param.support * param.support); }
    static const bool AlphaAlwaysDefined = true;
    static const bool Decaying = true;
};

struct RBFKernel_Cubed {
    struct Param {} param;
    double operator()       (double r) const { return r * r * r; }
    double derivative_first (double r) const { return 3 * r * r; }
    double derivative_second(double r) const { return 6 * r; }
    double alpha_singular() const { return 0; }
    static const bool AlphaAlwaysDefined = true;
    static const bool Decaying = false;
};

// gradient is undefined for singular (r=0) cases --> cannot be used with HermiteRBF |
//-----------------------------------------------------------------------------------+
struct RBFKernel_Identity {
    struct Param {} param;
    double operator()       (double r) const { return r; }
    double derivative_first (double r) const { return 1; }
    double derivative_second(double r) const { return 0; }
    double alpha_singular() const { throw std::logic_error("undefined!"); }
    static const bool AlphaAlwaysDefined = false;
    static const bool Decaying = false;
};

struct RBFKernel_SquaredLog {
    struct Param {} param;
    double operator()(double r) const {
        if (r == 0)
            return 0;
        return r * r * std::log(r);
    }
    double derivative_first(double r) const {
        if (r == 0)
            return 0;
        return 2 * r * std::log(r) + r;
    }
    double derivative_second(double r) const {
        if (r == 0)
            return -std::numeric_limits<double>::infinity();        // not sure if this treatment is correct...
        return 2 * std::log(r) + 3;
    }
    double alpha_singular() const { throw std::logic_error("undefined!"); }
    static const bool AlphaAlwaysDefined = false;
    static const bool Decaying = false;
};

// general template classes |
//--------------------------+

template <class _RBFKernel_Core>
struct RBFKernel_Univariate {
    typedef _RBFKernel_Core Core;
    
    Core core;
    
    double operator()(double r) const { return core(r); }
    double derivative_first(double r) const { return core.derivative_first(r); }
    double derivative_second(double r) const { return core.derivative_second(r); }
    double alpha(double r) const { return core.derivative_first(r) / r; }
    double beta(double r) const { return (r * core.derivative_second(r) - core.derivative_first(r)) / (r * r * r); }
};

template <int _Dim, class _RBFKernel_Core>
struct RBFKernel_Bivariate {
    enum { Dim  = _Dim };
    
    typedef Eigen::Matrix<double, Dim, 1> Point;
    typedef Eigen::Matrix<double, 1, Dim> Gradient;
    typedef Eigen::Matrix<double, Dim, Dim> Hessian;
    typedef RBFKernel_Univariate<_RBFKernel_Core> Univariate;
    
    static const bool GradientAlwaysDefined = Univariate::Core::AlphaAlwaysDefined;
    static const bool Decaying              = Univariate::Core::Decaying;
    
    Univariate univariate;
    
    double operator()(const Point& p0, const Point& p1) const { return univariate((p1 - p0).norm()); }
    Gradient gradient(const Point& p_variable, const Point& p_fixed) const {
        Point p_diff = p_variable - p_fixed;
        double r = p_diff.norm();
        if (r == 0)
            return Gradient::Zero();
        return univariate.alpha(r) * p_diff.transpose();
    }
    Hessian hessian(const Point& p_variable, const Point& p_fixed) const {
        Point p_diff = p_variable - p_fixed;
        double r = p_diff.norm();
        if (r == 0)
            return univariate.core.alpha_singular() * Hessian::Identity();
        return univariate.alpha(r) * Hessian::Identity() + univariate.beta(r) * p_diff * p_diff.transpose();
    }
    
    // easy access to kernel parameters
    typename Univariate::Core::Param& param() { return univariate.core.param; }
    const typename Univariate::Core::Param& param() const { return univariate.core.param; }
};

}

