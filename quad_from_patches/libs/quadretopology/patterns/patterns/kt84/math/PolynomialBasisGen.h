#pragma once

#include "BinomialCoefficient.h"
#include <Eigen/Core>

namespace kt84 {

template <int _DimIn, int _Degree>
struct PolynomialBasisGenT {
    enum { DimOut = PolynomialBasisGenT<_DimIn, _Degree - 1>::DimOut + BinomialCoefficient<_DimIn + _Degree - 1, _Degree>::Value };
    typedef Eigen::Matrix<double, _DimIn , 1> Point;
    typedef Eigen::Matrix<double, DimOut, 1> Basis;
    typedef Eigen::Matrix<double, DimOut, _DimIn> Gradient;
    
    static Basis basis(const Point& x) {
        const int N = BinomialCoefficient<_DimIn + _Degree - 1, _Degree>::Value;
        Eigen::Matrix<double, N, 1> result_partial;
        int index_out = 0;
        int index_in[_Degree];
        for (int i = 0; i < _Degree; ++i)
            index_in[i] = 0;
        while (true) {
            double d = 1;
            for (int i = 0; i < _Degree; ++i)
                d *= x[index_in[i]];
            result_partial[index_out++] = d;
            int is_complete = true;
            for (int i = _Degree - 1; i >= 0; --i) {
                if (index_in[i] == _DimIn - 1)
                    continue;
                ++index_in[i];
                for (int j = i + 1; j < _Degree; ++j)
                    index_in[j] = index_in[i];
                is_complete = false;
                break;
            }
            if (is_complete)
                break;
        }
        Basis result;
        result << PolynomialBasisGenT<_DimIn, _Degree - 1>::basis(x), result_partial;
        return result;
    }
    static Gradient gradient(const Point& x) {
        const int N = BinomialCoefficient<_DimIn + _Degree - 1, _Degree>::Value;
        Eigen::Matrix<double, N, _DimIn> result_partial;
        int index_out = 0;
        int index_in[_Degree];
        for (int i = 0; i < _Degree; ++i)
            index_in[i] = 0;
        while (true) {
            for (int i = 0; i < _DimIn; ++i) {
                int cnt = 0;
                double d = 1;
                for (int j = 0; j < _Degree; ++j) {
                    if (index_in[j] == i) {
                        ++cnt;
                        continue;
                    }
                    d *= x[index_in[j]];
                }
                result_partial(index_out, i) = cnt == 0 ? 0 : cnt * d * std::pow(x[i], cnt - 1);
            }
            ++index_out;
            int is_complete = true;
            for (int i = _Degree - 1; i >= 0; --i) {
                if (index_in[i] == _DimIn - 1)
                    continue;
                ++index_in[i];
                for (int j = i + 1; j < _Degree; ++j)
                    index_in[j] = index_in[i];
                is_complete = false;
                break;
            }
            if (is_complete)
                break;
        }
        Gradient result;
        result << PolynomialBasisGenT<_DimIn, _Degree - 1>::gradient(x), result_partial;
        return result;
    }
};

template <int _DimIn>
struct PolynomialBasisGenT<_DimIn, 0> {
    enum { DimOut = 1 };
    typedef Eigen::Matrix<double, _DimIn , 1> Point;
    typedef Eigen::Matrix<double, DimOut, 1> Basis;
    typedef Eigen::Matrix<double, DimOut, _DimIn> Gradient;
    
    static Basis basis(const Point& x) { return Basis::Constant(1); }
    static Gradient gradient(const Point& x) { return Gradient::Zero(); }
};

}

