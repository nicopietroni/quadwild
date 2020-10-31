#pragma once
#include <Eigen/Core>

namespace kt84 {
    template <typename _Scalar, int _Rows, int _Cols>
    struct Matrix11ToScalar {
        using Type = Eigen::Matrix<_Scalar, _Rows, _Cols>;
        static auto Zero() -> decltype(Type::Zero()) { return Type::Zero(); }
        static auto Constant(const _Scalar& value) -> decltype(Type::Constant(value)) { return Type::Constant(value); }
    };
    template <typename _Scalar>
    struct Matrix11ToScalar<_Scalar, 1, 1> {
        using Type = _Scalar;
        static _Scalar Zero() { return 0; }
        static _Scalar Constant(const _Scalar& value) { return value; }
    };
}
