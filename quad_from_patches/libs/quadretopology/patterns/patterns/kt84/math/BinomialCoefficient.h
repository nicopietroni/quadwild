#pragma once

namespace kt84 {

template <int N, int K>
struct BinomialCoefficient {
    static const int Value = BinomialCoefficient<N, K - 1>::Value * (N - K + 1) / K;
};

template <int N>
struct BinomialCoefficient<N, 0> {
    static const int Value = 1;
};

}

