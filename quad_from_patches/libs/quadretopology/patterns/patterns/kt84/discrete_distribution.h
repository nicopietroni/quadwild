#pragma once
#include <random>
#include <iterator>

// Workaround for the lacking ctor of (Iter, Iter) in the stupid VC11 implementation. Should be unnecessary in future...
namespace kt84 {

template <typename T = int>
struct discrete_distribution : public std::discrete_distribution<T> {
    discrete_distribution() : std::discrete_distribution<T>() {}    // Default ctor
    template<class InputIt>
    discrete_distribution(InputIt first, InputIt last)
        : std::discrete_distribution<T>(std::distance(first, last), 0, std::distance(first, last),
                                        [&](double i){
                                            auto temp = first;
                                            std::advance(temp, static_cast<int>(i));
                                            return *temp;
                                        })
    {}
    template<class UnaryOperation>
    discrete_distribution(std::size_t count, double xmin, double xmax, UnaryOperation unary_op)
        : std::discrete_distribution<T>(count, xmin, xmax, unary_op)
    {}
};

}
