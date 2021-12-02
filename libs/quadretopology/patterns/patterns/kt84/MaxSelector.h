#pragma once
#include <limits>

namespace kt84 {

template <typename Value, typename Score = double>
struct MaxSelector {
    Score score;
    Value value;
    MaxSelector(const Value& init_value = Value())
        : score(-std::numeric_limits<Score>::max())
        , value(init_value)
    {}
    bool update(Score score_, const Value& value_) {
        if (score_ > score) {
            score = score_;
            value = value_;
            return true;
        }
        return false;
    }
};

}
