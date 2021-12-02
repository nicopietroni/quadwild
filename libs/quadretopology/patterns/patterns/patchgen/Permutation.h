#pragma once
#include <Eigen/Core>

namespace patchgen {

    struct Permutation {
        int num_sides;  // number of sides
        int id;         // 0 <= id <=   num_sides-1 : permutation with original order and offset id
        // n <= id <= 2*num_sides-1 : permutation with reversed order and offset id-num_sides

        Permutation() {
            init(0);
        }
        void init(int num_sides_, int id_ = 0) {
            num_sides = num_sides_;
            id = id_;
        }
        bool is_flipped() const {
            return id >= num_sides;
        }
        bool operator==(const Permutation& rhs) const {
            return num_sides == rhs.num_sides && id == rhs.id;
        }
        // loop utility
        void next(bool is_forward = true) {
            id += is_forward ? 1 : -1;
        }
        bool is_valid() const {
            return 0 <= id && id < 2 * num_sides;
        }
        // accessor
        int operator[](int index) const {
            return id < num_sides ? ((index + id) % num_sides)
                : (num_sides - 1 - (index + id - num_sides) % num_sides);
        }
        Eigen::VectorXi operator()(const Eigen::VectorXi& l) const {
            Eigen::VectorXi l_permuted(num_sides);
            for (int i = 0; i < num_sides; ++i)
                l_permuted[i] = l[this->operator[](i)];
            return l_permuted;
        }
    };

}