#pragma once

#include <Eigen/Core>

namespace kt84 {

typedef Eigen::Matrix<double, 6, 1> PointNormal;
inline double pn_norm(const PointNormal& pn) { return pn.head(3).norm(); }
inline void   pn_normalize(PointNormal& pn) { pn.tail(3).normalize(); }
inline PointNormal pn_normalized(const PointNormal& pn) { auto temp = pn; temp.tail(3).normalize(); return temp; }

}
