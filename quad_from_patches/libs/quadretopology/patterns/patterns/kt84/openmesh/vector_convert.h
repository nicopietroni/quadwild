#pragma once

#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <Eigen/Core>
#include <utility>
#include <vector>

namespace kt84 {

// vector type conversion from Eigen to OpenMesh
template <int N>
inline OpenMesh::VectorT<double, N> e2o (const Eigen::Matrix<double, N, 1, 0, N, 1>& v) {
    OpenMesh::VectorT<double, N> result;
    for (int i = 0; i < N; ++i)
        result[i] = v[i];
    return result;
}
template <> inline OpenMesh::Vec2d e2o(const Eigen::Vector2d& v) { return OpenMesh::Vec2d(v[0], v[1]); }
template <> inline OpenMesh::Vec3d e2o(const Eigen::Vector3d& v) { return OpenMesh::Vec3d(v[0], v[1], v[2]); }
template <> inline OpenMesh::Vec4d e2o(const Eigen::Vector4d& v) { return OpenMesh::Vec4d(v[0], v[1], v[2], v[3]); }
template <int N>    // convert elements in std::pair
inline std::pair<OpenMesh::VectorT<double, N>, OpenMesh::VectorT<double, N>> e2o(const std::pair<Eigen::Matrix<double, N, 1, 0, N, 1>, Eigen::Matrix<double, N, 1, 0, N, 1>>& v_pair) {
    return std::make_pair(e2o(v_pair.first), e2o(v_pair.second));
}
template <int N>    // convert elements in std::vector
inline std::vector<OpenMesh::VectorT<double, N>> e2o(const std::vector<Eigen::Matrix<double, N, 1>>& v_list) {
    std::vector<OpenMesh::VectorT<double, N>> result;
    result.reserve(v_list.size());
    for (auto& v : v_list)
        result.push_back(e2o(v));
    return result;
}

// vector type conversion from OpenMesh to Eigen
template <int N>
inline Eigen::Matrix<double, N, 1, 0, N, 1> o2e(const OpenMesh::VectorT<double, N>& v) {
    Eigen::Matrix<double, N, 1, 0, N, 1> result;
    for (int i = 0; i < N; ++i)
        result[i] = v[i];
    return result;
}
template <> inline Eigen::Vector2d o2e(const OpenMesh::Vec2d& v) { return Eigen::Vector2d(v[0], v[1]); }
template <> inline Eigen::Vector3d o2e(const OpenMesh::Vec3d& v) { return Eigen::Vector3d(v[0], v[1], v[2]); }
template <> inline Eigen::Vector4d o2e(const OpenMesh::Vec4d& v) { return Eigen::Vector4d(v[0], v[1], v[2], v[3]); }
template <int N>    // convert elements in std::pair
inline std::pair<Eigen::Matrix<double, N, 1, 0, N, 1>, Eigen::Matrix<double, N, 1, 0, N, 1>> o2e(const std::pair<OpenMesh::VectorT<double, N>, OpenMesh::VectorT<double, N>>& v_pair) {
    return std::make_pair(o2e(v_pair.first), o2e(v_pair.second));
}
template <int N>    // convert elements in std::vector
inline std::vector<Eigen::Matrix<double, N, 1>> o2e(const std::vector<OpenMesh::VectorT<double, N>>& v_list) {
    std::vector<Eigen::Matrix<double, N, 1>> result;
    result.reserve(v_list.size());
    for (auto& v : v_list)
        result.push_back(o2e(v));
    return result;
}

}
