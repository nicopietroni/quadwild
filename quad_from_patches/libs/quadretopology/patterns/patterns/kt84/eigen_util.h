#pragma once
#include <cmath>
#include <algorithm>
#include <limits>
#include <Eigen/Geometry>
#include <boost/optional.hpp>
#include "BaryCoordT.h"
#include "util.h"
#include "RangeAdaptor.h"

namespace kt84 {
    namespace eigen_util {
        template <typename TMatrix>
        inline auto range_elements(TMatrix& m) -> RangeAdaptor<decltype(m.data())> {
            return make_RangeAdaptor(m.data(), m.data() + m.size());
        }
        template <typename TVector, typename TScalar>
        inline void push_back(TVector& v, TScalar s) {
            int n = v.size();
            v.conservativeResize(n + 1);
            v[n] = s;
        }
        template <typename TVector>
        inline void erase_at(TVector& v, int index) {
            int n = v.size();
            auto tmp(v);
            v.resize(n - 1);
            v << tmp.head(index), tmp.tail(n - 1 - index);
        }
        template <typename TVector>
        inline void swap_xy(TVector& v) {
            std::swap(v.x(), v.y());
        }
        template <typename TAlignedBox>
        inline void bbox_add_margin(TAlignedBox& bbox, double ratio, bool is_margin_square = false) {
            if (is_margin_square) {
                typedef typename TAlignedBox::VectorType VectorType;
                VectorType margin = VectorType::Ones();
                margin *= bbox.diagonal().norm() * ratio;
                bbox.max() += margin;
                bbox.min() -= margin;
            } else {
                bbox.extend((1 + ratio) * bbox.max() - ratio * bbox.min());
                bbox.extend((1 + ratio) * bbox.min() - ratio * bbox.max());
            }
        }
        inline Eigen::Vector2d bbox_bilinear(const Eigen::AlignedBox2d& bbox, double tx, double ty) {            // (0, 0) corresponds to bbox.min, (1, 1) corresponds to bbox.max
            double x = (1 - tx) * bbox.min().x() + tx * bbox.max().x();
            double y = (1 - ty) * bbox.min().y() + ty * bbox.max().x();
            return Eigen::Vector2d(x, y);
        }
        inline Eigen::Vector3d bbox_trilinear(const Eigen::AlignedBox3d& bbox, double tx, double ty, double tz) {            // (0, 0, 0) corresponds to bbox.min, (1, 1, 1) corresponds to bbox.max
            double x = (1 - tx) * bbox.min().x() + tx * bbox.max().x();
            double y = (1 - ty) * bbox.min().y() + ty * bbox.max().x();
            double z = (1 - tz) * bbox.min().z() + tz * bbox.max().z();
            return Eigen::Vector3d(x, y, z);
        }
        inline Eigen::Vector3d orientation_color(Eigen::Vector3d d) {
            auto cx = d.x() > 0 ? Eigen::Vector3d(1, 0, 0) : Eigen::Vector3d(0, 1, 1);
            auto cy = d.y() > 0 ? Eigen::Vector3d(0, 1, 0) : Eigen::Vector3d(1, 0, 1);
            auto cz = d.z() > 0 ? Eigen::Vector3d(0, 0, 1) : Eigen::Vector3d(1, 1, 0);
            d /= std::abs<double>(d.sum());
            d = d.cwiseAbs();
            return d.x() * cx + d.y() * cy + d.z() * cz;
        }
        inline Eigen::Vector3d heat_color(double t) {
            // t     | 0    | 0.25 | 0.5   | 0.75   | 1   |
            // color | blue | cyan | green | yellow | red |
            t = util::clamp(t, 0., 1.);
            Eigen::Vector3d colors[5] = {
                Eigen::Vector3d(0, 0, 1),
                Eigen::Vector3d(0, 1, 1),
                Eigen::Vector3d(0, 1, 0),
                Eigen::Vector3d(1, 1, 0),
                Eigen::Vector3d(1, 0, 0)
            };
            int i = t < 0.25 ? 0 : t < 0.5 ? 1 : t < 0.75 ? 2 : 3;
            double s = (t - i * 0.25) * 4;
            return (1 - s) * colors[i] + s * colors[i + 1];
        }
        template <int N>
        inline int closest_axis(const Eigen::Matrix<double, N, 1, 0, N, 1>& v) {
            int result = -1;
            double v_abs_max = 0;
            for (int i = 0; i < N; ++i) {
                double v_abs = std::abs(v[i]);
                if (v_abs_max < v_abs) {
                    v_abs_max = v_abs;
                    result = i;
                }
            }
            return result;
        }
        inline Eigen::Vector2d compute_gradient(const Eigen::Vector2d& x0, const Eigen::Vector2d& x1, const Eigen::Vector2d& x2, double y0, double y1, double y2) {
            /*
                a.x0 + b = y0
                a.x1 + b = y1
                a.x2 + b = y2
                -->
                a.(x1 - x0) = y1 - y0
                a.(x2 - x0) = y2 - y0
                -->
                |(x1 - x0)^T| * a = |y1 - y0|
                |(x2 - x0)^T|       |y2 - y0|
            */
            Eigen::Matrix2d A;
            A << Eigen::RowVector2d(x1 - x0),
                 Eigen::RowVector2d(x2 - x0);
            return A.inverse() * Eigen::Vector2d(y1 - y0, y2 - y0);
        }
        inline Eigen::Vector3d compute_gradient(const Eigen::Vector3d& x0, const Eigen::Vector3d& x1, const Eigen::Vector3d& x2, const Eigen::Vector3d& x3, double y0, double y1, double y2, double y3) {
            /*
                a.x0 + b = y0
                a.x1 + b = y1
                a.x2 + b = y2
                a.x3 + b = y3
                -->
                a.(x1 - x0) = y1 - y0
                a.(x2 - x0) = y2 - y0
                a.(x3 - x0) = y3 - y0
                -->
                |(x1 - x0)^T|       |y1 - y0|
                |(x2 - x0)^T| * a = |y2 - y0|
                |(x3 - x0)^T|       |y3 - y0|
            */
            Eigen::Matrix3d A;
            A << Eigen::RowVector3d(x1 - x0),
                 Eigen::RowVector3d(x2 - x0),
                 Eigen::RowVector3d(x3 - x0);
            return A.inverse() * Eigen::Vector3d(y1 - y0, y2 - y0, y3 - y0);
        }
        inline Eigen::Vector3d compute_gradient(const Eigen::Vector3d& x0, const Eigen::Vector3d& x1, const Eigen::Vector3d& x2, double y0, double y1, double y2) {
            /*
                Compute gradient restricted to the tangent vectors on the triangle x0-x1-x2.
                a.x0 + b = y0
                a.x1 + b = y1
                a.x2 + b = y2
                a.n      = 0                (n: normal)
                -->
                a.(x1 - x0) = y1 - y0
                a.(x2 - x0) = y2 - y0
                a.n         = 0
                -->
                |(x1 - x0)^T|       |y1 - y0|
                |(x2 - x0)^T| * a = |y2 - y0|
                | n       ^T|       |0      |
            */
            Eigen::Matrix3d A;
            A << Eigen::RowVector3d(x1 - x0),
                 Eigen::RowVector3d(x2 - x0),
                 Eigen::RowVector3d((x1 - x0).cross(x2 - x0));
            return A.inverse() * Eigen::Vector3d(y1 - y0, y2 - y0, 0);
        }
        inline Eigen::Vector2d rotate90(const Eigen::Vector2d& xy) { return Eigen::Vector2d(-xy[1], xy[0]); }
        inline double angle(const Eigen::Vector2d& d0, const Eigen::Vector2d& d1) { return std::atan2(rotate90(d0).dot(d1), d0.dot(d1)); }
        inline double angle(const Eigen::Vector3d& d0, const Eigen::Vector3d& d1) { return std::acos(d0.normalized().dot(d1.normalized())); }
        template <class T>
        inline void orthonormalize(const T& unit, T& p) {
            p -= unit.dot(p) * unit;
            p.normalize();
        }
        template <class T>
        inline bool project_to_line(const T& line_v0, const T& line_v1, const T& point, Eigen::Vector2d& t) {
            /*
                x0 := line_v0
                x1 := line_v1
                y := point
                compute t (which sums up to one) such that
                    | t[0] * x0 + t[1] * x1 - y |^2
                is minimized.
                ---------------------
                u := t[1]
                (1 - u) * x0 + u * x1 =~ y
                u =~ (y - x0).dot(x1 - x0) / (x1 - x0).squaredNorm()
            */
            
            double r = (line_v1 - line_v0).squaredNorm();
            if (r == 0)
                // degenerate
                return false;
            
            t[1] = (point - line_v0).dot(line_v1 - line_v0) / r;
            t[0] = 1 - t[1];
            
            return true;
        }
        
        template <class T>
        inline bool project_to_triangle(const T& triangle_v0, const T& triangle_v1, const T& triangle_v2, const T& point, Eigen::Vector3d& t) {
            /*
                x0 := triangle_v0
                x1 := triangle_v1
                x2 := triangle_v2
                y := point
                compute t (which sums up to one) such that
                    | t[0] * x0 + t[1] * x1 + t[2] * x2 - y |^2
                is minimized.
                ---------------------
                u := t[1]
                v := t[2]
                (1 - u - v) * x0 + u * x1 + v * x2 =~ y
                (x1-x0, x2-x0) * |u| =~ y-x0
                                 |v|
                |u| =~ |(x1-x0).squaredNorm(), (x1-x0).dot(x2-x0)   |^-1 * | (x1-x0).dot(y-x0) |
                |v|    |(x1-x0).dot(x2-x0)   , (x2-x0).squaredNorm()|      | (x2-x0).dot(y-x0) |
            */
            
            T d01 = triangle_v1 - triangle_v0;
            T d02 = triangle_v2 - triangle_v0;
            T d0p = point       - triangle_v0;
            Eigen::Matrix2d M;
            M <<
                d01.squaredNorm(), d01.dot(d02),
                d01.dot(d02)     , d02.squaredNorm();
            if (M.determinant() == 0)
                // degenerate
                return false;
            
            Eigen::Vector2d b;
            b <<
                d01.dot(d0p),
                d02.dot(d0p);
            
            Eigen::Vector2d uv = ((Eigen::Matrix2d)M.inverse()) * b;
            
            t[0] = 1 - uv[0] - uv[1];
            t[1] = uv[0];
            t[2] = uv[1];
            
            return true;
        }
        
        template <class T>
        inline boost::optional<double> distance_to_line(const T& line_v0, const T& line_v1, const T& point, bool do_clamp = false) {
            Eigen::Vector2d t;
            if (!project_to_line(line_v0, line_v1, point, t))
                // degenrate case
                return boost::none;
            
            if (do_clamp) {
                t[0] = util::clamp(t[0], 0.0, 1.0);
                t[1] = 1 - t[0];
            }
            
            return (t[0] * line_v0 + t[1] * line_v1 - point).norm();
        }
        
        template <class T>
        inline boost::optional<double> distance_to_triangle(const T& triangle_v0, const T& triangle_v1, const T& triangle_v2, const T& point, bool do_clamp = false) {
            Eigen::Vector3d t;
            if (!project_to_triangle(triangle_v0, triangle_v1, triangle_v2, point, t))
                // degenrate case
                return boost::none;
            
            if (do_clamp && t[0] < 0 || t[1] < 0 || t[2] < 0) {
                double d0 = *distance_to_line(triangle_v0, triangle_v1, point, true);
                double d1 = *distance_to_line(triangle_v1, triangle_v2, point, true);
                double d2 = *distance_to_line(triangle_v2, triangle_v0, point, true);
                return util::min(d0, d1, d2);
            }
            
            return (t[0] * triangle_v0 + t[1] * triangle_v1 + t[2] * triangle_v2 - point).norm();
        }
        inline double triangle_area(const Eigen::Vector2d& v0, const Eigen::Vector2d& v1, const Eigen::Vector2d& v2) {
            Eigen::Vector2d d1 = v1 - v0;
            Eigen::Vector2d d2 = v2 - v0;
            return d1.x() * d2.y() - d1.y() * d2.x();
        }
        inline double triangle_area(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) {
            Eigen::Vector3d d1 = v1 - v0;
            Eigen::Vector3d d2 = v2 - v0;
            return d1.cross(d2).norm() / 2;
        }
        template <class T>
        inline bool intersection(const T& line0_p0, const T& line0_p1, const T& line1_p0, const T& line1_p1, double& line0_coordinate, double& line1_coordinate) {
            /*
                notation:
                    v0 := line0_p0
                    v1 := line0_p1
                    w0 := line1_p0
                    w1 := line1_p1
                    s  := line0_coordinate
                    t  := line1_coordinate
                seek for (s, t) which minimizes:
                    |(1 - s) * v0 + s * v1 - (1 - t) * w0 - t * w1|^2
                least square sense:
                    | v1-v0, -w1+w0 | * |s| =~ -v0+w0
                                        |t|
            */
            T d0 = line0_p1 - line0_p0;
            T d1 = line1_p1 - line1_p0;
            T e  = line1_p0 - line0_p0;
            double d0d0 = d0.squaredNorm();
            double d0d1 = d0.dot(d1);
            double d1d1 = d1.squaredNorm();
            Eigen::Matrix2d A;
            A <<
                d0d0, -d0d1,
                -d0d1, d1d1;
            
            if (A.determinant() == 0)
                // two lines are parallel
                return false;
            
            Eigen::Vector2d b(d0.dot(e), -d1.dot(e));
            Eigen::Vector2d st = ((Eigen::Matrix2d)A.inverse()) * b;
            
            line0_coordinate = st[0];
            line1_coordinate = st[1];
            return true;
        }
        template <typename T>
        inline std::vector<T> eigen_vectorx_to_std_vector(const Eigen::Matrix<T, -1, 1>& eigen_vectorx) {
            int n = eigen_vectorx.rows();
            std::vector<T> result(n);
            for (int i = 0; i < n; ++i)
                result[i] = eigen_vectorx[i];
            return result;
        }
        template <typename T>
        inline Eigen::Matrix<T, -1, 1> std_vector_to_eigen_vectorx(const std::vector<T>& std_vector) {
            int n = std_vector.size();
            Eigen::Matrix<T, -1, 1> result = Eigen::Matrix<T, -1, 1>::Zero(n);
            for (int i = 0; i < n; ++i)
                result[i] = std_vector[i];
            return result;
        }
        template <typename TVector>
        inline int max_axis(const TVector& v) {
            int index_max = -1;
            typename TVector::Scalar value_max = -std::numeric_limits<typename TVector::Scalar>::max();
            for (int i = 0; i < v.size(); ++i) {
                if (value_max < v[i]) {
                    value_max = v[i];
                    index_max = i;
                }
            }
            return index_max;
        }
        template <typename TVector>
        inline int min_axis(const TVector& v) {
            int index_min = -1;
            typename TVector::Scalar value_min = std::numeric_limits<typename TVector::Scalar>::max();
            for (int i = 0; i < v.size(); ++i) {
                if (value_min < v[i]) {
                    value_min = v[i];
                    index_min = i;
                }
            }
            return index_min;
        }
        template <typename TVector>
        inline void rotate(TVector& v) {
            TVector tmp(v);
            int n = v.size();
            for (int i = 0; i < n; ++i)
                v[i] = tmp[(i + 1) % n];
        }
        template <typename TVector>
        inline void reverse(TVector& v) {
            TVector tmp(v);
            int n = v.size();
            for (int i = 0; i < n; ++i)
                v[i] = tmp[n - 1 - i];
        }
    }
}
