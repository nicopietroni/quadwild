#pragma once

#include <GL/glew.h>
#include <vector>
#include <tuple>
#include <Eigen/Core>
#include "../geometry/PointNormal.h"

namespace kt84 {

namespace phong_tessellation {
    enum struct Mode {
        TRIANGLES = 0,
        LINES,
        TRIANGLE_FAN,
        LINE_LOOP,
        LINE_STRIP
    };
    
    inline int& subdiv() {
        static int subdiv = 3;
        return subdiv;
    }
    inline double& weight() {
        static double weight = 0.7;
        return weight;
    }
    
    namespace internal {
        inline Mode& mode() {
            static Mode mode = Mode::TRIANGLES;
            return mode;
        }
        inline bool& enabled() {
            static bool enabled = true;
            return enabled;
        }
        inline std::vector<std::tuple<bool, int, double>>& config_stack() {
            static std::vector<std::tuple<bool, int, double>> config_stack;
            return config_stack;
        }
        inline std::vector<PointNormal>& pointnormals() {
            static std::vector<PointNormal> pointnormals;
            return pointnormals;
        }
        inline PointNormal project(const PointNormal& p, const PointNormal& q) {
            PointNormal result;
            result.head(3) = q.head(3) - (q - p).head(3).dot(p.tail(3)) * p.tail(3);
            result.tail(3) = p.tail(3);
            return result;
        }
        inline PointNormal interpolate(const Eigen::Vector3d& baryCoord, const PointNormal& pn0, const PointNormal& pn1, const PointNormal& pn2) {
            PointNormal result = baryCoord[0] * pn0 + baryCoord[1] * pn1 + baryCoord[2] * pn2;
            pn_normalize(result);
            return result;
        }
        inline PointNormal interpolate(double u, const PointNormal& pn0, const PointNormal& pn1) {
            PointNormal result = (1 - u) * pn0 + u * pn1;
            pn_normalize(result);
            return result;
        }
        inline void draw_triangle_sub(int i, int j, const PointNormal& pn0, const PointNormal& pn1, const PointNormal& pn2) {
            Eigen::Vector3d baryCoord;
            baryCoord[1] = i / static_cast<double>(subdiv());
            baryCoord[2] = j / static_cast<double>(subdiv());
            baryCoord[0] = 1.0 - baryCoord[1] - baryCoord[2];
            PointNormal p  = interpolate(baryCoord, pn0, pn1, pn2);
            PointNormal q0 = project(pn0, p);
            PointNormal q1 = project(pn1, p);
            PointNormal q2 = project(pn2, p);
            PointNormal q  = interpolate(baryCoord, q0, q1, q2);
            PointNormal r = (1.0 - weight()) * p + weight() * q;
            pn_normalize(r);
            glNormal3dv(&r[3]);
            glVertex3dv(&r[0]);
        }
        inline void draw_line_sub(int i, const PointNormal& pn0, const PointNormal& pn1) {
            double u = i / static_cast<double>(subdiv());
            PointNormal p  = interpolate(u, pn0, pn1);
            PointNormal q0 = project(pn0, p);
            PointNormal q1 = project(pn1, p);
            PointNormal q  = interpolate(u, q0, q1);
            PointNormal r = (1.0 - weight()) * p + weight() * q;
            pn_normalize(r);
            glNormal3dv(&r[3]);
            glVertex3dv(&r[0]);
        }
    }
    
    inline void begin(Mode mode_) { internal::mode() = mode_; }
    
    inline void vertex(const PointNormal& pn) { internal::pointnormals().push_back(pn); }
    inline void vertex(const Eigen::Vector3d& point, const Eigen::Vector3d& normal) {
        PointNormal pn;
        pn << point, normal;
        vertex(pn);
    }
    
    inline void enable () { internal::enabled() = true; }
    inline void disable() { internal::enabled() = false; }

    // config stack
    inline void push_config() {
        internal::config_stack().push_back(std::make_tuple(internal::enabled(), subdiv(), weight()));
    }
    inline void pop_config () {
        std::tie(internal::enabled(), subdiv(), weight()) = internal::config_stack().back();
        internal::config_stack().pop_back();
    }
    
    inline void draw_triangle(const PointNormal& pn0, const PointNormal& pn1, const PointNormal& pn2) {
        glBegin(GL_TRIANGLES);
        if (internal::enabled()) {
            for (int i = 0; i < subdiv(); ++i) {
                for (int j = 0; j < subdiv() - 1 - i; ++j) {
                    internal::draw_triangle_sub(i    , j    , pn0, pn1, pn2);
                    internal::draw_triangle_sub(i + 1, j    , pn0, pn1, pn2);
                    internal::draw_triangle_sub(i    , j + 1, pn0, pn1, pn2);
                    internal::draw_triangle_sub(i    , j + 1, pn0, pn1, pn2);
                    internal::draw_triangle_sub(i + 1, j    , pn0, pn1, pn2);
                    internal::draw_triangle_sub(i + 1, j + 1, pn0, pn1, pn2);
                }
                internal::draw_triangle_sub(i    , subdiv()     - i, pn0, pn1, pn2);
                internal::draw_triangle_sub(i    , subdiv() - 1 - i, pn0, pn1, pn2);
                internal::draw_triangle_sub(i + 1, subdiv() - 1 - i, pn0, pn1, pn2);
            }
        } else {
            glNormal3dv(&pn0[3]);    glVertex3dv(&pn0[0]);
            glNormal3dv(&pn1[3]);    glVertex3dv(&pn1[0]);
            glNormal3dv(&pn2[3]);    glVertex3dv(&pn2[0]);
        }
        glEnd();
    }
    
    inline void draw_line(const PointNormal& pn0, const PointNormal& pn1) {
        glBegin(GL_LINE_STRIP);
        if (internal::enabled()) {
            for (int i = 0; i < subdiv(); ++i) {
                internal::draw_line_sub(i    , pn0, pn1);
                internal::draw_line_sub(i + 1, pn0, pn1);
            }
        } else {
            glNormal3dv(&pn0[3]);    glVertex3dv(&pn0[0]);
            glNormal3dv(&pn1[3]);    glVertex3dv(&pn1[0]);
        }
        glEnd();
    }
    
    inline void end() {
        auto& mode = internal::mode();
        auto& pointnormals = internal::pointnormals();
        auto n = pointnormals.size();
        
        if (mode == Mode::TRIANGLES) {
            if (n < 3) return;
            size_t m = n / 3;
            for (size_t i = 0; i < m; ++i)
                draw_triangle(pointnormals[3 * i], pointnormals[3 * i + 1], pointnormals[3 * i + 2]);
        
        } else if (mode == Mode::LINES) {
            if (n < 2) return;
            size_t m = n / 2;
            for (size_t i = 0; i < m; ++i)
                draw_line(pointnormals[2 * i], pointnormals[2 * i + 1]);
        
        } else if (mode == Mode::TRIANGLE_FAN) {
            if (n < 3) return;
            for (size_t i = 1; i < n - 1; ++i)
                draw_triangle(pointnormals[0], pointnormals[i], pointnormals[i + 1]);
        
        } else if (mode == Mode::LINE_LOOP) {
            if (n < 2) return;
            for (size_t i = 0; i < n; ++i)
                draw_line(pointnormals[i], pointnormals[(i + 1) % n]);
        
        } else if (mode == Mode::LINE_STRIP) {
            if (n < 2) return;
            for (size_t i = 0; i < n - 1; ++i)
                draw_line(pointnormals[i], pointnormals[i + 1]);
        }
        
        pointnormals .clear();
    }
};

}

