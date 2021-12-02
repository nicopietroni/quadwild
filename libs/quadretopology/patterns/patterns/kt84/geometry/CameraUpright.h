#pragma once

#include "Camera.h"
#include "../util.h"
#include <Eigen/Geometry>

namespace kt84 {

struct CameraUpright : public Camera {
    double r, theta, phi;
    
    CameraUpright()
        : r(1)
        , theta(0)
        , phi(0)
    {}
    
    void init(const Eigen::Vector3d& eye, const Eigen::Vector3d& center_, const Eigen::Vector3d& up) {
        center = center_;
        Eigen::Vector3d v = eye - center;
        r = v.norm();
        double r_xz = std::sqrt(v.x() * v.x() + v.z() * v.z());
        
        if (r_xz == 0) {
            bool is_y_positive = v.y() > 0;
            theta = std::atan2(up.x(), up.z()) + (is_y_positive ? util::pi() : 0);
            phi   = (is_y_positive ? 0.5 : -0.5) * util::pi();
        } else {
            theta = std::atan2(v.x(), v.z());
            phi   = std::atan(v.y() / r_xz);
        }
        
        if (up.y() < 0) {
            theta += util::pi();
            phi    = util::pi() - phi;
        }
    }
    
    Eigen::Vector3d get_up() const {
        return Eigen::Vector3d(
            -sin(theta) * sin(phi),
            cos(phi),
            -cos(theta) * sin(phi));
    }
    
    Eigen::Vector3d get_eye() const {
        return center + r * Eigen::Vector3d(
            sin(theta) * cos(phi),
            sin(phi),
            cos(theta) * cos(phi));
    }
    
    void mouse_move(int x, int y) {
        if (auto_flip_y) y = height - y;
        const int viewport_size = (width + height) / 2;
        
        if (drag_mode == DragMode::NONE)
            return;
        
        Eigen::Vector2i pos(x, y);
        switch (drag_mode) {
        case DragMode::ROTATE:
            {
                theta -= (2 * util::pi() * (pos.x() - prev_pos.x())) / viewport_size;
                phi   += (2 * util::pi() * (pos.y() - prev_pos.y())) / viewport_size;
            }
            break;
        case DragMode::PAN:
            {
                Eigen::Vector3d right(cos(theta), 0, -sin(theta));
                center -= (r * (pos.x() - prev_pos.x()) / viewport_size) * right;
                center += (r * (pos.y() - prev_pos.y()) / viewport_size) * get_up();
            }
            break;
        case DragMode::ZOOM:
            {
                r *= 1 - (pos.x() - prev_pos.x() + pos.y() - prev_pos.y()) / static_cast<double>(viewport_size);
            }
            break;
        }
        prev_pos = pos;
    }
    
    void update_center(const Eigen::Vector3d& center_new) {
        // find the point closest to center_new between center and eye. update center and r accordingly.
        // center + t * center_to_eye = center_new
        auto d = center_to_eye();
        double t = d.dot(center_new - center) / d.squaredNorm();
        center += t * d;
        r *= 1 - t;
    }
    void snap_to_canonical() {
        auto snap_to_pi_2 = [] (double& radian) {
            while (radian < 0)
                radian += 2 * util::pi();
            
            double t = radian / (0.5 * util::pi());
            int n = static_cast<int>(t);
            if (t - n > 0.5)
                ++n;
            
            radian = n * 0.5 * util::pi();
        };
        
        snap_to_pi_2(theta);
        snap_to_pi_2(phi  );
    }
};

}
