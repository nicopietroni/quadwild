#pragma once

#include "Camera.h"
#include "../eigen_util.h"
#include <Eigen/Geometry>

namespace kt84 {

struct CameraFree : public Camera {
    Eigen::Vector3d eye, up;
    
    CameraFree()
        : eye(Eigen::Vector3d::UnitZ())
        , up(Eigen::Vector3d::UnitY())
    {}
    
    void init(const Eigen::Vector3d& eye_, const Eigen::Vector3d& center_, const Eigen::Vector3d& up_) {
        eye    = eye_   ;
        center = center_;
        up     = up_    ;
    }
    Eigen::Vector3d get_eye() const { return eye; }
    Eigen::Vector3d get_up () const { return up ; }
    
    void mouse_move(int x, int y) {
        if (auto_flip_y) y = height - y;
        const int viewport_size = (width + height) / 2;
        
        if (drag_mode == DragMode::NONE)
            return;
        
        Eigen::Vector2i pos(x, y);
        static const double PI = 2.0 * std::asin(1.0);
        switch (drag_mode) {
        case DragMode::ROTATE:
            {
                double theta_x = (2 * PI * (pos.x() - prev_pos.x())) / viewport_size;
                double theta_y = (2 * PI * (pos.y() - prev_pos.y())) / viewport_size;
                
                Eigen::AngleAxisd rot_hrz(-theta_x, up  );
                auto center_to_eye_new = rot_hrz * center_to_eye();
                
                Eigen::Vector3d left = center_to_eye_new.cross(up).normalized();
                Eigen::AngleAxisd rot_vrt(-theta_y, left);
                center_to_eye_new = rot_vrt * center_to_eye_new;
                up = rot_vrt * up;
                
                eye = center + center_to_eye_new;
            }
            break;
        case DragMode::PAN:
            {
                double len = center_to_eye().norm();
                
                Eigen::Vector3d trans_x(eye_to_center().cross(up));
                Eigen::Vector3d trans_y(up);
                trans_x.normalize();
                trans_y.normalize();
                trans_x *= -len * (pos.x() - prev_pos.x()) / viewport_size;
                trans_y *= -len * (pos.y() - prev_pos.y()) / viewport_size;
                
                eye    += trans_x + trans_y;
                center += trans_x + trans_y;
            }
            break;
        case DragMode::ZOOM:
            {
                auto delta = eye_to_center() * (pos.y() - prev_pos.y()) / static_cast<double>(viewport_size);
                eye -= delta;
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
    }
    void snap_to_canonical() {
        auto d = center_to_eye();
        double dr = d.norm();
        int d_axis = eigen_util::closest_axis(d);
        d[(d_axis + 1) % 3] = d[(d_axis + 2) % 3] = 0;
        d[d_axis] = d[d_axis] < 0 ? -dr : dr;
        eye = center + d;
        int up_axis = eigen_util::closest_axis(up);
        up[(up_axis + 1) % 3] = up[(up_axis + 2) % 3] = 0;
        up[up_axis] = up[up_axis] < 0 ? -1 : 1;
    }
};

}
