#pragma once

#include <Eigen/Core>

namespace kt84 {

struct Camera {
    enum class DragMode {
        NONE = 0,
        ROTATE,
        PAN,
        ZOOM,
    };
    
    DragMode drag_mode;
    Eigen::Vector3d center;
    Eigen::Vector2i prev_pos;
    int width;
    int height;
    bool auto_flip_y;
    
    Camera()
        : drag_mode(DragMode::NONE)
        , center(0, 0, 0)
        , prev_pos(0, 0)
        , width (0)
        , height(0)
        , auto_flip_y(true)
    {}
    virtual ~Camera() {};
    
    // abstract member function
    virtual void init(const Eigen::Vector3d& eye_, const Eigen::Vector3d& center_, const Eigen::Vector3d& up_) = 0;
    virtual Eigen::Vector3d get_eye() const = 0;
    virtual Eigen::Vector3d get_up () const = 0;
    virtual void mouse_move(int x, int y) = 0;
    virtual void update_center(const Eigen::Vector3d& center_new) {}
    virtual void snap_to_canonical() {}
    
    // non-abstract member function
    void reshape(int width_, int height_) {
        width  = width_;
        height = height_;
    }
    
    inline Eigen::Vector3d center_to_eye() const { return get_eye() - center; }
    inline Eigen::Vector3d eye_to_center() const { return -center_to_eye(); }
    
    void mouse_down(int x, int y, DragMode drag_mode_) {
        if (auto_flip_y) y = height - y;
        prev_pos << x, y;
        drag_mode = drag_mode_;
    }
    
    void mouse_up() { drag_mode = DragMode::NONE; }
    
};

}
