#pragma once
#include <GL/glew.h>
#include <boost/optional.hpp>
#include "../eigen_def.h"

namespace kt84 {

namespace graphics_util {
    enum struct ParamName {
        #if defined(KT84_USE_OPENGL_4)
        #include "ParamName4.inl"
        #elif defined(KT84_USE_OPENGL_3_3)
        #include "ParamName3.3.inl"
        #else
        #include "ParamName2.1.inl"
        #endif
    };
    inline int    glGet1i(ParamName pname) { int    value; glGetIntegerv(static_cast<GLenum>(pname), &value); return value; }
    inline float  glGet1f(ParamName pname) { float  value; glGetFloatv  (static_cast<GLenum>(pname), &value); return value; }
    inline double glGet1d(ParamName pname) { double value; glGetDoublev (static_cast<GLenum>(pname), &value); return value; }
    template <class EigenType> inline EigenType glGetXi(ParamName pname) { EigenType value; glGetIntegerv(static_cast<GLenum>(pname), value.data()); return value; }
    template <class EigenType> inline EigenType glGetXf(ParamName pname) { EigenType value; glGetFloatv  (static_cast<GLenum>(pname), value.data()); return value; }
    template <class EigenType> inline EigenType glGetXd(ParamName pname) { EigenType value; glGetDoublev (static_cast<GLenum>(pname), value.data()); return value; }
    
    inline Eigen::Matrix4d get_modelview_matrix() {
        return glGetXd<Eigen::Matrix4d>(ParamName::MODELVIEW_MATRIX);
    }
    inline Eigen::Matrix4d get_projection_matrix() {
        return glGetXd<Eigen::Matrix4d>(ParamName::PROJECTION_MATRIX);
    }
    inline Eigen::Vector4i get_viewport() {
        return glGetXi<Eigen::Vector4i>(ParamName::VIEWPORT);
    }
    inline float read_depth(int x, int y) {
        float z = 0;
    	glReadPixels(x, y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &z);
        return z;
    }
    inline Eigen::Vector3d unproject(const Eigen::Vector3d& win_xyz,
                                     const Eigen::Matrix4d& modelview_matrix,
                                     const Eigen::Matrix4d& projection_matrix,
                                     const Eigen::Vector4i& viewport)
    {
        Eigen::Vector3d obj_xyz;
    	gluUnProject( win_xyz.x(),  win_xyz.y(),  win_xyz.z(),     modelview_matrix.data(), projection_matrix.data(), viewport.data(),
                     &obj_xyz.x(), &obj_xyz.y(), &obj_xyz.z());
        return obj_xyz;
    }
    inline Eigen::Vector3d unproject(const Eigen::Vector3d& win_xyz) {
        return unproject(win_xyz, get_modelview_matrix(), get_projection_matrix(), get_viewport());
    }
    inline boost::optional<Eigen::Vector3d> read_and_unproject(int x, int y) {
        auto modelview_matrix  = get_modelview_matrix();
        auto projection_matrix = get_projection_matrix();
        auto viewport          = get_viewport();
        
        float z = read_depth(x, y);
        if (z == 1.0f)
            return boost::none;
        
        return unproject(Eigen::Vector3d(x, y, z));
    }
    inline Eigen::Vector3d project(const Eigen::Vector3d& obj_xyz,
                                   const Eigen::Matrix4d& modelview_matrix,
                                   const Eigen::Matrix4d& projection_matrix,
                                   const Eigen::Vector4i& viewport)
    {
        Eigen::Vector3d win_xyz;
    	gluProject( obj_xyz.x(),  obj_xyz.y(),  obj_xyz.z(),    modelview_matrix.data(), projection_matrix.data(), viewport.data(),
                   &win_xyz.x(), &win_xyz.y(), &win_xyz.z());
        return win_xyz;
    }
    inline Eigen::Vector3d project(const Eigen::Vector3d& obj_xyz) {
        return project(obj_xyz, get_modelview_matrix(), get_projection_matrix(), get_viewport());
    }

    inline void glColor3f (const Eigen::Vector3f & color) { ::glColor3fv (&color[0]); }
    inline void glColor3d (const Eigen::Vector3d & color) { ::glColor3dv (&color[0]); }
    //inline void glColor3b (const Eigen::Vector3c & color) { ::glColor3bv (&color[0]); }       // Weird error: cannot convert parameter 1 from 'const char *' to 'const GLbyte *'
    inline void glColor3i (const Eigen::Vector3i & color) { ::glColor3iv (&color[0]); }
    inline void glColor3ub(const Eigen::Vector3uc& color) { ::glColor3ubv(&color[0]); }
    inline void glColor3ui(const Eigen::Vector3ui& color) { ::glColor3uiv(&color[0]); }
    inline void glColor4f (const Eigen::Vector4f & color) { ::glColor4fv (&color[0]); }
    inline void glColor4d (const Eigen::Vector4d & color) { ::glColor4dv (&color[0]); }
    //inline void glColor4b (const Eigen::Vector4c & color) { ::glColor4bv (&color[0]); }       // Weird error: cannot convert parameter 1 from 'const char *' to 'const GLbyte *'
    inline void glColor4i (const Eigen::Vector4i & color) { ::glColor4iv (&color[0]); }
    inline void glColor4ub(const Eigen::Vector4uc& color) { ::glColor4ubv(&color[0]); }
    inline void glColor4ui(const Eigen::Vector4ui& color) { ::glColor4uiv(&color[0]); }
    inline void glColor4f (const Eigen::Vector3f& color_rgb, float  color_alpha) { ::glColor4f(color_rgb[0], color_rgb[1], color_rgb[2], color_alpha); }
    inline void glColor4d (const Eigen::Vector3d& color_rgb, double color_alpha) { ::glColor4d(color_rgb[0], color_rgb[1], color_rgb[2], color_alpha); }
    inline void glVertex2f(const Eigen::Vector2f& vertex) { ::glVertex2fv(&vertex[0]); }
    inline void glVertex2d(const Eigen::Vector2d& vertex) { ::glVertex2dv(&vertex[0]); }
    inline void glVertex3f(const Eigen::Vector3f& vertex) { ::glVertex3fv(&vertex[0]); }
    inline void glVertex3d(const Eigen::Vector3d& vertex) { ::glVertex3dv(&vertex[0]); }
    inline void glNormal3f(const Eigen::Vector3f& normal) { ::glNormal3fv(&normal[0]); }
    inline void glNormal3d(const Eigen::Vector3d& normal) { ::glNormal3dv(&normal[0]); }
    inline void glTexCoord2f(const Eigen::Vector2f& texCoord) { ::glTexCoord2fv(&texCoord[0]); }
    inline void glTexCoord2d(const Eigen::Vector2d& texCoord) { ::glTexCoord2dv(&texCoord[0]); }
    inline void glTexCoord3f(const Eigen::Vector3f& texCoord) { ::glTexCoord3fv(&texCoord[0]); }
    inline void glTexCoord3d(const Eigen::Vector3d& texCoord) { ::glTexCoord3dv(&texCoord[0]); }
    inline void glTexCoord4f(const Eigen::Vector4f& texCoord) { ::glTexCoord4fv(&texCoord[0]); }
    inline void glTexCoord4d(const Eigen::Vector4d& texCoord) { ::glTexCoord4dv(&texCoord[0]); }
    inline void glTexCoord4f(const Eigen::Vector3f& texCoord_str, float  texCoord_q) { ::glTexCoord4f(texCoord_str[0], texCoord_str[1], texCoord_str[2], texCoord_q); }
    inline void glTexCoord4d(const Eigen::Vector3d& texCoord_str, double texCoord_q) { ::glTexCoord4d(texCoord_str[0], texCoord_str[1], texCoord_str[2], texCoord_q); }
    inline void glScalef(const Eigen::Vector3f& v) { ::glScalef(v[0], v[1], v[2]); }
    inline void glScaled(const Eigen::Vector3d& v) { ::glScaled(v[0], v[1], v[2]); }
    inline void glScalef(float  s) { ::glScalef(s, s, s); }
    inline void glScaled(double s) { ::glScaled(s, s, s); }
    inline void glTranslatef(const Eigen::Vector3f& v) { ::glTranslatef(v[0], v[1], v[2]); }
    inline void glTranslated(const Eigen::Vector3d& v) { ::glTranslated(v[0], v[1], v[2]); }
    inline void gluLookAt(const Eigen::Vector3d& eye, const Eigen::Vector3d& center, const Eigen::Vector3d& up) { ::gluLookAt(eye.x(), eye.y(), eye.z(), center.x(), center.y(), center.z(), up.x(), up.y(), up.z()); }
#ifdef OPENMESH_VECTOR_HH
    // for OpenMesh::VectorT
    inline void glColor3f (const OpenMesh::Vec3f& color) { ::glColor3fv (&color[0]); }
    inline void glColor3d (const OpenMesh::Vec3d& color) { ::glColor3dv (&color[0]); }
    inline void glColor3i (const OpenMesh::Vec3i& color) { ::glColor3iv (&color[0]); }
    inline void glColor3ub(const OpenMesh::Vec3uc& color) { ::glColor3ubv(&color[0]); }
    inline void glColor3ui(const OpenMesh::Vec3ui& color) { ::glColor3uiv(&color[0]); }
    inline void glColor4f (const OpenMesh::Vec4f & color) { ::glColor4fv (&color[0]); }
    inline void glColor4d (const OpenMesh::Vec4d & color) { ::glColor4dv (&color[0]); }
    inline void glColor4i (const OpenMesh::Vec4i & color) { ::glColor4iv (&color[0]); }
    inline void glColor4ub(const OpenMesh::Vec4uc& color) { ::glColor4ubv(&color[0]); }
    inline void glColor4ui(const OpenMesh::Vec4ui& color) { ::glColor4uiv(&color[0]); }
    inline void glColor4f (const OpenMesh::Vec3f& color_rgb, float  color_alpha) { ::glColor4f(color_rgb[0], color_rgb[1], color_rgb[2], color_alpha); }
    inline void glColor4d (const OpenMesh::Vec3d& color_rgb, double color_alpha) { ::glColor4d(color_rgb[0], color_rgb[1], color_rgb[2], color_alpha); }
    inline void glVertex2f(const OpenMesh::Vec2f& vertex) { ::glVertex2fv(&vertex[0]); }
    inline void glVertex2d(const OpenMesh::Vec2d& vertex) { ::glVertex2dv(&vertex[0]); }
    inline void glVertex3f(const OpenMesh::Vec3f& vertex) { ::glVertex3fv(&vertex[0]); }
    inline void glVertex3d(const OpenMesh::Vec3d& vertex) { ::glVertex3dv(&vertex[0]); }
    inline void glNormal3f(const OpenMesh::Vec3f& normal) { ::glNormal3fv(&normal[0]); }
    inline void glNormal3d(const OpenMesh::Vec3d& normal) { ::glNormal3dv(&normal[0]); }
    inline void glTexCoord2f(const OpenMesh::Vec2f& texCoord) { ::glTexCoord2fv(&texCoord[0]); }
    inline void glTexCoord2d(const OpenMesh::Vec2d& texCoord) { ::glTexCoord2dv(&texCoord[0]); }
    inline void glTexCoord3f(const OpenMesh::Vec3f& texCoord) { ::glTexCoord3fv(&texCoord[0]); }
    inline void glTexCoord3d(const OpenMesh::Vec3d& texCoord) { ::glTexCoord3dv(&texCoord[0]); }
#endif

}

}
