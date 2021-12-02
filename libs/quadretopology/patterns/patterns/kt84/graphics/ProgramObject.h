#pragma once

#include "ShaderObject.h"

namespace kt84 {

struct ProgramObject {
    GLuint	handle;
    
    ProgramObject()
        : handle(0)
    {}
    void init() {
        handle = glCreateProgramObjectARB();
    }
    void attach(const ShaderObject& s) {
	    glAttachObjectARB(handle, s.handle);
    }
    void link(void) const {
	    glLinkProgramARB(handle);
	    GLint result;
	    glGetObjectParameterivARB(handle, GL_OBJECT_LINK_STATUS_ARB, &result);
	    if (result == GL_FALSE) {
		    int	length;
		    glGetObjectParameterivARB(handle, GL_OBJECT_INFO_LOG_LENGTH_ARB, &length);
		    if (length > 0) {
			    int	l;
                std::string info_log;
                info_log.resize(length);
			    glGetInfoLogARB(handle, length, &l, &info_log[0]);
                std::cerr << info_log << std::endl;
		    }
	    }
    }
    void enable() const { glUseProgramObjectARB(handle); }
    static void disable() { glUseProgramObjectARB(0); }

    GLint get_uniform_location(const std::string& name) const {
        GLint ul = glGetUniformLocationARB(handle, name.c_str());
        if (ul == -1) {
            std::cerr << "no such uniform named " << name << std::endl;
        }
        return ul;
    }

    GLint get_attrib_location(const std::string& name) const {
        GLint al = glGetAttribLocationARB(handle, name.c_str());
        if (al == -1) {
            std::cerr << "no such attrib named " << name << std::endl;
        }
        return al;
    }

#pragma warning(disable: 4244)
    
    // glUniform<N><T>
    template <typename T>
    void set_uniform_1(const std::string& name, T v0) const {
        auto& t = typeid(T);
        GLint location = get_uniform_location(name);
        if      (t == typeid(GLint   )) glUniform1i (location, v0);
        else if (t == typeid(GLuint  )) glUniform1ui(location, v0);
        else if (t == typeid(GLfloat )) glUniform1f (location, v0);
        else if (t == typeid(GLdouble)) glUniform1d (location, v0);
    }
    template <typename T>
    void set_uniform_2(const std::string& name, T v0, T v1) const {
        auto& t = typeid(T);
        GLint location = get_uniform_location(name);
        if      (t == typeid(GLint   )) glUniform2i (location, v0, v1);
        else if (t == typeid(GLuint  )) glUniform2ui(location, v0, v1);
        else if (t == typeid(GLfloat )) glUniform2f (location, v0, v1);
        else if (t == typeid(GLdouble)) glUniform2d (location, v0, v1);
    }
    template <typename T>
    void set_uniform_3(const std::string& name, T v0, T v1, T v2) const {
        auto& t = typeid(T);
        GLint location = get_uniform_location(name);
        if      (t == typeid(GLint   )) glUniform3i (location, v0, v1, v2);
        else if (t == typeid(GLuint  )) glUniform3ui(location, v0, v1, v2);
        else if (t == typeid(GLfloat )) glUniform3f (location, v0, v1, v2);
        else if (t == typeid(GLdouble)) glUniform3d (location, v0, v1, v2);
    }
    template <typename T>
    void set_uniform_4(const std::string& name, T v0, T v1, T v2, T v3) const {
        auto& t = typeid(T);
        GLint location = get_uniform_location(name);
        if      (t == typeid(GLint   )) glUniform4i (location, v0, v1, v2, v3);
        else if (t == typeid(GLuint  )) glUniform4ui(location, v0, v1, v2, v3);
        else if (t == typeid(GLfloat )) glUniform4f (location, v0, v1, v2, v3);
        else if (t == typeid(GLdouble)) glUniform4d (location, v0, v1, v2, v3);
    }
    template <class TVector> void set_uniform_2( const std::string& name, const TVector& v ) const { set_uniform_2(name, v[0], v[1]); }
    template <class TVector> void set_uniform_3( const std::string& name, const TVector& v ) const { set_uniform_3(name, v[0], v[1], v[2]); }
    template <class TVector> void set_uniform_4( const std::string& name, const TVector& v ) const { set_uniform_4(name, v[0], v[1], v[2], v[3]); }
    
    // glUniform<N><T>v
    template <typename T>
    void set_uniform_1v( const std::string& name, GLuint count, const T *v ) const {
        auto& t = typeid(T);
        GLint location = get_uniform_location(name);
        if      (t == typeid(GLint   )) glUniform1iv (location, count, v);
        else if (t == typeid(GLuint  )) glUniform1uiv(location, count, v);
        else if (t == typeid(GLfloat )) glUniform1fv (location, count, v);
        else if (t == typeid(GLdouble)) glUniform1dv (location, count, v);
    }
    template <typename T>
    void set_uniform_2v( const std::string& name, GLuint count, const T *v ) const {
        auto& t = typeid(T);
        GLint location = get_uniform_location(name);
        if      (t == typeid(GLint   )) glUniform2iv (location, count, v);
        else if (t == typeid(GLuint  )) glUniform2uiv(location, count, v);
        else if (t == typeid(GLfloat )) glUniform2fv (location, count, v);
        else if (t == typeid(GLdouble)) glUniform2dv (location, count, v);
    }
    template <typename T>
    void set_uniform_3v( const std::string& name, GLuint count, const T *v ) const {
        auto& t = typeid(T);
        GLint location = get_uniform_location(name);
        if      (t == typeid(GLint   )) glUniform3iv (location, count, v);
        else if (t == typeid(GLuint  )) glUniform3uiv(location, count, v);
        else if (t == typeid(GLfloat )) glUniform3fv (location, count, v);
        else if (t == typeid(GLdouble)) glUniform3dv (location, count, v);
    }
    template <typename T>
    void set_uniform_4v( const std::string& name, GLuint count, const T *v ) const {
        auto& t = typeid(T);
        GLint location = get_uniform_location(name);
        if      (t == typeid(GLint   )) glUniform4iv (location, count, v);
        else if (t == typeid(GLuint  )) glUniform4uiv(location, count, v);
        else if (t == typeid(GLfloat )) glUniform4fv (location, count, v);
        else if (t == typeid(GLdouble)) glUniform4dv (location, count, v);
    }
    
    // glUniformMatrix<N><T>v
    template <typename T>
    void set_uniform_matrix_2v( const std::string& name, GLuint count, GLboolean transpose, const T *v ) const {
        auto& t = typeid(T);
        GLint location = get_uniform_location(name);
        if      (t == typeid(GLfloat )) glUniformMatrix2fv(location, count, transpose, v);
        else if (t == typeid(GLdouble)) glUniformMatrix2dv(location, count, transpose, v);
    }
    template <typename T>
    void set_uniform_matrix_3v( const std::string& name, GLuint count, GLboolean transpose, const T *v ) const {
        auto& t = typeid(T);
        GLint location = get_uniform_location(name);
        if      (t == typeid(GLfloat )) glUniformMatrix3fv(location, count, transpose, v);
        else if (t == typeid(GLdouble)) glUniformMatrix3dv(location, count, transpose, v);
    }
    template <typename T>
    void set_uniform_matrix_4v( const std::string& name, GLuint count, GLboolean transpose, const T *v ) const {
        auto& t = typeid(T);
        GLint location = get_uniform_location(name);
        if      (t == typeid(GLfloat )) glUniformMatrix4fv(location, count, transpose, v);
        else if (t == typeid(GLdouble)) glUniformMatrix4dv(location, count, transpose, v);
    }
    template <class TMatrix> void set_uniform_matrix_2( const std::string& name, const TMatrix& m) const { set_uniform_matrix_2v( name, 1, GL_TRUE, &m(0, 0)); }
    template <class TMatrix> void set_uniform_matrix_3( const std::string& name, const TMatrix& m) const { set_uniform_matrix_3v( name, 1, GL_TRUE, &m(0, 0)); }
    template <class TMatrix> void set_uniform_matrix_4( const std::string& name, const TMatrix& m) const { set_uniform_matrix_4v( name, 1, GL_TRUE, &m(0, 0)); }
    
    // glVertexattrib<N><T>
    template <typename T>
    void set_attrib_1(const std::string& name, T v0) const {
        auto& t = typeid(T);
        GLint location = get_attrib_location(name);
        if      (t == typeid(GLfloat )) glVertexAttrib1f(location, v0);
        else if (t == typeid(GLdouble)) glVertexAttrib1d(location, v0);
        else if (t == typeid(GLshort )) glVertexAttrib1s(location, v0);
    }
    template <typename T>
    void set_attrib_2(const std::string& name, T v0, T v1) const {
        auto& t = typeid(T);
        GLint location = get_attrib_location(name);
        if      (t == typeid(GLfloat )) glVertexAttrib2f(location, v0, v1);
        else if (t == typeid(GLdouble)) glVertexAttrib2d(location, v0, v1);
        else if (t == typeid(GLshort )) glVertexAttrib2s(location, v0, v1);
    }
    template <typename T>
    void set_attrib_3(const std::string& name, T v0, T v1, T v2) const {
        auto& t = typeid(T);
        GLint location = get_attrib_location(name);
        if      (t == typeid(GLfloat )) glVertexAttrib3f(location, v0, v1, v2);
        else if (t == typeid(GLdouble)) glVertexAttrib3d(location, v0, v1, v2);
        else if (t == typeid(GLshort )) glVertexAttrib3s(location, v0, v1, v2);
    }
    template <typename T>
    void set_attrib_4(const std::string& name, T v0, T v1, T v2, T v3) const {
        auto& t = typeid(T);
        GLint location = get_attrib_location(name);
        if      (t == typeid(GLfloat )) glVertexAttrib4f(location, v0, v1, v2, v3);
        else if (t == typeid(GLdouble)) glVertexAttrib4d(location, v0, v1, v2, v3);
        else if (t == typeid(GLshort )) glVertexAttrib4s(location, v0, v1, v2, v3);
    }
    template <class TVector> void set_attrib_2( const std::string& name, const TVector& v ) const { set_attrib_2( name, v[0], v[1] ); }
    template <class TVector> void set_attrib_3( const std::string& name, const TVector& v ) const { set_attrib_3( name, v[0], v[1], v[2] ); }
    template <class TVector> void set_attrib_4( const std::string& name, const TVector& v ) const { set_attrib_4( name, v[0], v[1], v[2], v[3] ); }
    
    // glVertexattrib<N><T>v
    template <typename T>
    void set_attrib_1v( const std::string& name, const T *v ) const {
        auto& t = typeid(T);
        GLint location = get_attrib_location(name);
        if      (t == typeid(GLfloat )) glVertexAttrib1fv(location, v);
        else if (t == typeid(GLdouble)) glVertexAttrib1dv(location, v);
        else if (t == typeid(GLshort )) glVertexAttrib1sv(location, v);
    }
    template <typename T>
    void set_attrib_2v( const std::string& name, const T *v ) const {
        auto& t = typeid(T);
        GLint location = get_attrib_location(name);
        if      (t == typeid(GLfloat )) glVertexAttrib2fv(location, v);
        else if (t == typeid(GLdouble)) glVertexAttrib2dv(location, v);
        else if (t == typeid(GLshort )) glVertexAttrib2sv(location, v);
    }
    template <typename T>
    void set_attrib_3v( const std::string& name, const T *v ) const {
        auto& t = typeid(T);
        GLint location = get_attrib_location(name);
        if      (t == typeid(GLfloat )) glVertexAttrib3fv(location, v);
        else if (t == typeid(GLdouble)) glVertexAttrib3dv(location, v);
        else if (t == typeid(GLshort )) glVertexAttrib3sv(location, v);
    }
    template <typename T>
    void set_attrib_4v( const std::string& name, const T *v ) const {
        auto& t = typeid(T);
        GLint location = get_attrib_location(name);
        if      (t == typeid(GLfloat )) glVertexAttrib4fv(location, v);
        else if (t == typeid(GLdouble)) glVertexAttrib4dv(location, v);
        else if (t == typeid(GLshort )) glVertexAttrib4sv(location, v);
    }
#pragma warning(default: 4244)
};

}

