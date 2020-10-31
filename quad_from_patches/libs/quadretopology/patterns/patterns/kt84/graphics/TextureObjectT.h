#pragma once
#include <GL/glew.h>

namespace kt84 {

struct FramebufferObject;

template <GLenum TTarget = GL_TEXTURE_2D>
struct TextureObjectT {
    static const GLenum Target = TTarget;
    
    GLuint  handle        ;
    GLsizei width         ;
    GLsizei height        ;
    GLenum  format        ;
    GLint   internalformat;
    
    TextureObjectT()
        : handle        (0)
        , width         (0)
        , height        (0)
        , format        (0)
        , internalformat(0)
    {}
    
    void init() { glGenTextures(1, &handle); }
    
    void bind() const { glBindTexture(TTarget, handle); }
    
    static void unbind() { glBindTexture(TTarget, 0); }
    
    static void enable () { glEnable (TTarget); }
    
    static void disable() { glDisable(TTarget); }
    
    void allocate(GLsizei width_, GLsizei height_, GLenum  format_ = GL_RGBA, GLint   internalformat_ = GL_RGBA) {
        width          = width_;
        height         = height_;
        format         = format_;
        internalformat = internalformat_;
        glTexImage2D(TTarget, 0, internalformat, width, height, 0, format, GL_UNSIGNED_BYTE, 0);
    }
    
    void copy_cpu2gpu(GLenum type, const GLvoid* pixels) {
        glTexImage2D(TTarget, 0, internalformat, width, height, 0, format, type, pixels);
    }
    void copy_cpu2gpu(FramebufferObject& framebuffer, const GLvoid* pixels);
    
    void copy_gpu2cpu(GLenum type, GLvoid* pixels) {
        glGetTexImage(TTarget, 0, format, type, pixels);
    }
    void copy_gpu2cpu(FramebufferObject& framebuffer, GLvoid* pixels);
    
    void copy_gpu2gpu(TextureObjectT<TTarget>& dst);
    
    void set_wrap (GLenum param) const {
        glTexParameteri(TTarget, GL_TEXTURE_WRAP_S, param);
        glTexParameteri(TTarget, GL_TEXTURE_WRAP_T, param);
    }
    void set_filter(GLenum param) const {
        glTexParameteri(TTarget, GL_TEXTURE_MIN_FILTER, param);
        glTexParameteri(TTarget, GL_TEXTURE_MAG_FILTER, param);
    }
    void set_env(GLint param) const {
        glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, param);
    }
    void set_default_param() const {
        set_wrap(GL_REPEAT);
        set_filter(GL_LINEAR);
        set_env(GL_REPLACE);
    }
    
    bool is_compatible(const TextureObjectT<TTarget>& rhs) const {
        return
            width          == rhs.width          &&
            height         == rhs.height         &&
            format         == rhs.format         &&
            internalformat == rhs.internalformat;
    }
};

typedef TextureObjectT<GL_TEXTURE_2D> TextureObject;

}
