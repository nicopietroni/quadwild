#pragma once
#include <gl/glew.h>

namespace kt84 {

struct RenderbufferObject {
    GLuint  handle        ;
    GLsizei width         ;
    GLsizei height        ;
    GLenum  internalformat;
    
    RenderbufferObject()
        : handle        (0)
        , width         (0)
        , height        (0)
        , internalformat(0)
    {}
    
    void init() { glGenRenderbuffersEXT(1, &handle); }
    
    void bind() { glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, handle); }
    
    static void unbind() { glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0); }
    
    void allocate(GLsizei width_, GLsizei height_, GLenum internalformat_) {
        width          = width_ ;
        height         = height_;
        internalformat = internalformat_;
        glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, internalformat, width, height);
    }
};

}
