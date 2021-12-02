#pragma once
#include <gl/glew.h>
#include "TextureObjectT.h"
#include "RenderbufferObject.h"

namespace kt84 {

struct FramebufferObject {
    GLuint handle;
    
    FramebufferObject()
        : handle(0)
    {}
    
    void init() { glGenFramebuffersEXT(1, &handle); }
    
    void bind() const { glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, handle); }
    
    static void unbind() { glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0); }
    
    bool is_bound() const {
        GLint current_handle;
        glGetIntegerv(GL_FRAMEBUFFER_BINDING_EXT, &current_handle);
        return handle == current_handle;
    }
    
    template <GLenum Target>
    static void attach_texture(GLenum attachment, const TextureObjectT<Target>& texture) {
        glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, attachment, Target, texture.handle, 0);
    }
    template <GLenum Target>
    static void detach_texture(GLenum attachment, const TextureObjectT<Target>& texture) {
        glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, attachment, Target, 0, 0);
    }
    
    static void attach_renderbuffer(GLenum attachment, const RenderbufferObject& renderbuffer) {
        glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, attachment, GL_RENDERBUFFER_EXT, renderbuffer.handle);
    }
    static void detach_renderbuffer(GLenum attachment) {
        glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, attachment, GL_RENDERBUFFER_EXT, 0);
    }
};

}
