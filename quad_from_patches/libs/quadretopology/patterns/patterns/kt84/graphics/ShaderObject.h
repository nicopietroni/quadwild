#pragma once

#include <GL/glew.h>
#include <iostream>
#include <fstream>
#include <sstream>

namespace kt84 {

struct ShaderObject {
    enum class Type {
        FRAGMENT_SHADER,
        VERTEX_SHADER  ,
        GEOMETRY_SHADER
    };
    
    GLuint	handle;
    std::string filename;
    
    ShaderObject()
        : handle(0)
    {}
    void init(Type type) {
        handle = glCreateShaderObjectARB(
            type == Type::FRAGMENT_SHADER ? GL_FRAGMENT_SHADER_ARB :
            type == Type::VERTEX_SHADER   ? GL_VERTEX_SHADER_ARB   :
            GL_GEOMETRY_SHADER_ARB);
    }
    void set_source( const std::string& source) {
        int len = static_cast<int>(source.length());
        auto source_ptr = source.c_str();
        glShaderSourceARB(handle, 1, &source_ptr, &len);
    }
    bool set_source_from_file(const std::string& filename_) {
        std::ifstream f_in(filename_.c_str(), std::ios::binary);
        if (f_in.fail()) {
            std::cerr << "cannot open file: " << filename_ << std::endl;
            return false;
        }
        filename = filename_;
        
        std::ostringstream oss;
        oss << f_in.rdbuf();
        
        set_source(oss.str());
        
        f_in.close();
        return true;
    }
    void compile() {
        // compile
        glCompileShaderARB( handle );
        
        // check errors
        GLint	result;
        glGetObjectParameterivARB( handle, GL_OBJECT_COMPILE_STATUS_ARB, &result );
        if (result == GL_FALSE ) {
            int	length;
            glGetObjectParameterivARB( handle, GL_OBJECT_INFO_LOG_LENGTH_ARB, &length );
            if ( length > 0 ) {
                int	l;
                std::string info_log;
                info_log.resize(length);
                glGetInfoLogARB(handle, length, &l, &info_log[0]);
                std::cerr << "(filename:" << filename << ") " << info_log << std::endl;
            }
        }
    }
};

}

