#pragma once

#include <GL/glew.h>

namespace kt84 {

struct DisplayList {
    unsigned int list;
    bool is_valid;
    
    DisplayList()
        : list()
        , is_valid()
    {}
    void invalidate() { is_valid = false; }
    template <typename TRenderFunc>
    void render(TRenderFunc renderFunc) {
        if (!is_valid) {
            if (!list)
                list = glGenLists(1);
            
            glNewList(list, GL_COMPILE);
            renderFunc();
            glEndList();
            
            is_valid = true;
        }
        glCallList(list);
    }
};

}
