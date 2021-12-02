#pragma once
#if defined(__APPLE__) && defined(__MACH__)
#   include <GLUT/glut.h>
#else
#   include <GL/glut.h>
#endif

namespace kt84 {
    namespace glut_util {
        namespace defaultcb {
            inline void display_pre() {
                glPushAttrib(GL_ALL_ATTRIB_BITS);
                glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
            }
            inline void display_post() {
                glPopAttrib();
#ifdef TW_INCLUDED
                TwDraw();
#endif
                glutSwapBuffers();
                glutReportErrors();
            }
            inline void reshape(int width, int height) {
                glViewport(0, 0, width, height);
#ifdef TW_INCLUDED
                TwWindowSize(width, height);
#endif
            }
            inline void keyboard(unsigned char key, int x, int y) {
#ifdef TW_INCLUDED
                if (TwEventKeyboardGLUT(key, x, y)) return glutPostRedisplay();
#endif
            }
            inline void special(int key, int x, int y) {
#ifdef TW_INCLUDED
                if (TwEventSpecialGLUT(key, x, y)) return glutPostRedisplay();
#endif
            }
            inline void mouse(int glut_button, int state, int x, int y) {
#ifdef TW_INCLUDED
                if (TwEventMouseButtonGLUT(glut_button, state, x, y)) return glutPostRedisplay();
#endif
            }
            inline void motion(int x, int y) {
#ifdef TW_INCLUDED
                if (TwEventMouseMotionGLUT(x, y)) return glutPostRedisplay();
#endif
            }
        }
        inline void init(
            int argc, char* argv[],
            unsigned int display_mode = GLUT_DOUBLE | GLUT_RGBA,
            double screen_window_size_ratio = 0.5,
            bool is_window_centered = true,
            const char* window_title = "glut_window",
            void (*display_pre)()                                   = defaultcb::display_pre,
            void (*display_main)()                                  = [] () {},
            void (*display_post)()                                  = defaultcb::display_post,
            void (*reshape)(int width, int height)                  = defaultcb::reshape,
            void (*keyboard)(unsigned char key, int x, int y)       = defaultcb::keyboard,
            void (*keyboardup)(unsigned char key, int x, int y)     = nullptr,
            void (*special)(int key, int x, int y)                  = defaultcb::special,
            void (*mouse)(int glut_button, int state, int x, int y) = defaultcb::mouse,
            void (*motion)(int x, int y)                            = defaultcb::motion, 
            void (*passive_motion)(int x, int y)                    = defaultcb::motion,
            void (*idle)()                                          = nullptr)
        {
            glutInit(&argc, argv);
            
            // Create window
            glutInitDisplayMode (display_mode);
            int screen_width  = glutGet(GLUT_SCREEN_WIDTH );
            int screen_height = glutGet(GLUT_SCREEN_HEIGHT);
            int window_width  = static_cast<int>(screen_width  * screen_window_size_ratio);
            int window_height = static_cast<int>(screen_height * screen_window_size_ratio);
            glutInitWindowSize(window_width, window_height);
            if (is_window_centered) {
                int screen_center_x = screen_width  / 2;
                int screen_center_y = screen_height / 2;
                glutInitWindowPosition(screen_center_x - window_width / 2, screen_center_y - window_height / 2);
            }
            glutCreateWindow(window_title);
            
            // Register callback functions
            static auto static_display_pre  = display_pre;
            static auto static_display_main = display_main;
            static auto static_display_post = display_post;
            glutDisplayFunc([] () {
                static_display_pre();
                static_display_main();
                static_display_post();
            });
            glutReshapeFunc(reshape);
            glutKeyboardFunc(keyboard);
            glutKeyboardUpFunc(keyboardup);
            glutSpecialFunc(special);
            glutMouseFunc(mouse);
            glutMotionFunc(motion);
            glutPassiveMotionFunc(passive_motion);
            glutIdleFunc(idle);
        }
    }
}
