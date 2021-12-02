LIBS_DIR = $$PWD/../../../../
message($$PWD)
message($$LIBS_DIR)
# VCGlib
contains(LIBS_EXTERNAL, vcg) {
	isEmpty(VCG_DIR) {
                VCG_DIR = $$LIBS_DIR/vcglib
	}
	INCLUDEPATH *= $$VCG_DIR
	DEPENDPATH  *= $$VCG_DIR

	contains(LIBS_EXTERNAL, glew) {
		# trackball
		SOURCES += $$VCG_DIR/wrap/gui/trackmode.cpp
		SOURCES += $$VCG_DIR/wrap/gui/trackball.cpp
	}
	# ply format
	SOURCES += $$VCG_DIR/wrap/ply/plylib.cpp
}

# VCGlibfork
contains(LIBS_EXTERNAL, vcgfork) {
        isEmpty(VCG_DIR) {
                VCG_DIR = C:\Users\thoma\Documents\GitHub\vcglib
        }
        INCLUDEPATH *= $$VCG_DIR
        DEPENDPATH  *= $$VCG_DIR

        contains(LIBS_EXTERNAL, glew) {
                # trackball
                SOURCES += $$VCG_DIR/wrap/gui/trackmode.cpp
                SOURCES += $$VCG_DIR/wrap/gui/trackball.cpp
        }
        # ply format
        SOURCES += $$VCG_DIR/wrap/ply/plylib.cpp
}

# libigl
contains(LIBS_EXTERNAL, igl) {
	# libigl dep
	LIBS_EXTERNAL *= eigen
	isEmpty(LIBIGL_DIR) {
		LIBIGL_DIR = $$LIBS_DIR/libigl
	}
	INCLUDEPATH *= $$LIBIGL_DIR/include
	DEPENDPATH  *= $$LIBIGL_DIR/include
}

#cinolib
contains(LIBS_EXTERNAL, cino) {
        # libigl dep
        LIBS_EXTERNAL *= eigen
        isEmpty(CINOLIB_DIR) {
                CINOLIB_DIR = $$LIBS_DIR/cinolib
        }
        INCLUDEPATH *= $$CINOLIB_DIR/include
        DEPENDPATH  *= $$CINOLIB_DIR/include
        message($$INCLUDEPATH)
}

# Eigen
contains(LIBS_EXTERNAL, eigen) {
	isEmpty(EIGEN_DIR) {
		EIGEN_DIR = $$LIBS_DIR/eigen
	}
	# To use standard Eigen includes
	INCLUDEPATH *= $$EIGEN_DIR
	DEPENDPATH  *= $$EIGEN_DIR
}

# GLEW
contains(LIBS_EXTERNAL, glew) {
	isEmpty(GLEW_DIR) {
                GLEW_DIR = $$LIBS_DIR/glew-2.0.0
	}

        #DEFINES     *= GLEW_STATIC
	INCLUDEPATH *= $$GLEW_DIR/include
	DEPENDPATH  *= $$GLEW_DIR/include
        win32 {
            LIBS *= $$GLEW_DIR\lib\Release\x64\glew32.lib
            LIBS *= -L$$GLEW_DIR\bin\Release\x64
            LIBS *= -lopengl32 -lGLU32
        }

        SOURCES *= $$GLEW_DIR/src/glew.c

	linux {
		LIBS    *= -lGLU
	}
}

#GLFW
contains (LIBS_EXTERNAL, glfw) {
    isEmpty(GLFW_DIR) {
        GLFW_DIR = $$LIBS_DIR/glfw
    }

    INCLUDEPATH *= $$GLFW_DIR/include
    DEPENDPATH  *= $$GLFW_DIR/include

    win32 {
        LIBS *= -L$$GLFW_DIR/lib
        LIBS *= -lglfw3
        LIBS *= -lopengl32
    }
}

#imgui
contains (LIBS_EXTERNAL, imgui) {
    isEmpty (IMGUI_DIR) {
        IMGUI_DIR = $$LIBS_DIR/imgui
    }

    INCLUDEPATH *= $$IMGUI_DIR
    DEPENDPATH  *= $$IMGUI_DIR

    SOURCES *= \
        $$IMGUI_DIR\imgui.cpp \
        $$IMGUI_DIR\imgui_draw.cpp \
        $$IMGUI_DIR\imgui_demo.cpp \
        $$IMGUI_DIR\imgui_widgets.cpp \
        $$IMGUI_DIR\imgui_impl_glfw.cpp \
        $$IMGUI_DIR\imgui_impl_opengl3.cpp

    HEADERS *= \
        $$IMGUI_DIR\imgui.h \
        $$IMGUI_DIR\imgui_internal.h \
        $$IMGUI_DIR\imgui_impl_glfw.h \
        $$IMGUI_DIR\imgui_impl_opengl3.h

    LIBS *= -lgdi32



    contains(LIBS_EXTERNAL, glew) {
        DEFINES *= IMGUI_IMPL_OPENGL_LOADER_GLEW
    }
}

#nanogui
contains (LIBS_EXTERNAL, nanogui) {
    LIBS_EXTERNAL *= eigen
    isEmpty(NANOGUI_DIR) {
        NANOGUI_DIR = $$LIBS_DIR/nanogui
    }

    DEFINES *= NANOGUI_GLAD NANOGUI_SHARED

    INCLUDEPATH *= $$NANOGUI_DIR/include
    DEPENDPATH  *= $$NANOGUI_DIR/include

    INCLUDEPATH *= $$NANOGUI_DIR/ext/nanovg/src
    DEPENDPATH  *= $$NANOGUI_DIR/ext/nanovg/src

    INCLUDEPATH *= $$NANOGUI_DIR/ext/glad/include
    DEPENDPATH  *= $$NANOGUI_DIR/ext/glad/include

#    SOURCES *= $$NANOGUI_DIR/ext/glad/src/glad.c

    INCLUDEPATH *= $$NANOGUI_DIR/ext/glfw/include
    DEPENDPATH  *= $$NANOGUI_DIR/ext/glfw/include


    win32 {
        CONFIG(release, debug|release) {
            message("NANOGUI release")
            LIBS *= $$NANOGUI_DIR/lib/Release/nanogui.lib
        }

        CONFIG(debug, debug|release) {
            message("NANOGUI debug")
            LIBS *= $$NANOGUI_DIR/lib/Debug/nanogui.lib
        }
    }

}
# anttweakbar
contains(LIBS_EXTERNAL, anttweakbar) {
	isEmpty(ANTTWEAKBAR_DIR) {
		ANTTWEAKBAR_DIR  = $$LIBS_DIR/AntTweakBar1.16
	}

	contains(LIBS_EXTERNAL , vcg) {
		SOURCES         *= $$VCG_DIR/wrap/qt/anttweakbarMapperNew.cpp
	}

	INCLUDEPATH     *= $$ANTTWEAKBAR_DIR/include
	
	# Awful problem with windows..
	win32 {
                #DEFINES += NOMINMAX
                LIBS    += -L$$ANTTWEAKBAR_DIR/lib -lAntTweakBar64
	}

	linux {
		LIBS    *= -L$$ANTTWEAKBAR_DIR/lib -lAntTweakBar
	}
	
	mac {
		# include specific anttweakbar version using different std library
		macx-g++49 {
			ANTTWEAKLIB = libAntTweakBar.dylib
		} else {
			contains(CONFIG, c++11) {
				ANTTWEAKLIB = libAntTweakBar_libc++.dylib
			} else {
				ANTTWEAKLIB = libAntTweakBar.dylib
			}
		}
		LIBS *= $$ANTTWEAKBAR_DIR/lib/$$ANTTWEAKLIB

#		equals(TEMPLATE, app) {
			QMAKE_POST_LINK += "cp -P $$ANTTWEAKBAR_DIR/lib/$$ANTTWEAKLIB $$DESTDIR/libAntTweakBar.dylib ;"
			QMAKE_POST_LINK += "install_name_tool -change ../lib/libAntTweakBar.dylib ./libAntTweakBar.dylib $$DESTDIR/$$TARGET ;"
#		}
	}
}

# Clipper
contains(LIBS_EXTERNAL, clipper) {
	isEmpty(CLIPPER_DIR) {
		CLIPPER_DIR = $$LIBS_DIR/clipper_ver6.4.2/cpp
	}
	INCLUDEPATH *= $$CLIPPER_DIR
	DEPENDPATH  *= $$CLIPPER_DIR
	SOURCES     *= $$CLIPPER_DIR/clipper.cpp
}

# Triangle
contains(LIBS_EXTERNAL, triangle) {
	isEmpty(TRIANGLE_DIR) {
		TRIANGLE_DIR = $$LIBS_DIR/triangle
	}
	INCLUDEPATH *= $$TRIANGLE_DIR
	DEPENDPATH  *= $$TRIANGLE_DIR
	SOURCES     *= $$TRIANGLE_DIR/triangle.c
}

# clang OpenMP
contains(LIBS_EXTERNAL, omp) {
    isEmpty(OMP_DIR) {
#		-I/opt/local/include/libomp -L/opt/local/lib/libomp -fopenmp
        OMP_DIR = /opt/local
    }
	win32 {
		QMAKE_CXXFLAGS *= -openmp
	}
	macx {
		INCLUDEPATH    *= $$OMP_DIR/include/libomp
		LIBS           *= -L$$OMP_DIR/lib/libomp

		QMAKE_CXXFLAGS *= -fopenmp
		QMAKE_LFLAGS   *= -fopenmp
	}

        unix:!macx {
                QMAKE_CXXFLAGS *= -fopenmp
                QMAKE_LFLAGS += -fopenmp -lgomp
        }
}

#
# ShapeOp
contains(LIBS_EXTERNAL, shapeop) {
	isEmpty(SHAPEOP_DIR) {
		SHAPEOP_DIR = $$LIBS_DIR/shapeop010/libShapeOp
	}
	DEFINES     += SHAPEOP_HEADER_ONLY=1
	INCLUDEPATH *= $$SHAPEOP_DIR/src
	DEPENDPATH  *= $$SHAPEOP_DIR/src

	contains(LIBS_EXTERNAL, omp) {
		DEFINES *= OMP_NESTED=TRUE
		DEFINES += SHAPEOP_OPENMP
		DEFINES *= EIGEN_DONT_PARALLELIZE
	}
}

# OpenVDB
contains(LIBS_EXTERNAL, openvdb) {
	isEmpty(OPENVDB_DIR) {
		OPENVDB_DIR = $$LIBS_DIR/openvdb
    }


    win32 {
        OPENVDB_DIR = $$LIBS_DIR/libs/include
        QMAKE_CXXFLAGS *= -DOPENVDB_3_ABI_COMPATIBLE -DOPENVDB_DLL -DOPENEXR_DLL -DWIN32 -D_WINDOWS -D_WIN320
        QMAKE_CXXFLAGS *= -DNOMINMAX
        INCLUDEPATH *= C:\local\boost_1_63_0
 Release:
{
    message("Openvdb in release mode..")
    LIBS += -LC:\local\boost_1_63_0\lib64-msvc-14.0 \
            -L$$LIBS_DIR/libs/lib -lHalf \
            -L$$LIBS_DIR/libs/lib/zlib -lzlib \
            -L$$LIBS_DIR/libs/lib/GLFW -lglfw3 \
            -L$$LIBS_DIR/libs/lib/blosc -lblosc \
            -L$$LIBS_DIR/libs/lib/glew -lglew32 \
            -L$$LIBS_DIR/libs/lib/tbb -ltbb \
            -L$$LIBS_DIR/libs/lib/tbb -ltbbmalloc \
            -L$$LIBS_DIR/libs/lib/tbb -ltbbproxy \
            -L$$LIBS_DIR/libs/lib -lIex \
            -L$$LIBS_DIR/libs/lib/ -lopenvdb \
}
    Debug: {
        message("Openvdb in debug mode..")
     LIBS += -LC:\local\boost_1_63_0\lib64-msvc-14.0 \
            -L$$LIBS_DIR/libs/lib -lHalfd \
            -L$$LIBS_DIR/libs/lib/zlib -lzlibd \
            -L$$LIBS_DIR/libs/lib/GLFW -lglfw3d \
            -L$$LIBS_DIR/libs/lib/blosc -lbloscd \
            -L$$LIBS_DIR/libs/lib/glew -lglew32 \
            -L$$LIBS_DIR/libs/lib/tbb -ltbb_debug \
            -L$$LIBS_DIR/libs/lib/tbb -ltbbmalloc_debug \
            -L$$LIBS_DIR/libs/lib/tbb -ltbbproxy_debug \
            -L$$LIBS_DIR/libs/lib -lIexd \
            -L$$LIBS_DIR/libs/lib/ -lopenvdbd \
}
    }
        INCLUDEPATH *= $$OPENVDB_DIR
        DEPENDPATH  *= $$OPENVDB_DIR
#        QMAKE_CXXFLAGS *= -DOPENVDB_3_ABI_COMPATIBLE

	macx {
        INCLUDEPATH *= $$OPENVDB_DIR
        DEPENDPATH  *= $$OPENVDB_DIR

		LIBS *= -L$$OPENVDB_DIR/build/openvdb
		LIBS *= /Users/mal/devel/openvdb/build/openvdb/libopenvdb.a

		INCLUDEPATH *= /opt/local/include

		LIBS *= -lz

		LIBS *= -L/opt/local/lib
		LIBS *= -lboost_iostreams-mt -lboost_system-mt -lboost_thread-mt

		LIBS *= -lblosc

		LIBS *= -ltbb -ltbbmalloc -ltbbmalloc_proxy
		LIBS *= -lHalf
		LIBS *= -lIex
		LIBS *= -lIlmThread
		LIBS *= -lImath
	}
        linux {
#                INCLUDEPATH *= $$OPENVDB_DIR
#                DEPENDPATH  *= $$OPENVDB_DIR

        LIBS *= -L/usr/local/lib
#        LIBS *= /usr/local/lib/libopenvdb.a

                INCLUDEPATH *= /usr/local/include
                DEPENDPATH  *= /usr/local/include

                LIBS *= -lz
            message("linking openvdb [linux]..")
#                LIBS *= -L/opt/local/lib
                LIBS *= -lboost_iostreams -lboost_system -lboost_thread
                LIBS *= -lblosc
                LIBS *= -ltbb -ltbbmalloc -ltbbmalloc_proxy
                LIBS *= -lHalf
                LIBS *= -lIex
                LIBS *= -lIlmThread
                LIBS *= -lImath
                LIBS *= -lopenvdb
        }
}

#GUROBI
contains (LIBS_EXTERNAL, gurobi) {
win32 {
        LIBS += -L$$LIBS_DIR/gurobi9/win64/lib

        CONFIG(debug, debug|release){
            message(debug)
            LIBS += -lgurobi_c++mdd2015 -lgurobi90
        } else {
            message(release)
            LIBS += -lgurobi_c++md2015 -lgurobi90
        }

    INCLUDEPATH += $$LIBS_DIR/gurobi9/win64/include
    DEPENDPATH += $$LIBS_DIR/gurobi9/win64/include
}
}
# CGAL
contains(LIBS_EXTERNAL, cgal) {

win32 {

    # gmp
    LIBS += C:\local\CGAL-4.11\auxiliary\gmp\lib\libgmp-10.lib
    # mpfr
    LIBS += C:\local\CGAL-4.11\auxiliary\gmp\lib\libmpfr-4.lib
    # CGAL
    INCLUDEPATH += C:/local/CGAL-4.11/include
    INCLUDEPATH += C:\local\CGAL-4.11\auxiliary\gmp\include

    Release: {
        LIBS +=  \
                -LC:\local\CGAL-4.11\build\lib\ -lCGAL_Core-vc140-mt-4.11 \
                -LC:\local\CGAL-4.11\build\lib\ -lCGAL-vc140-mt-4.11 \

    }
    Debug: {
        LIBS += \
                C:\local\CGAL-4.11\build\lib\CGAL_Core-vc140-mt-gd-4.11.lib \
                C:\local\CGAL-4.11\build\lib\CGAL-vc140-mt-gd-4.11.lib \
    }

    Release: LIBS += -LC:/local/boost_1_63_0/lib64-msvc-14.0/ -lboost_thread-vc140-mt-1_63
    Debug: LIBS += -LC:/local/boost_1_63_0/lib64-msvc-14.0/ -lboost_thread-vc140-mt-gd-1_63
    INCLUDEPATH += C:\local\boost_1_63_0
    INCLUDEPATH += C:/local/boost_1_63_0/lib64-msvc-14.0
    DEPENDPATH += C:/local/boost_1_63_0/lib64-msvc-14.0
}

	macx {
		isEmpty(CGAL_DIR) {
			CGAL_DIR = /opt/local
		}

#		LIBS_EXTERNAL *= boost
		INCLUDEPATH *= $$CGAL_DIR/include
		LIBS *= -L$$CGAL_DIR/lib

		LIBS *= -L/opt/local/lib
		LIBS *= -lboost_iostreams-mt -lboost_system-mt -lboost_thread-mt
		LIBS *= -lCGAL -lCGAL_Core
		LIBS *= -lmpfr
		LIBS *= -lgmp
		LIBS  *= -frounding-math
	}
	linux {
		isEmpty(CGAL_DIR) {
			CGAL_DIR = /opt/local
		}

		INCLUDEPATH *= $$CGAL_DIR/include
		LIBS *= -L$$CGAL_DIR/lib

		LIBS *= -L/opt/local/lib
		LIBS *= -lboost_iostreams -lboost_system -lboost_thread
		LIBS *= -lCGAL -lCGAL_Core
		LIBS *= -lmpfr
		LIBS *= -lgmp
		LIBS *= -frounding-math
	}
}

