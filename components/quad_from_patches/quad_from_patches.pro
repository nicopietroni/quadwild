############################ TARGET AND FLAGS ############################

#App config
TARGET = quad_from_patches
TEMPLATE = app
CONFIG += c++11
CONFIG -= app_bundle

#Debug/release optimization flags
CONFIG(debug, debug|release){
    DEFINES += DEBUG
}
CONFIG(release, debug|release){
    DEFINES -= DEBUG
    #just uncomment next line if you want to ignore asserts and got a more optimized binary
    CONFIG += FINAL_RELEASE
}

#Final release optimization flag
FINAL_RELEASE {
    unix:!macx{
        QMAKE_CXXFLAGS_RELEASE -= -g -O2
        QMAKE_CXXFLAGS += -O3 -DNDEBUG
    }
}

macx {
    QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.13
    QMAKE_MAC_SDK = macosx10.13
}


############################ LIBRARIES ############################

#Setting library paths
include(../../libs/libs.pri)

#Configuration
include(configuration.pri)

#Quad retopology
include($$QUADRETOPOLOGY_PATH/quadretopology.pri)

#Libigl
INCLUDEPATH += $$LIBIGL_PATH/include/
QMAKE_CXXFLAGS += -isystem $$LIBIGL_PATH/include/

#Vcglib
INCLUDEPATH += $$VCGLIB_PATH

#Eigen
INCLUDEPATH += $$EIGEN_PATH

#Boost
INCLUDEPATH += $$BOOST_PATH

#Gurobi
INCLUDEPATH += $$GUROBI_PATH/include
LIBS += -L$$GUROBI_PATH/lib -l$$GUROBI_COMPILER -l$$GUROBI_LIB
DEFINES += GUROBI_DEFINED

#Parallel computation
unix:!mac {
    QMAKE_CXXFLAGS += -fopenmp
    LIBS += -fopenmp
}
#macx{
#    QMAKE_CXXFLAGS += -Xpreprocessor -fopenmp -lomp -I/usr/local/include
#    QMAKE_LFLAGS += -lomp
#    LIBS += -L /usr/local/lib/usr/local/lib/libomp.dylib
#}

win32{
    DEFINES += NOMINMAX # Awful problem with windows..
    DEFINES *= _USE_MATH_DEFINES
    DEFINES *= _SCL_SECURE_NO_DEPRECATE
    QMAKE_CXXFLAGS *= /bigobj
}

############################ PROJECT FILES ############################

SOURCES +=  \
    load_save.cpp \
    main.cpp \
    quad_from_patches.cpp \
    quad_mesh_tracer.cpp

HEADERS += \
    load_save.h \
    mesh_types.h \
    quad_from_patches.h \
    quad_mesh_tracer.h \
    smooth_mesh.h

#Vcg ply (needed to save ply files)
HEADERS += \
    $$VCGLIB_PATH/wrap/ply/plylib.h
SOURCES += \
    $$VCGLIB_PATH/wrap/ply/plylib.cpp

