############################ PROJECT FILES ############################

include(../libs/libs.pri)
include($$QUADRETOPOLOGY_PATH/quadretopology.pri)

SOURCES += \
    functions.cpp \
    quadwild.cpp

HEADERS += \
    functions.h

############################ TARGET ############################

#App config
TARGET = quadwild
CONFIG += c++11
CONFIG -= app_bundle

macx {
    QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.13
    QMAKE_MAC_SDK = macosx10.13
}


DEFINES += GLEW_STATIC
DEFINES += INCLUDE_TEMPLATES


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

############################ INCLUDES ############################

MESHFIELD_PATH = ../components/field_computation
MESHTRACE_PATH =../components/field_tracing
QUADRANGULATE_PATH = ../components/quad_from_patches

HEADERS += \
    $$LIBIGL_PATH/include/igl/principal_curvature.h \
    $$QUADRANGULATE_PATH/load_save.h \
    $$QUADRANGULATE_PATH/quad_from_patches.h \
    $$VCGLIB_PATH/wrap/ply/plylib.h

SOURCES += \
    $$LIBIGL_PATH/include/igl/principal_curvature.cpp \
    $$QUADRANGULATE_PATH/load_save.cpp \
    $$QUADRANGULATE_PATH/quad_from_patches.cpp \
    $$VCGLIB_PATH/wrap/ply/plylib.cpp

#vcglib
INCLUDEPATH += $$VCGLIB_PATH

#eigen
INCLUDEPATH += $$EIGEN_PATH

#tracing
INCLUDEPATH += $$MESHFIELD_PATH
INCLUDEPATH += $$MESHTRACE_PATH
INCLUDEPATH += $$QUADRANGULATE_PATH

contains(DEFINES, COMISO_FIELD) {
    #comiso
    LIBS += -L$$COMISO_PATH/build/Build/lib/CoMISo/ -lCoMISo
    INCLUDEPATH += $$COMISO_PATH/..

    #gmm (we have to use comiso gmm)
    INCLUDEPATH += $$GMM_PATH/include

    HEADERS += \
        $$LIBIGL_PATH/include/igl/copyleft/comiso/nrosy.h
    SOURCES += \
        $$LIBIGL_PATH/include/igl/copyleft/comiso/nrosy.cpp
}

#libigl
INCLUDEPATH += $$LIBIGL_PATH/include
QMAKE_CXXFLAGS += -isystem $$LIBIGL_PATH/include/

#Boost
INCLUDEPATH += $$BOOST_PATH

#Gurobi
INCLUDEPATH += $$GUROBI_PATH/include
LIBS += -L$$GUROBI_PATH/lib -l$$GUROBI_COMPILER -l$$GUROBI_LIB
DEFINES += GUROBI_DEFINED

#Parallel computation (just in release)
unix:!mac {
    QMAKE_CXXFLAGS += -fopenmp
    LIBS += -fopenmp
}
#macx{
#    QMAKE_CXXFLAGS += -Xpreprocessor -fopenmp -lomp -I/usr/local/include
#    QMAKE_LFLAGS += -lomp
#    LIBS += -L /usr/local/lib /usr/local/lib/libomp.dylib
#}

win32{
    DEFINES += NOMINMAX # Awful problem with windows..
    DEFINES *= _USE_MATH_DEFINES
    DEFINES *= _SCL_SECURE_NO_DEPRECATE
    QMAKE_CXXFLAGS *= /bigobj
}

