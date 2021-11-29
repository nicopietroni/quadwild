############################ PROJECT FILES ############################

include(libs.pri)
include($$QUADRETOPOLOGY_PATH/quadretopology.pri)

#HEADERS =
#    parafashion.h

SOURCES = \
    quadrangulate_commandline.cpp \

############################ TARGET ############################

#App config
TARGET = quadrangulate_commandline
CONFIG += c++11
CONFIG -= app_bundle

macx {
    QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.13
    QMAKE_MAC_SDK = macosx10.13
}


DEFINES += GLEW_STATIC
DEFINES += INCLUDE_TEMPLATES
DEFINES += COMISO_FIELD


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

#vcglib
INCLUDEPATH += $$VCGLIBPATH

#eigen
INCLUDEPATH += $$EIGENPATH

#tracing
INCLUDEPATH += $$MESHFIELDPATH
INCLUDEPATH += $$MESHTRACEPATH
INCLUDEPATH += $$QUADRANGULATEPATH

#libigl
#INCLUDEPATH += $$LIBIGLPATH/include
HEADERS += \
#    $$LIBIGLPATH/include/igl/principal_curvature.h \
#    $$LIBIGLPATH/include/igl/copyleft/comiso/nrosy.h \
#    local_para_smooth.h

SOURCES += \
    $$LIBIGLPATH/include/igl/principal_curvature.cpp \
    $$LIBIGLPATH/include/igl/copyleft/comiso/nrosy.cpp \
    $$QUADRANGULATEPATH/load_save.cpp \
    $$QUADRANGULATEPATH/quad_from_patches.cpp

#AntTweakBar
#INCLUDEPATH += $$ANTTWEAKBARPATH/include
#LIBS += -L$$ANTTWEAKBARPATH/lib -lAntTweakBar

win32{ # Awful problem with windows..
    DEFINES += NOMINMAX
}

#glew
#LIBS += -lGLU
#INCLUDEPATH += $$GLEWPATH/include
#SOURCES += $$GLEWPATH/src/glew.c

#comiso
LIBS += -L$$COMISOPATH/build/Build/lib/CoMISo/ -lCoMISo
INCLUDEPATH += $$COMISOPATH/..

#gmm (we have to use comiso gmm)
INCLUDEPATH += $$GMMPATH/include

#libigl
INCLUDEPATH += $$LIBIGLPATH/include
QMAKE_CXXFLAGS += -isystem $$LIBIGL_PATH/include/

#Boost
INCLUDEPATH += $$BOOST_PATH

#Gurobi
INCLUDEPATH += $$GUROBI_PATH/include
LIBS += -L$$GUROBI_PATH/lib -lgurobi_g++4.2 -lgurobi90
DEFINES += GUROBI_DEFINED


SOURCES += $$VCGLIBPATH/wrap/ply/plylib.cpp


# Mac specific Config required to avoid to make application bundles
macx{
#    CONFIG -= app_bundle
#    LIBS += $$ANTTWEAKBARPATH/lib/libAntTweakBar.dylib
    QMAKE_POST_LINK +="cp -P ../../../code/lib/AntTweakBar1.16/lib/libAntTweakBar.dylib . ; "
    QMAKE_POST_LINK +="install_name_tool -change ../lib/libAntTweakBar.dylib ./libAntTweakBar.dylib $$TARGET ; "
    QMAKE_POST_LINK +="install_name_tool -change libCoMISo.dylib $$COMISOPATH/build/Build/lib/CoMISo/libCoMISo.dylib $$TARGET ;"
    DEPENDPATH += .
}



