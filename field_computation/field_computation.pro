############################ PROJECT FILES ############################

include(libs.pri)

#\macx: QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.13
#QMAKE_MAC_SDK = macosx10.13
#CONFIG += c++11

INCLUDEPATH += $$VCGLIBPATH
INCLUDEPATH += $$GLEWPATH/include
INCLUDEPATH += $$ANTPATH/include
INCLUDEPATH += $$EIGENPATH
INCLUDEPATH += $$LIBIGLPATH

HEADERS       = glwidget.h \
                triangle_mesh_type.h \
                $$LIBIGLPATH/igl/principal_curvature.h

SOURCES       = glwidget.cpp \
                main.cpp \
    		
QT           += opengl

DEFINES += GLEW_STATIC
DEFINES += INCLUDE_TEMPLATES
DEFINES += COMISO_FIELD

SOURCES += $$GLEWPATH/src/glew.c


INCLUDEPATH += $$COMISOPATH/gmm/include
INCLUDEPATH += $$COMISOPATH/Solver
INCLUDEPATH += $$COMISOPATH/..

SOURCES += $$VCGLIBPATH/wrap/ply/plylib.cpp
SOURCES += $$VCGLIBPATH/wrap/gui/trackball.cpp
SOURCES += $$VCGLIBPATH/wrap/gui/trackmode.cpp
SOURCES += $$VCGLIBPATH/wrap/qt/anttweakbarMapperNew.cpp
SOURCES += $$LIBIGLPATH/igl/principal_curvature.cpp
SOURCES += $$LIBIGLPATH/igl/copyleft/comiso/nrosy.cpp

############################ TARGET ############################

#App config
TARGET = field_computation

TEMPLATE        = app
CONFIG         += c++11
CONFIG         -= app_bundle

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

# Awful problem with windows..
win32{
  DEFINES += NOMINMAX
  LIBS +=$$ANTPATH/lib/AntTweakBar.lib
}

mac{
# Mac specific Config required to avoid to make application bundles
  CONFIG -= app_bundle
  LIBS +=$$ANTPATH/lib/libAntTweakBar.dylib
  QMAKE_POST_LINK +="cp -P ../../../code/lib/AntTweakBar1.16/lib/libAntTweakBar.dylib . ; "
  QMAKE_POST_LINK +="install_name_tool -change ../lib/libAntTweakBar.dylib ./libAntTweakBar.dylib $$TARGET ; "
  QMAKE_POST_LINK +="install_name_tool -change libCoMISo.dylib $$COMISOPATH/build/Build/lib/CoMISo/libCoMISo.dylib $$TARGET ;"
  LIBS += -L $$COMISOPATH/build/Build/lib/CoMISo/ -lCoMISo
  INCLUDEPATH += $$COMISOPATH/build/Build/lib/CoMISo
  DEPENDPATH += $$COMISOPATH/build/Build/lib/CoMISo
  DEPENDPATH += .
}

