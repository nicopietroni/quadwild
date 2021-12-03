############################ PROJECT FILES ############################

include(../../libs/libs.pri)

HEADERS = \
    glwidget.h \
    triangle_mesh_type.h

SOURCES = \
    glwidget.cpp \
    main.cpp

HEADERS += \
    fields/field_smoother.h \
    fields/n_polyvector.h \
    fields/polyroots.h
SOURCES += \
    fields/n_polyvector.cpp \
    fields/polyroots.cpp

DEFINES += GLEW_STATIC
DEFINES += INCLUDE_TEMPLATES
DEFINES += COMISO_FIELD


############################ TARGET ############################

#App config
TARGET = field_computation

TEMPLATE = app
CONFIG += qt
CONFIG += c++11
CONFIG -= app_bundle

QT += core gui opengl xml widgets

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
INCLUDEPATH += $$VCGLIB_PATH
SOURCES += $$VCGLIB_PATH/wrap/ply/plylib.cpp
SOURCES += $$VCGLIB_PATH/wrap/gui/trackball.cpp
SOURCES += $$VCGLIB_PATH/wrap/gui/trackmode.cpp
SOURCES += $$VCGLIB_PATH/wrap/qt/anttweakbarMapperNew.cpp

#eigen
INCLUDEPATH += $$EIGEN_PATH

#libigl
INCLUDEPATH += $$LIBIGL_PATH/include
HEADERS += \
    $$LIBIGL_PATH/include/igl/principal_curvature.h
SOURCES += \
    $$LIBIGL_PATH/include/igl/principal_curvature.cpp

#AntTweakBar
INCLUDEPATH += $$ANTTWEAKBAR_PATH/include
LIBS += -L$$ANTTWEAKBAR_PATH/lib -lAntTweakBar
win32{ # Awful problem with windows..
    DEFINES += NOMINMAX
}

#glew
LIBS += -lGLU
INCLUDEPATH += $$GLEW_PATH/include
SOURCES += $$GLEW_PATH/src/glew.c

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

# Mac specific Config required to avoid to make application bundles
macx{
#    CONFIG -= app_bundle
#    LIBS += $$ANTTWEAKBAR_PATH/lib/libAntTweakBar.dylib
    QMAKE_POST_LINK +="cp -P ../../../code/lib/AntTweakBar1.16/lib/libAntTweakBar.dylib . ; "
    QMAKE_POST_LINK +="install_name_tool -change ../lib/libAntTweakBar.dylib ./libAntTweakBar.dylib $$TARGET ; "
#    QMAKE_POST_LINK +="install_name_tool -change libCoMISo.dylib $$COMISO_PATH/build/Build/lib/CoMISo/libCoMISo.dylib $$TARGET ;"
    DEPENDPATH += .
#    DEPENDPATH += $$COMISO_PATH/build/Build/lib/CoMISo
#    INCLUDEPATH += $$COMISO_PATH/build/Build/lib/CoMISo
}

#Parallel computation
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
