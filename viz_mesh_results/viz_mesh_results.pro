############################ TARGET AND FLAGS ############################

#App config
TARGET = viz_mesh_result
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

#Setting library paths and configuration
include(configuration.pri)

#Libigl
#INCLUDEPATH += $$LIBIGL_PATH/include/
#QMAKE_CXXFLAGS += -isystem $$LIBIGL_PATH/include/

#Vcglib
INCLUDEPATH += $$VCGLIB_PATH

#Eigen
INCLUDEPATH += $$EIGEN_PATH



############################ PROJECT FILES ############################

SOURCES +=  \
    main.cpp \

HEADERS += \
    mesh_types.h \

#Vcg ply (needed to save ply files)
HEADERS += \
    $$VCGLIB_PATH/wrap/ply/plylib.h
SOURCES += \
    $$VCGLIB_PATH/wrap/ply/plylib.cpp
