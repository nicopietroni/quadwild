QT += core

DEFINES *= QT_NO_OPENGL_ES_2
DEFINES *= QT_NO_KEYWORDS

CONFIG += c++11

TARGET = polyMetrics

CONFIG -= app_bundle
CONFIG += console


CONFIG(release, debug|release): DEFINES += NDEBUG

win32 {
    QMAKE_CXXFLAGS_DEBUG *= -bigobj
}

LIBS_EXTERNAL = vcg eigen omp json
include(../libs.pri)

INCLUDEPATH *= ./external/matplotlib-cpp
DEPENDPATH  *= ./external/matplotlib-cpp
INCLUDEPATH *= ./external/nlohmann

INCLUDEPATH *= C:\Python27\include
LIBS *= -LC:\Python27\libs -lpython27

SOURCES *= \
    main.cpp

HEADERS += \
    mesh_def.h \
    matplotlibcpp.h





