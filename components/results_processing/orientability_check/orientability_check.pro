QT += core

DEFINES *= QT_NO_OPENGL_ES_2
DEFINES *= QT_NO_KEYWORDS

CONFIG += c++11

TARGET = orientability_check

CONFIG -= app_bundle
CONFIG += console


CONFIG(release, debug|release): DEFINES += NDEBUG

LIBS_EXTERNAL = vcg eigen omp json
include(../libs.pri)

SOURCES *= \
    main.cpp

HEADERS += \
    mesh_def.h





