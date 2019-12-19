#-------------------------------------------------
#
# Project created by QtCreator 2013-03-17T15:36:47
#
#-------------------------------------------------
\macx: QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.13
QMAKE_MAC_SDK = macosx10.13
CONFIG += c++11
QT       += core

QT       -= gui

TARGET = test
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

VCGLIBDIR = ../../../vcglib
EIGENLIB  = ../lib/libigl/external/nanogui/ext/eigen

#VCGLIBDIR = ../../../vcglib
#GLEWDIR   = ../../../code/lib/glew
#ANTDIR    = ../../../code/lib/AntTweakBar1.14

INCLUDEPATH += $$VCGLIBDIR
INCLUDEPATH += $$EIGENLIB
#INCLUDEPATH += ../lib/libigl/include
#INCLUDEPATH += ../lib/libigl/external/nanogui/ext/eigen/
#INCLUDEPATH += libigl-master/include/
#INCLUDEPATH += libigl-master/external/eigen/
#INCLUDEPATH += libigl-master/external/eigen/Eigen
#INCLUDEPATH += libigl-master/external/eigen/Eigen

#INCLUDEPATH += /opt/local/include/
SOURCES += quad_from_pathes.cpp
SOURCES += $$VCGLIBDIR/wrap/ply/plylib.cpp


#LIBS += -L/opt/local/lib -lcgal -lgmp -lmpfr -lboost_system-mt -lboost_thread-mt

DEPENDPATH += .
