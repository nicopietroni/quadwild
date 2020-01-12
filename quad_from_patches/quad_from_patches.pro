############################ PROJECT FILES ############################

include(libs.pri)

SOURCES +=  \
    includes/charts.cpp \
    includes/convert.cpp \
    includes/ilp.cpp \
    includes/mapping.cpp \
    includes/patterns.cpp \
    includes/utils.cpp \
    load_save.cpp \
    main.cpp \
    quad_from_patches.cpp

HEADERS += \
    includes/charts.h \
    includes/common.h \
    includes/convert.h \
    includes/ilp.h \
    includes/mapping.h \
    includes/orient_faces.h \
    includes/patterns.h \
    includes/utils.h \
    load_save.h \
    mesh_types.h \
    quad_from_patches.h

#DEFINES += SAVEMESHESFORDEBUG


############################ TARGET ############################

#App config
TARGET = quad_from_patches

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


############################ INCLUDES ############################

#Patterns
include($$PATTERNSPATH/patterns.pri)
INCLUDEPATH += $$PATTERNSPATH

#libigl
INCLUDEPATH += $$LIBIGLPATH/include/
QMAKE_CXXFLAGS += -isystem $$LIBIGLPATH/include/

#vcglib
INCLUDEPATH += $$VCGLIBPATH
DEFINES += NOCOMISO

#eigen
INCLUDEPATH += $$EIGENPATH

#gurobi
INCLUDEPATH += $$GUROBIPATH/include
LIBS += -L$$GUROBIPATH/lib -lgurobi_g++5.2 -lgurobi81
DEFINES += GUROBI_DEFINED

#lpsolve
LIBS += -llpsolve55
LIBS += -lm -ldl
INCLUDEPATH += $$LPSOLVEPATH
HEADERS += \
    $$LPSOLVEPATH/lp_lib.h \
