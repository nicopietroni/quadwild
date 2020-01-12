LPSOLVER_PATH = $$PWD/patterns/ktmethod/lp_solve
LPSOLVER_PATH0 = $$PWD/patterns/ktmethod/lp_solve/bfp/
LPSOLVER_PATH1 = $$PWD/patterns/ktmethod/lp_solve/colamd/

win32:CONFIG(release, debug|release): LIBS += -L$$LPSOLVER_PATH/release/ -llpsolve55
else:win32:CONFIG(debug, debug|release): LIBS += -L$$LPSOLVER_PATH/debug/ -llpsolve55
else:unix: LIBS += -L$$LPSOLVER_PATH/ -llpsolve55

DEPENDPATH += $$LPSOLVER_PATH

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$LPSOLVER_PATH/release/liblpsolve55.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$LPSOLVER_PATH/debug/liblpsolve55.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$LPSOLVER_PATH/release/lpsolve55.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$LPSOLVER_PATH/debug/lpsolve55.lib
else:unix: PRE_TARGETDEPS += $$LPSOLVER_PATH/liblpsolve55.a

LIBS += -ldl -lm
LIBS += -L$$LPSOLVER_PATH -llpsolve55
INCLUDEPATH += $$LPSOLVER_PATH $$EIGEN_PATH

QMAKE_CXXFLAGS += -std=c++11 -fpermissive

SOURCES += \
    $$PWD/patterns/myutils.cpp \
    $$PWD/patterns/ktmethod/patchgen/get_boundary_geometry.cpp \
    $$PWD/patterns/ktmethod/patchgen/get_default_parameter.cpp \
    $$PWD/patterns/ktmethod/patchgen/get_param_str.cpp \
    $$PWD/patterns/ktmethod/patchgen/extradefinition.cpp

HEADERS += \
    $$LPSOLVER_PATH/lp_lib.h \
    $$PWD/patterns/laplacianreconstruction.h \
    $$PWD/patterns/meshtypes.h \
    $$PWD/patterns/myutils.h \
    $$PWD/patterns/patchg.h \
    $$PWD/patterns/ktmethod/patchgen/Permutation.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern_all.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern_6_3.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern_6_2.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern_6_1.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern_6_0.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern_5_4.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern_5_3.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern_5_2.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern_5_1.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern_5_0.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern_4_4.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern_4_3.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern_4_2.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern_4_1.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern_4_0.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern_3_3.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern_3_2.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern_3_1.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern_3_0.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern_2_1.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern_2_0.h \
    $$PWD/patterns/ktmethod/patchgen/Pattern.h \
    $$PWD/patterns/ktmethod/patchgen/PatchParam.h \
    $$PWD/patterns/ktmethod/patchgen/ILP.h \
    $$PWD/patterns/ktmethod/patchgen/generate_subtopology.h \
    $$PWD/patterns/ktmethod/patchgen/decl.h \
    $$PWD/patterns/ktmethod/patchgen/edgeloop.h
