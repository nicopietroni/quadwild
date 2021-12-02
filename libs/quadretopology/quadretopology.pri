INCLUDEPATH += $$PWD

#Include patterns
include($$PWD/patterns/patterns.pri)

SOURCES += \
        $$PWD/quadretopology/includes/qr_charts.cpp \
        $$PWD/quadretopology/includes/qr_convert.cpp \
        $$PWD/quadretopology/includes/qr_ilp.cpp \
        $$PWD/quadretopology/includes/qr_patterns.cpp \
        $$PWD/quadretopology/includes/qr_mapping.cpp \
        $$PWD/quadretopology/includes/qr_utils.cpp \
        $$PWD/quadretopology/quadretopology.cpp

HEADERS += \
        $$PWD/quadretopology/includes/qr_convert.h \
        $$PWD/quadretopology/includes/qr_charts.h \
        $$PWD/quadretopology/includes/qr_convert.h \
        $$PWD/quadretopology/includes/qr_ilp.h \
        $$PWD/quadretopology/includes/qr_parameters.h \
        $$PWD/quadretopology/includes/qr_patterns.h \
        $$PWD/quadretopology/includes/qr_mapping.h \
        $$PWD/quadretopology/includes/qr_utils.h \
        $$PWD/quadretopology/quadretopology.h
