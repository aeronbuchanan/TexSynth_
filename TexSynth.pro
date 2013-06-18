TEMPLATE = app
CONFIG += console
CONFIG -= qt
CONFIG += x11
CONFIG += thread

QMAKE_CXXFLAGS += -std=c++11

INCLUDEPATH += "../../CImg-1.5.4/"

SOURCES += main.cpp

HEADERS += \
    texSynth.h \
	texSynth.hpp \
    patch.h \
    seam.h \
    vecn.h \
    common.h \
    table.h \
    circularSeam.h

