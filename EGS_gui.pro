# GUI main and other files
SOURCES += \
        gui/main.cpp \
        gui/maindialog.cpp \
        gui/qcustomplot.cpp

HEADERS += \
        gui/maindialog.h \
        gui/qcustomplot.h

FORMS += \
        gui/maindialog.ui


# Core files
include(EGS.pri)

# C++14
CONFIG += c++14
QMAKE_CXXFLAGS += -std=c++14

# High warnings levels
QMAKE_CXXFLAGS += -Wall -Wextra -Wnon-virtual-dtor -pedantic -Werror
QMAKE_CXXFLAGS_RELEASE += -O3

# Allow debug and release mode
CONFIG += debug_and_release

# In release mode, remove asserts
CONFIG(release, debug|release) {

  DEFINES += NDEBUG

}

# Qt
QT       += core gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

# Boost.Test
LIBS += -lboost_unit_test_framework
LIBS += -L"/usr/local/Cellar/boost/1.70.0/lib"
INCLUDEPATH += "/usr/local/Cellar/boost/1.70.0/include"
