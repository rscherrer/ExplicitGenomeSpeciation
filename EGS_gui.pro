# GUI main and other files
SOURCES += \
        gui/main.cpp \
        gui/maindialog.cpp

HEADERS += \
        gui/maindialog.h

FORMS += \
        gui/maindialog.ui


# Core files
include(EGS.pri)

# C++14
CONFIG += c++17
QMAKE_CXXFLAGS += -std=c++17

# High warnings levels
QMAKE_CXXFLAGS += -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -Werror

# Allow debug and release mode
CONFIG += debug_and_release

# In release mode, turn on profiling
CONFIG(release, debug|release) {

  DEFINES += NDEBUG

  # gprof
  QMAKE_CXXFLAGS += -pg
  QMAKE_LFLAGS += -pg
}


# Qt
QT       += core gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets


LIBS += -L"/usr/local/Cellar/boost/1.70.0/lib"
INCLUDEPATH += "/usr/local/Cellar/boost/1.70.0/include"