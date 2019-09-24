# This script compiles the EGS simulation program into an executable
# Qt5's qmake build system is used to produce the makefiles from the project file
# Then make calls GNU's GCC to compile

module load Qt5
qmake ./EGS.pro
make --silent release
