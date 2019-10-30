#ifndef EXPLICITGENOMESPECIATION_TESTUTILITIES_H
#define EXPLICITGENOMESPECIATION_TESTUTILITIES_H

#include "library/Types.h"
#include <fstream>
#include <vector>
#include <iostream>

namespace tst
{

    void makeValidParamFile();
    void makeValidParamFile2();
    void makeInvalidParamName();
    void makeInvalidParamValue();
    void makeInvalidParamValue2();
    void makeParamFileWithArchitecture();
    void makeParamFileWithMissingArchitecture();
    vecDbl readfile(const std::string&);

}

#endif
