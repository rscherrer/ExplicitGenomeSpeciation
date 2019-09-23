#include "library/Param.h"
#include "library/GenArch.h"
#include "library/Deme.h"
#include "library/Random.h"
#include "library/doMain.h"
#include <cassert>
#include <iostream>
#include <chrono>
#include <vector>
#include <string>


int main(int argc, char * argv[])
{

    // Convert arguments into a vector of strings
    const vecStrings args(argv, argv + argc);

    // Run the program
    return doMain(args);

}
