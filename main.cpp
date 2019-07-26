//#include "OutputFile.h"
#include "library/ParameterSet.h"
#include "library/GeneticArchitecture.h"
#include "library/Population.h"
#include "library/Random.h"
#include "library/doMain.h"
#include <cassert>
#include <iostream>
#include <chrono>
#include <vector>
#include <string>


/// Program to run an individual-based simulation of a speciation event with
/// explicit genomic features
int main(int argc, char * argv[])
{

    // Convert arguments into a vector of strings
    const std::vector<std::string> args(argv, argv + argc);

    // Run the program
    return doMain(args);

}
