#include "ParameterSet.h"
#include <iostream>
#include <chrono>
#include <sstream>


/// Function to create a default seed based on what time it is
size_t ParameterSet::makeDefaultSeed()
{
    return static_cast<size_t>(std::chrono::high_resolution_clock::now().
     time_since_epoch().count());
}






