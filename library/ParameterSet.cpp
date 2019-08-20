#include "ParameterSet.h"
#include <iostream>
#include <chrono>
#include <sstream>

void ParameterSet::capEdges()
{
    for (size_t trait = 0u; trait < 3u; ++trait) {
        const size_t n = nLociPerTrait[trait];

        // Number of edges in a compete graph with N vertices
        const size_t emax = n * (n - 1u) / 2u;

        // Cap the number of edges
        if (nEdgesPerTrait[trait] > emax) nEdgesPerTrait[trait] = emax;
    }
}

ParameterSet::ParameterSet()
{

    capEdges();

}

/// Function to create a default seed based on what time it is
size_t ParameterSet::makeDefaultSeed()
{
    return static_cast<size_t>(std::chrono::high_resolution_clock::now().
     time_since_epoch().count());
}






