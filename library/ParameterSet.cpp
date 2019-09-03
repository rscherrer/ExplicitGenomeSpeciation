#include "ParameterSet.h"
#include <iostream>
#include <chrono>
#include <sstream>
#include <cassert>

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

ParameterSet::ParameterSet() : seed(makeDefaultSeed())
{

    capEdges();

    assert(dispersalRate >= 0.0);
    assert(birthRate >= 0.0);
    assert(habitatSymmetry >= 0.0);
    assert(habitatSymmetry <= 1.0);
    assert(survivalProb >= 0.0);
    assert(survivalProb <= 1.0);
    assert(ecoSelCoeff >= 0.0);
    assert(matePreferenceStrength >= 0.0);
    assert(mateEvaluationCost >= 0.0);
    assert(maxResourceCapacity >= 0.0);
    assert(maxResourceGrowth >= 0.0);
    assert(nEcoLoci > 1u);
    assert(nMatLoci > 1u);
    assert(nNtrLoci > 1u);
    assert(nChromosomes > 0u);
    assert(nLoci > 5u);
    assert(freqSNP >= 0.0);
    assert(freqSNP <= 1.0);
    assert(mutationRate >= 0.0);
    assert(genomeLength >= 0.0);
    assert(recombinationRate >= 0.0);
    for (size_t i = 0u; i < 3u; ++i) {
        assert(skewnesses[i] >= 0.0);
        assert(scaleA[i] >= 0.0);
        assert(scaleD[i] >= 0.0);
        assert(scaleI[i] >= 0.0);
        assert(scaleE[i] >= 0.0);
    }
    assert(effectSizeShape >= 0.0);
    assert(effectSizeScale >= 0.0);
    assert(interactionWeightShape >= 0.0);
    assert(interactionWeightScale >= 0.0);
    assert(dominanceVariance >= 0.0);

}

/// Function to create a default seed based on what time it is
size_t ParameterSet::makeDefaultSeed()
{
    return static_cast<size_t>(std::chrono::high_resolution_clock::now().
     time_since_epoch().count());
}






