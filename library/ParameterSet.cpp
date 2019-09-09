#include "ParameterSet.h"
#include <iostream>
#include <chrono>
#include <sstream>
#include <cassert>
#include <cstdint>

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

namespace fnv
{
  constexpr uint64_t _(uint64_t h, const char* s)
  {
    return (*s == 0) ? h :
      _((h * 1099511628211ull) ^ static_cast<uint64_t>(*s), s+1);
  }
}

constexpr uint64_t _(const char* s)
{
  return fnv::_(14695981039346656037ull, s);
}

uint64_t _(const std::string& s)
{
  return fnv::_(14695981039346656037ull, s.data());
}

void ParameterSet::readParams(std::ifstream file)
{

    // Here goes several rounds of input reading
    // Evaluate input
    // Update the corresponding parameter with the next value
    // Use while file >> string to keep going through the parameters

    std::string input;
    while (file >> input) {

        switch (_(input)) {

        case _("maxResourceCapacity"): file >> maxResourceCapacity; break;
        case _("maxResourceGrowth"): file >> maxResourceGrowth; break;
        case _("habitatSymmetry"): file >> habitatSymmetry; break;
        case _("ecoSelCoeff"): file >> ecoSelCoeff; break;
        case _("initialPopSize"): file >> initialPopSize; break;
        case _("dispersalRate"): file >> dispersalRate; break;
        case _("birthRate"): file >> birthRate; break;
        case _("survivalProb"): file >> survivalProb; break;
        case _("matePreferenceStrength"):
            file >> matePreferenceStrength; break;
        case _("mateEvalutationCost"): file >> mateEvaluationCost; break;
        case _("nEcoLoci"): file >> nEcoLoci; break;
        case _("nMatLoci"): file >> nMatLoci; break;
        case _("nNtrLoci"): file >> nNtrLoci; break;
        case _("nEcoEdges"): file >> nEcoEdges; break;
        case _("nMatEdges"): file >> nMatEdges; break;
        case _("nNtrEdges"): file >> nNtrEdges; break;
        case _("nChromosomes"): file >> nChromosomes; break;
        case _("mutationRate"): file >> mutationRate; break;
        case _("recombinationRate"): file >> recombinationRate; break;
        case _("genomeLength"): file >> genomeLength; break;
        case _("freqSNP"): file >> freqSNP; break;
        case _("isFemaleHeterogamy"): file >> isFemaleHeterogamy; break;
        case _("isGeneticArchitecture"): file >> isGenerateArchitecture; break;
        case _("architectureFileName"): file >> architectureFileName; break;
        case _("scaleA"):
            for (size_t i = 0u; i < 3u; ++i) file >> scaleA[i];
            break;
        case _("scaleD"):
            for (size_t i = 0u; i < 3u; ++i) file >> scaleD[i];
            break;
        case _("scaleI"):
            for (size_t i = 0u; i < 3u; ++i) file >> scaleI[i];
            break;
        case _("scaleE"):
            for (size_t i = 0u; i < 3u; ++i) file >> scaleE[i];
            break;
        case _("skewnesses"):
            for (size_t i = 0u; i < 3u; ++i) file >> skewnesses[i];
            break;
        case _("effectSizeShape"): file >> effectSizeShape; break;
        case _("effectSizeScale"): file >> effectSizeScale; break;
        case _("interactionWeightShape"):
            file >> interactionWeightShape; break;
        case _("interactionWeightScale"):
            file >> interactionWeightScale; break;
        case _("dominanceVariance"): file >> dominanceVariance; break;
        case _("tBurnIn"): file >> tBurnIn; break;
        case _("tEndSim"): file >> tEndSim; break;
        case _("tSave"): file >> tSave; break;
        case _("tiny"): file >> tiny; break;
        case _("seed"): file >> seed; break;
        case _("record"): file >> record; break;

        default:
            throw std::runtime_error("Invalid parameter name"); break;


        }

    }

}





