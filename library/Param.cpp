#include "Param.h"
#include <iostream>
#include <chrono>
#include <sstream>
#include <cassert>
#include <cstdint>

void Param::capEdges()
{
    for (size_t trait = 0u; trait < 3u; ++trait) {
        const size_t n = nLociPerTrait[trait];

        // Number of edges in a compete graph with N vertices
        const size_t emax = n * (n - 1u) / 2u;

        // Cap the number of edges
        if (nEdgesPerTrait[trait] > emax) nEdgesPerTrait[trait] = emax;
    }
}

Param::Param() : seed(makeDefaultSeed())
{

    capEdges();
    checkParams();

}

/// Function to create a default seed based on what time it is
size_t Param::makeDefaultSeed()
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

void Param::readParams(std::ifstream &file)
{

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
        case _("seed"): file >> seed; break;
        case _("record"): file >> record; break;

        default:
            throw std::runtime_error("Invalid parameter name: " + input); break;

        }
    }

    // Now update interactive parameters
    nLoci = nEcoLoci + nMatLoci + nNtrLoci;
    nLociPerTrait = { nEcoLoci, nMatLoci, nNtrLoci };
    nEdgesPerTrait = { nEcoEdges, nMatEdges, nNtrEdges };

    capEdges();

    checkParams();

    std::cout << "Parameters were read in succesfully.\n";

}

void Param::checkParams()
{
    std::string msg = "No error detected";

    if (dispersalRate < 0.0)
        msg = "Dispersal rate should be positive";
    if (dispersalRate > 1.0)
        msg = "Dispersal rate should be at most one";
    if (birthRate < 0.0)
        msg = "Birth rate should be positive";
    if (habitatSymmetry < 0.0)
        msg = "Habitat symmetry should be positive";
    if (habitatSymmetry > 1.0)
        msg = "Habitat symmetry should be at most one";
    if (survivalProb < 0.0)
        msg = "Survival probability should be positive";
    if (survivalProb > 1.0)
        msg = "Survival probability should be at most one";
    if (ecoSelCoeff < 0.0)
        msg = "Selection coefficient should be positive";
    if (matePreferenceStrength < 0.0)
        msg = "Mate preference strength should be positive";
    if (mateEvaluationCost < 0.0)
        msg = "Mate evaluation cost should be positive";
    if (maxResourceCapacity < 0.0)
        msg = "Maximum resource capacity should be positive";
    if (maxResourceGrowth < 0.0)
        msg = "Maximum resource growth should be positive";
    if (nEcoLoci <= 1u)
        msg = "Numer of ecological loci should be at least two";
    if (nMatLoci <= 1u)
        msg = "Number of mating loci should be at least two";
    if (nNtrLoci <= 1u)
        msg = "Number of neutral loci should be at least two";
    if (nChromosomes == 0u)
        msg = "Number of chromosomes should be at least one";
    if (nLoci <= 5u)
        msg = "Total number of loci should be at least six";
    if (freqSNP < 0.0)
        msg = "Frequency of SNPs should be positive";
    if (freqSNP > 1.0)
        msg = "Frequency of SNPs should be at most one";
    if (mutationRate < 0.0)
        msg = "Mutation rate should be positive";
    if (mutationRate > 1.0)
        msg = "Mutation rate should be at most one";
    if (genomeLength < 0.0)
        msg = "Genome length should be positive";
    if (recombinationRate < 0.0)
        msg = "Recombination rate should be positive";
    for (size_t i = 0u; i < 3u; ++i) {
        if (skewnesses[i] < 0.0)
            msg = "Skewness should be positive";
        if (scaleA[i] < 0.0)
            msg = "Additive scaling should be positive";
        if (scaleD[i] < 0.0)
            msg = "Dominance scaling should be positive";
        if (scaleI[i] < 0.0)
            msg = "Interaction scaling should be positive";
        if (scaleE[i] < 0.0)
            msg = "Environmental scaling should be positive";
    }
    if (effectSizeShape < 0.0)
        msg = "Effect size shape should be positive";
    if (effectSizeScale < 0.0)
        msg = "Effect size scale should be positive";
    if (interactionWeightShape < 0.0)
        msg = "Interaction weight shape should be positive";
    if (interactionWeightScale < 0.0)
        msg = "Interaction weight scale should be positive";
    if (dominanceVariance < 0.0)
        msg = "Dominance variance should be positive";

    if(msg != "No error detected")
        throw std::runtime_error(msg);

}

void Param::setMatePreferenceStrength(const double &alpha)
{
    matePreferenceStrength = alpha;
}

void Param::setNLociPerTrait(const vecUns &locipertrait)
{
    nLociPerTrait = locipertrait;
    nLoci = nLociPerTrait[0u] + nLociPerTrait[1u] + nLociPerTrait[2u];
    capEdges();
}

void Param::setNEdgesPerTrait(const vecUns &edgespertrait)
{
    nEdgesPerTrait = edgespertrait;
    capEdges();
}

void Param::setInteractionWeightScale(const double &x)
{
    interactionWeightScale = x;
}




