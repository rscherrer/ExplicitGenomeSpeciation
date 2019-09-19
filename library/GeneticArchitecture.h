#ifndef EXPLICITGENOMESPECIATION_GENETICARCHITECTURE_H
#define EXPLICITGENOMESPECIATION_GENETICARCHITECTURE_H

#include "ParameterSet.h"
#include "Random.h"
#include "Genome.h"
#include "Network.h"
#include "types.h"
#include <vector>
#include <cstddef>
#include <cassert>
#include <stddef.h>


typedef std::pair<size_t, size_t> Edge;
typedef std::vector<Network> MultiNet;


class GeneticArchitecture {

public:

    GeneticArchitecture(const ParameterSet &pars) :
        nChromosomes(pars.getNChromosomes()),
        nLoci(pars.getNLoci()),
        nLociPerTrait(pars.getNLociPerTrait()),
        nEdgesPerTrait(pars.getNEdgesPerTrait()),
        skewnesses(pars.getSkewnesses()),
        effectSizeShape(pars.getEffectSizeShape()),
        effectSizeScale(pars.getEffectSizeScale()),
        interactionWeightShape(pars.getInteractionWeightShape()),
        interactionWeightScale(pars.getInteractionWeightScale()),
        dominanceVariance(pars.getDominanceVariance()),
        snpFreq(pars.getSNPFreq()),
        femHeterogamy(pars.getIsFemaleHeterogamy()),
        recombinationRate(pars.getRecombinationRate()),
        scaleA(pars.getScaleA()),
        scaleD(pars.getScaleD()),
        scaleI(pars.getScaleI()),
        scaleE(pars.getScaleE()),
        genome(makeGenome()),
        networks(makeNetworks()),
        traits(makeEncodedTraits()),
        locations(makeLocations()),
        dominances(makeDominances())
    {
        assert(networks.size() == 3u);
    }

    Genome getGenome() const { return genome; }
    MultiNet getNetworks() const { return networks; }

    size_t nChromosomes;
    size_t nLoci;
    vecUns nLociPerTrait;
    vecUns nEdgesPerTrait;
    vecDbl skewnesses;
    double effectSizeShape;
    double effectSizeScale;
    double interactionWeightShape;
    double interactionWeightScale;
    double dominanceVariance;
    double snpFreq;
    bool femHeterogamy;
    double recombinationRate;

    vecDbl scaleA;
    vecDbl scaleD;
    vecDbl scaleI;
    vecDbl scaleE;

    Genome genome;
    MultiNet networks;

    vecUns traits;
    vecDbl locations;
    vecDbl dominances;


private:

    vecDbl makeChromosomes();
    Genome makeGenome();
    MultiNet makeNetworks();

    vecUns makeEncodedTraits();
    vecDbl makeLocations();
    vecDbl makeDominances();

};




#endif
