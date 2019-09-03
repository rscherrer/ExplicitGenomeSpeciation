#ifndef EXPLICITGENOMESPECIATION_GENETICARCHITECTURE_H
#define EXPLICITGENOMESPECIATION_GENETICARCHITECTURE_H

#include "ParameterSet.h"
#include "Random.h"
#include "Genome.h"
#include "Network.h"
#include <vector>
#include <cstddef>
#include <stddef.h>


typedef std::pair<size_t, size_t> Edge;
typedef std::vector<Network> MultiNet;
typedef std::vector<size_t> uVector;
typedef std::vector<double> dVector;

/// A container for all constant genetic features
class GeneticArchitecture {

public:

    GeneticArchitecture(const ParameterSet&);

    /// Getters
    Genome getGenome() const { return genome; }
    MultiNet getNetworks() const { return networks; }

private:

    size_t nChromosomes;
    size_t nLoci;
    uVector nLociPerTrait;
    uVector nEdgesPerTrait;
    dVector skewnesses;
    double effectSizeShape;
    double effectSizeScale;
    double interactionWeightShape;
    double interactionWeightScale;
    double dominanceVariance;

    Genome genome;
    MultiNet networks;

    // A set of vectors of loci underlying each trait
    //std::vector<std::vector<size_t> > traitUnderlyingLoci;

    /// Makers
    dVector makeChromosomes();
    Genome makeGenome();
    MultiNet makeNetworks();

};



#endif
