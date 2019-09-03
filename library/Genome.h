#ifndef EXPLICITGENOMESPECIATION_GENOME_H
#define EXPLICITGENOMESPECIATION_GENOME_H

#include <vector>
#include <stddef.h>

typedef std::vector<double> dVector;
typedef std::vector<size_t> uVector;

/// A container of constant genetic features across the genome
struct Genome
{

    Genome(const uVector&, const size_t&, const size_t&,
     const double&, const double&, const double&);

    size_t nloci;
    dVector chromosomes;

    uVector traits;
    dVector locations;
    dVector effects;
    dVector dominances;

    // A vector of locus epistatic interaction partners and weights
    // std::vector<std::vector<std::pair<size_t, double> > > interactions;

    // Member functions
    dVector makeChromosomes(const size_t&);
    uVector makeEncodedTraits(const uVector&);
    dVector makeLocations();
    dVector makeEffects(const double&, const double&);
    dVector makeDominances(const double& = 1.0);

    void setLocationsEffectSizesAndDominance(const size_t&,
     const double&, const double&);


};

#endif
