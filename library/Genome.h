#ifndef EXPLICITGENOMESPECIATION_GENOME_H
#define EXPLICITGENOMESPECIATION_GENOME_H

typedef std::vector<double> Vector;
typedef std::vector<size_t> Uector;

/// A container of constant genetic features across the genome
struct Genome
{

    Genome(const std::vector<size_t>&, const size_t&, const size_t&,
     const double&, const double&, const double&);

    size_t nloci;
    Vector chromosomes;

    Uector traits;
    Vector locations;
    Vector effects;
    Vector dominances;

    // A vector of locus epistatic interaction partners and weights
    // std::vector<std::vector<std::pair<size_t, double> > > interactions;

    // Member functions
    Vector makeChromosomes(const size_t&);
    Uector makeEncodedTraits(const Uector&);
    Vector makeLocations();
    Vector makeEffects(const double&, const double&);
    Vector makeDominances(const double& = 1.0);

    void setLocationsEffectSizesAndDominance(const size_t&,
     const double&, const double&);


};

#endif
