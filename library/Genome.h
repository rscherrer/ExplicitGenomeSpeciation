#ifndef EXPLICITGENOMESPECIATION_GENOME_H
#define EXPLICITGENOMESPECIATION_GENOME_H

/// A container of constant genetic features across the genome
struct Genome
{

    Genome(const std::vector<size_t>&, const size_t&, const size_t&,
     const double&, const double&);

    size_t nloci;
    std::vector<double> chromosomes;

    std::vector<size_t> traits;
    std::vector<double> locations;
    std::vector<double> effects;
    std::vector<double> dominances;

    // A vector of locus epistatic interaction partners and weights
    // std::vector<std::vector<std::pair<size_t, double> > > interactions;

    // Member functions
    std::vector<double> makeChromosomes(const size_t&);
    std::vector<size_t> makeEncodedTraits(const std::vector<size_t>&);
    void setLocationsEffectSizesAndDominance(const size_t&,
     const double&, const double&);

};

#endif
