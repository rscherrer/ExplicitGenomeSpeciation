#ifndef EXPLICITGENOMESPECIATION_GENOME_H
#define EXPLICITGENOMESPECIATION_GENOME_H

/// A container of constant genetic features across the genome
struct Genome
{

    Genome(const size_t&, const std::vector<size_t>&, const size_t&,
     const double&, const double&);

    std::vector<size_t> traits;
    std::vector<double> locations;
    std::vector<double> effects;
    std::vector<double> dominances;

    // A vector of locus epistatic interaction partners and weights
    // std::vector<std::vector<std::pair<size_t, double> > > interactions;

    // Member functions
    std::vector<size_t> makeEncodedTraits(const size_t&,
     const std::vector<size_t>&) const noexcept;
    void setLocationsEffectSizesAndDominance(const size_t&, const size_t&,
     const double&, const double&);

};

#endif
