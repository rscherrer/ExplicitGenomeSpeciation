#ifndef EXPLICITGENOMESPECIATION_GENETICARCHITECTURE_H
#define EXPLICITGENOMESPECIATION_GENETICARCHITECTURE_H

#include "ParameterSet.h"
#include <vector>
#include <list>


/// An edge is a pair of interacting loci
typedef std::pair<size_t, size_t> Edge;


/// A container of constant genetic features across the genome
struct Genome
{

    Genome(const size_t&, const std::vector<size_t>&, const size_t&, const double&, const double&);

    // A vector of locus encoded traits
    std::vector<size_t> encodedTraits;

    // A vector of locus locations
    std::vector<double> locations;

    // A vector of locus effect sizes
    std::vector<double> effectSizes;

    // A vector of locus dominance coefficients
    std::vector<double> dominanceCoeffs;

    // A vector of locus epistatic interaction partners and weights
    std::vector<std::vector<std::pair<size_t, double> > > interactions;

    // All the above vectors should have the same length, use assert for that

};


/// A container for a gene regulatory network
struct Network
{

    Network(const size_t&, const size_t&, const double&);

    // A map of interacting genes
    std::vector<Edge> map;

    // Number of vertices, edges and skewness
    size_t nVertices;
    size_t nEdges;
    double skewness;

    // Functions to make the network
    std::vector<Edge> makeNetwork(const size_t&, const double&, size_t) const noexcept;
    void initializeNetwork(std::vector<Edge>&, size_t&, std::vector<size_t>&) const noexcept;
    void growNetwork(std::vector<Edge>&, size_t&, std::vector<size_t>&, const size_t&, const double&) const noexcept;
    void sortNetwork(std::vector<Edge>&, const std::vector<size_t>&, const size_t&) const noexcept;

};


/// A container for all constant genetic features of the species we are simulating
class GeneticArchitecture {

private:

    size_t nChromosomes;
    size_t nTraits;
    size_t nLoci;
    double effectSizeShape;
    double effectSizeScale;
    std::vector<size_t> nLociPerTrait;
    std::vector<size_t> nEdgesPerTrait;
    std::vector<double> skewnesses;

    // A vector of chromosome sizes
    std::vector<double> chromosomeSizes;

    // A map of the gene regulatory network for each trait
    std::vector<Network> traitNetworkMaps;

    // A series of layers of locus-specific features across the genome
    Genome genome;

    // A set of vectors of loci underlying each trait
    //std::vector<std::vector<size_t> > traitUnderlyingLoci;

    /// Makers
    std::vector<double> makeChromosomeSizes() const noexcept;
    std::vector<Network> makeTraitNetworkMaps() const noexcept;
    Genome makeGenome() const noexcept;

public:

    /// Constructor
    GeneticArchitecture(const ParameterSet&);

    /// Getters
    std::vector<double> getChromosomeSizes() const { return chromosomeSizes; }
    std::vector<Network> getTraitNetworkMaps() const { return traitNetworkMaps; }

    /*
    // High-level functions
    void generateGeneticArchitecture(const ParameterSet&);
    void loadGeneticArchitecture(const ParameterSet&);
    void storeGeneticArchitecture(const ParameterSet&);

    // Low-level functions

    // Setting
    void createRecombinationMap(const ParameterSet&);
    void setChromosomeSizes(const size_t&); // obsolete
    void sampleGeneLocations(const ParameterSet&);
    void assignPhenotypes(const ParameterSet&);
    void sampleEffectSizes(const ParameterSet&);
    void sampleDominanceCoeff(const ParameterSet&);
    void makeRegulatoryNetworks(const ParameterSet&);
    void preferentialAttachmentNetwork(const size_t&, size_t&, const double&);
    void initializeNetwork(size_t&, std::vector<size_t>&);
    void growNetwork(const size_t&, size_t&, std::vector<size_t>&, const double&);
    void sortNetwork(const size_t&, const std::vector<size_t>&);
    void sampleInteractions(const ParameterSet&, const size_t&, const size_t&);  // Per phenotypic character

    // Reading
    bool validateArchitecture(std::ifstream&, const ParameterSet&);
    void loadChromosomeSizes(std::ifstream&);
    void loadLocusConstants(std::ifstream&, const ParameterSet&);
    void loadEpistaticInteractions(std::ifstream&);

    // Writing
    void writeChromosomeSizes(std::ofstream&);
    void writeLocusConstants(std::ofstream&, const ParameterSet&);
    void writeEpistaticInteractions(std::ofstream&, const ParameterSet&);

    std::vector<double> getGenomicLocations(const size_t&);
    std::vector<double> getChromosomeSizes(const size_t&);
    std::vector<size_t> getEncodedTraits(const size_t&, const std::vector<size_t>&);
    std::vector<double> getEffectSizes(const size_t&, const double&, const double&, const size_t&);
    std::vector<double> getDominanceCoeffs(const size_t&, const size_t&);
     */

};



#endif //EXPLICITGENOMESPECIATION_GENETICARCHITECTURE_H