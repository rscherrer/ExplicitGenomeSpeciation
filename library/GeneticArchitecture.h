#ifndef EXPLICITGENOMESPECIATION_GENETICARCHITECTURE_H
#define EXPLICITGENOMESPECIATION_GENETICARCHITECTURE_H

#include "ParameterSet.h"
#include <vector>
#include <list>

typedef std::pair<size_t, size_t> Edge;  // A network edge is a pair of locus indices

class GeneticArchitecture {

private:



public:


    // Constructor
    // It is a bad habit to have declared but non initialized members of a class
    // It is recommended to initialize the members using the member initialization list of the constructor
    // But genetic architecture is a large class with many members (many genes...)
    // So it would take a lot of computation to initialize the members by generating a new architecture,
    // If afterwards I end up erasing it all to replace it with the architecture that is given by some input file
    // A better alternative might be to have different constructors, one for generating a new architecture,
    // and the other to load an architecture from a file
    // Each constructor will have its own member initialization list
    // NB: functions to initialize variables can be used in the initialization list, initialization needs not be
    // assigning single values to the variables

    // Default constructor to generate a new architecture
    GeneticArchitecture(const size_t&, const size_t&);

    // Constructor to load an architecture from a file
    GeneticArchitecture(const std::string&);

    // Fields

    // A vector of chromosome sizes
    std::vector<double> chromosomeSizes;

    // A vector of genomic locations for all loci
    std::vector<double> locusGenomicLocations;

    // A vector of traits encoded by all loci
    std::vector<size_t> locusEncodedTraits;

    // A vector of effet sizes for all loci
    std::vector<double> locusEffectSizes;

    // A vector of dominance coefficients for all loci
    std::vector<double> locusDominanceCoeffs;

    // A vector of epistatic interactions for all loci (interacting partner + interaction weight)
    std::vector<std::vector<std::pair<size_t, double> > > locusInteractions;

    // Three vectors indicating what loci are underlying each trait
    std::vector<std::vector<size_t> > networkVertices;

    // A vector of locus indices
    std::vector<size_t> loci;

    // A vector of edges of the gene regulatory network (edges are pairs of locus indices)
    std::vector<Edge> edges;

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
    std::vector<double> getEffectSizes(const size_t&, const double&, const double&);

};

// Accessory functions
bool edgeCompare (const Edge&, const Edge&);






/// Constructor to generate a new genetic architecture
GeneticArchitecture::GeneticArchitecture(const size_t &nChromosomes, const size_t &nLoci, const size_t &nTraits,
        const std::vector<size_t> &nVertices, const double &shapeEffectSizes, const double &scaleEffectSizes) :
chromosomeSizes(getChromosomeSizes(nChromosomes)),
locusGenomicLocations(getGenomicLocations(nLoci)),
locusEncodedTraits(getEncodedTraits(nTraits, nVertices)),
locusEffectSizes(getEffectSizes(nLoci, shapeEffectSizes, scaleEffectSizes))//,
//locusDominanceCoeffs(),
//locusInteractions()

//networkVertices(getNetworkVertices()),
//loci(getLoci()),
//edges(getEdges())
{

}


/// Constructor to load a genetic architecture from a file
GeneticArchitecture::GeneticArchitecture(const std::string &architectureFileName)
{

}

#endif //EXPLICITGENOMESPECIATION_GENETICARCHITECTURE_H