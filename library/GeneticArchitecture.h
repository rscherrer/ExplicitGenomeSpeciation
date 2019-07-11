#ifndef EXPLICITGENOMESPECIATION_GENETICARCHITECTURE_H
#define EXPLICITGENOMESPECIATION_GENETICARCHITECTURE_H

#include "ParameterSet.h"
#include <vector>
#include <list>

typedef std::pair<size_t, size_t> Edge;  // A network edge is a pair of locus indices

class GeneticArchitecture {

private:

    // A vector of chromosome sizes
    std::vector<double> chromosomeSizes;

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
    GeneticArchitecture(const size_t&);

    // Constructor to load an architecture from a file
    GeneticArchitecture(const std::string&);

    // Fields

    // A vector of genomic locations for all loci
    std::vector<double> locusGenomicLocations;

    // A vector of traits encoded by all loci
    std::vector<size_t> locusEncodedTraits;

    // A vector of effet sizes for all loci
    std::vector<double> locusEffectSizes;

    // A vector of dominance coefficients for all loci
    std::vector<double> locusDominanceCoeffs;

    // A vector of epistatic interactions for all loci (interacting partner + interaction weight)
    std::vector<Edge> locusInteractions;

    std::vector<std::vector<size_t> > networkVertices;
    std::vector<size_t> loci;
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

};

// Accessory functions
bool edgeCompare (const Edge&, const Edge&);


/// Function to get a vector of chromosome lengths
std::vector<double> getChromosomeSizes(const size_t &nChromosomes)
{
    std::vector<double> chromosomeSizes;

    // Chromosomes all have the same size
    for (size_t i = 0u; i < nChromosomes; ++i) {
        chromosomeSizes.push_back((i + 1.0) / nChromosomes);
    }

    return chromosomeSizes;
}




/// Constructor to generate a new genetic architecture
GeneticArchitecture::GeneticArchitecture(const size_t &nChromosomes) //:
//chromosomeSizes(getChromosomeSizes(nChromosomes))//,
//locusGenomicLocations(getGenomicLocations()),
//locusEncodedTraits(getEncodedTraits()),
//locusEffectSizes(),
//locusDominanceCoeffs(),
//locusInteractions()
//locusConstants(getLocusConstants()),
//networkVertices(getNetworkVertices()),
//loci(getLoci()),
//edges(getEdges())
{
    chromosomeSizes = getChromosomeSizes(nChromosomes);
}


/// Constructor to load a genetic architecture from a file
GeneticArchitecture::GeneticArchitecture(const std::string &architectureFileName)
{

}

#endif //EXPLICITGENOMESPECIATION_GENETICARCHITECTURE_H