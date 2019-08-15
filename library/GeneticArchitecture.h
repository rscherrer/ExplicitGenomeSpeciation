#ifndef EXPLICITGENOMESPECIATION_GENETICARCHITECTURE_H
#define EXPLICITGENOMESPECIATION_GENETICARCHITECTURE_H

#include "ParameterSet.h"
#include "Random.h"
#include "Genome.h"
#include "Network.h"
#include <vector>


typedef std::pair<size_t, size_t> Edge;
typedef std::vector<Network> MultiNet;

/// A container for all constant genetic features
class GeneticArchitecture {

public:

    GeneticArchitecture(const ParameterSet&);

    /// Getters
    std::vector<double> getChromosomes() const { return chromosomeSizes; }
    Genome getGenome() const { return genome; }
    MultiNet getNetworks() const { return traitNetworks; }

private:

    size_t nChromosomes;
    size_t nLoci;
    std::vector<size_t> nLociPerTrait;
    std::vector<size_t> nEdgesPerTrait;
    std::vector<double> skewnesses;
    double effectSizeShape;
    double effectSizeScale;
    double interactionWeightShape;
    double interactionWeightScale;

    std::vector<double> chromosomeSizes;
    Genome genome;
    std::vector<Network> traitNetworks;

    // A set of vectors of loci underlying each trait
    //std::vector<std::vector<size_t> > traitUnderlyingLoci;

    /// Makers
    std::vector<double> makeChromosomes();
    Genome makeGenome();
    std::vector<Network> makeNetworks();

};

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
void growNetwork(const size_t&, size_t&, std::vector<size_t>&,
const double&);
void sortNetwork(const size_t&, const std::vector<size_t>&);
void sampleInteractions(const ParameterSet&, const size_t&, const size_t&);
// Per phenotypic character

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
std::vector<size_t> getEncodedTraits(const size_t&,
const std::vector<size_t>&);
std::vector<double> getEffectSizes(const size_t&, const double&,
const double&, const size_t&);
std::vector<double> getDominanceCoeffs(const size_t&, const size_t&);
 */

#endif //EXPLICITGENOMESPECIATION_GENETICARCHITECTURE_H
