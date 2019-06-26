#ifndef EXPLICITGENOMESPECIATION_GENETICARCHITECTURE_H
#define EXPLICITGENOMESPECIATION_GENETICARCHITECTURE_H

#include "ParameterSet.h"
#include <vector>
#include <list>

typedef std::pair<size_t, size_t> Edge;  // A network edge is a pair of locus indices

class GeneticArchitecture {

public:

    struct LocusConstants {

        size_t trait;
        size_t chromosome;
        double location;
        double effectSize;
        double dominanceCoeff;
        std::list<std::pair<size_t, double> > neighbors;

    };

    // Fields
    std::vector<double> chromosomeSizes;
    std::vector<LocusConstants> locusConstants;
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
    void setChromosomeSizes(const size_t&);
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


#endif //EXPLICITGENOMESPECIATION_GENETICARCHITECTURE_H