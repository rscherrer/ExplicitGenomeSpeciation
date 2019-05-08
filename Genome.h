//
// Created by p278834 on 8-5-2019.
//

#ifndef EXPLICITGENOMESPECIATION_GENOME_H
#define EXPLICITGENOMESPECIATION_GENOME_H

#include <vector>
#include <array>
#include <list>
#include <set>
#include "ParameterSet.h"
#include "Population.h"

bool edgeCompare (const Edge &x, const Edge &y) {
    if(x.first == y.first) return (x.second < y.second);
    else return (x.first < y.first);
}

class Genome {

public:

    // A Character object is locus-specific and population-wide
    struct Character
    {
        size_t character, linkageGroup;
        double location, effectSize, dominanceCoeff, avgEffectOfSubstitution,
                varD, F_it, F_is, F_st, P_st, G_st, Q_st, C_st;
        std::array<double, 3u> alleleFrequency, meanEffect, varP, varG, varA, varI;
        std::list<std::pair<size_t, double> > edges;
    };

    typedef std::pair<size_t, size_t> Edge;

    // These are genome specific
    std::vector<double> chromosomeSize; // size nchromosomes - 1
    std::vector<std::set<size_t> > vertices; // size nCharacter
    std::vector<Character> characterLocus; // size nLoci

    // Setters
    void generateGeneticArchitecture(const ParameterSet&);
    void storeGeneticArchitecture(const std::string&, const ParameterSet&);
    void loadGeneticArchitecture(const std::string&, const ParameterSet&, const Population&);

    void preferentialAttachmentNetwork(const size_t n, size_t e, const double exponent, std::vector<Edge> &edges);


private:


};


#endif //EXPLICITGENOMESPECIATION_GENOME_H
