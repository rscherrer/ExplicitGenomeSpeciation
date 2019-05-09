//
// Created by p278834 on 9-5-2019.
//

#include "GeneticArchitecture.h"
#include <iostream>
#include <algorithm>
#include <fstream>
#include "random.h"
#include "square.h"

typedef std::pair<size_t, size_t> Edge;  // A network edge is a pair of locus indices

GeneticArchitecture::GeneticArchitecture(const ParameterSet &parameters)
{
    if (parameters.isGenerateArchitecture) {
        generateGeneticArchitecture(parameters);
        storeGeneticArchitecture(parameters);
    }
    else {
        loadGeneticArchitecture(parameters);
    }
}

void GeneticArchitecture::generateGeneticArchitecture(const ParameterSet& parameters)
{
    std::clog << "Generating a new genetic architecture\n";

    // Recombination map
    std::clog << "  Placing loci across the genome.";

    // Chromosome all have the same size
    for (size_t i = 0u; i < parameters.nChromosomes - 1u; ++i) {
        chromosomeSizes[i] = (i + 1.0) / parameters.nChromosomes;
    }

    // Genes are distributed uniformly across genome
    std::vector<double> geneLocations;
    for (size_t i = 0u; i < parameters.nLoci; ++i) {
        geneLocations[i] = rnd::uniform();
    }
    std::sort(geneLocations.begin(), geneLocations.end());
    for (size_t i = 0u, lg = 0u; i < parameters.nLoci; ++i) {
        while(lg < parameters.nChromosomes - 1u && geneLocations[i] > chromosomeSizes[lg]) ++lg;
        locusConstants[i].linkageGroup = lg;
        locusConstants[i].location = geneLocations[i];
    }
    std::clog << "..done\n";

    // Randomize gene sequence and assign loci to characters
    std::clog << "  Assigning loci to characters.";
    std::vector<size_t> seq;
    for (size_t i = 0u; i < parameters.nLoci; ++i) {
        seq[i] = i;
    }
    std::shuffle(seq.begin(), seq.end(), rnd::rng);
    std::vector<size_t> nVtx {parameters.nEcoLoci, parameters.nMatLoci, parameters.nNtrLoci};
    for (size_t crctr = 0u, k = 0u; crctr < parameters.nCharacter; ++crctr) {
        for (size_t j = 0u; j < nVtx[crctr]; ++j, ++k) {
            networkVertices[crctr].insert(seq[k]);
            locusConstants[seq[k]].character = crctr;
        }
    }
    std::clog << "..done\n";

    // Assign additive and dominance effects
    std::clog << "  Sampling additive and dominance effects.";
    for (size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr) {
        double sumaa = 0.0, sumhh = 0.0;
        for (size_t i : networkVertices[crctr]) {
            locusConstants[i].effectSize = std::gamma_distribution<double>(parameters.alphaAdditive, 1.0)(rnd::rng);
            if(rnd::bernoulli(0.5)) locusConstants[i].effectSize *= -1.0;
            sumaa += sqr(locusConstants[i].effectSize);
            locusConstants[i].dominanceCoeff = fabs(rnd::normal(0.0, 1.0));
            sumhh += sqr(locusConstants[i].dominanceCoeff);
        }
        sumaa = sumaa > 0.0 ? sqrt(sumaa) : 1.0;
        sumhh = sumhh > 0.0 ? sqrt(sumhh) : 1.0;
        for (size_t i : networkVertices[crctr]) {
            locusConstants[i].effectSize /= sumaa;
            locusConstants[i].dominanceCoeff /=sumhh;
        }
    }
    std::clog << "..done\n";

    // Make a regulatory network
    std::clog << "  Creating gene interaction network.";
    std::vector<size_t> nEdg {parameters.nEcoInteractions, parameters.nMatInteractions, parameters.nNtrInteractions};

    // For each phenotypic character
    for (size_t crctr = 0u, offset = 0u; crctr < parameters.nCharacter; ++crctr) {

        std::vector<Edge> edges = preferentialAttachmentNetwork(nVtx[crctr], nEdg[crctr], parameters.networkSkewness);

        // Map network vertices to loci
        double sumww = 0.0;
        for (Edge &edge : edges) {

            size_t i = seq[offset + edge.first];
            size_t j = seq[offset + edge.second];
            if (!(locusConstants[i].character == crctr &&
                    locusConstants[j].character == crctr)) {
                throw std::logic_error("Invalid epistatic interaction in Individual::generateGeneticArchitecture()");
            }

            // Assign interaction effect sizes
            double wij = std::gamma_distribution<double>(parameters.alphaInteraction, 1.0)(rnd::rng);
            if (rnd::bernoulli(0.5)) {
                wij *= -1.0;
            }
            locusConstants[i].edges.emplace_back(std::make_pair(j, wij));
            sumww += wij * wij;
        }

        // Normalise interaction weights
        sumww = sumww > 0.0 ? sqrt(sumww) : 1.0;
        for (size_t i : networkVertices[crctr]) {
            for (std::pair<size_t, double> &edge : locusConstants[i].edges) {
                edge.second /= sumww;
            }
        }

        // Update the offset
        offset += nVtx[crctr];
    }

    std::clog << ".done\n";

}

std::vector<Edge> GeneticArchitecture::preferentialAttachmentNetwork(const size_t &nVertices, size_t &nEdges, const double &exponent)
{

    if (!(nVertices > 1u && exponent > 0.0)) {
        throw std::runtime_error("Invalid parameters in preferentialAttachmentNetwork()");
    }

    std::vector<Edge> edges;

    if (!edges.empty()) {
        edges.clear();
    }
    if (nEdges == 0u) {
        return edges;
    }

    // Create initial network
    edges.emplace_back(Edge {0u, 1u});
    std::vector<size_t> degrees(nVertices, 0u);
    degrees[0u] = degrees[1u] = 1u;
    --nEdges;

    // Grow network by linking preferentially to well-connected nodes
    std::clog << ':';
    for (size_t i = 2u; i < nVertices && nEdges > 0u; ++i) {
        std::vector<double> w(i);
        for (size_t j = 0u; j < i; ++j) {
            w[j] = pow(degrees[j], exponent);
        }
        size_t ki = (i == nVertices - 1u ? nEdges : rnd::binomial(nEdges, 1.0 / (nVertices - i)));

        while (ki) {

            // Sample neighbours without replacement
            double sumw = 0.0;
            for (size_t j = 0u; j < i; ++j) {
                sumw += w[j];
            }
            if (sumw < 1.0) {
                break;
            }
            std::discrete_distribution<size_t> prob(w.begin(), w.end());
            size_t j = prob(rnd::rng);
            edges.emplace_back(Edge {i, j});
            w[j] = 0.0;
            ++degrees[j];
            ++degrees[i];
            --ki;
            --nEdges;
        }
    }

    // Relabel node indices after sorting with respect to degree
    std::clog << ':';
    std::vector<size_t> ranks(nVertices, 0u);
    for (size_t i = 0u; i < nVertices - 1u; ++i) {
        for (size_t j = i + 1u; j < nVertices; ++j) {
            if (degrees[j] > degrees[i]) {
                ++ranks[i];
            }
            else {
                ++ranks[j];
            }
        }
    }
    for (Edge &edge : edges) {
        edge.first = ranks[edge.first];
        edge.second = ranks[edge.second];
        if (edge.first > edge.second) {
            std::swap(edge.first, edge.second);
        }
    }
    std::sort(edges.begin(), edges.end(), edgeCompare);
    std::clog << ':';

    return edges;
}

void GeneticArchitecture::loadGeneticArchitecture(const ParameterSet& parameters)
{
    std::clog << "Loading genetic architecture from file " << parameters.architectureFilename << '\n';

    // Open genetic architecture file
    std::ifstream ifs(parameters.architectureFilename);
    if (!ifs.is_open()) {
        throw std::runtime_error("Unable to open file in loadGeneticArchitecture()");
    }
    std::clog << "  Validation.";

    size_t tmpEcoLoci, tmpMatLoci, tmpNtrLoci, tmpEcoInteractions, tmpMatInteractions, tmpChromosomes;
    ifs >> tmpEcoLoci >> tmpMatLoci >> tmpNtrLoci >> tmpEcoInteractions >> tmpMatInteractions >> tmpChromosomes;

    // Validation
    bool isCompatible = tmpEcoLoci == parameters.nEcoLoci && 
                        tmpMatLoci == parameters.nMatLoci && 
                        tmpNtrLoci == parameters.nNtrLoci &&
                        tmpEcoInteractions == parameters.nEcoInteractions && 
                        tmpMatInteractions == parameters.nMatInteractions &&
                        tmpChromosomes == parameters.nChromosomes;
    if (!isCompatible) {
        throw std::logic_error("Genetic architecture in file is incompatible with current settings in loadGeneticArchitecture()");
    }
    std::clog << "..done\n";

    // Chromosome sizes
    std::clog << "  Loading recombination map.";
    for(double &x : chromosomeSizes) ifs >> x;
    std::clog << "..done\n";

    // Locus properties
    std::clog << "  Loading locus properties.";
    for (size_t i = 0u, lg = 0u; i < parameters.nLoci; ++i) {
        size_t j;
        ifs >> j;
        ifs >> locusConstants[j].character
            >> locusConstants[j].location
            >> locusConstants[j].effectSize
            >> locusConstants[j].dominanceCoeff;
        networkVertices[locusConstants[j].character].insert(j);

        while (lg < parameters.nChromosomes - 1u && locusConstants[j].location > chromosomeSizes[lg]) {
            ++lg;
        }
        locusConstants[i].linkageGroup = lg;
    }
    std::clog << "..done\n";

    // Epistatic interactions
    std::clog << "  Loading epistatic interactions.";
    size_t i, j;
    double w;
    while (ifs >> i >> j >> w) {
        locusConstants[i].edges.emplace_back(std::make_pair(j, w));
    }
    std::clog << "..done\n";
}

void GeneticArchitecture::storeGeneticArchitecture(const ParameterSet& parameters)
{
    // Open genetic architecture file
    std::ofstream ofs(parameters.architectureFilename);
    if (!ofs.is_open()) {
        throw std::runtime_error("unable to open file in storeGeneticArchitecture()");
    }

    ofs << parameters.nEcoLoci << '\n'
        << parameters.nMatLoci << '\n'
        << parameters.nNtrLoci << '\n'
        << parameters.nEcoInteractions << '\n'
        << parameters.nMatInteractions << '\n'
        << parameters.nChromosomes << '\n';

    // Chromosome sizes
    for (double x : chromosomeSizes) {
        ofs << x << '\t';
    }
    ofs << '\n';

    // Gene location, effect size and dominance coefficient
    for (size_t i = 0u; i < parameters.nLoci; ++i) {
        ofs << i << '\t'
            << locusConstants[i].character << '\t'
            << locusConstants[i].location << '\t'
            << locusConstants[i].effectSize << '\t'
            << locusConstants[i].dominanceCoeff << '\n';   
    }

    // Epistatic interactions
    for (size_t i = 0u; i < parameters.nLoci; ++i) {
        for (const std::pair<size_t, double> &edge : locusConstants[i].edges) {
            ofs << i << '\t' << edge.first << '\t' << edge.second << '\n';
        }
    }
    
}