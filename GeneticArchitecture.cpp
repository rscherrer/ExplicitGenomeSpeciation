#include "GeneticArchitecture.h"
#include "ParameterSet.h"
#include "random.h"
#include "utils.h"
#include <vector>
#include <iostream>
#include <algorithm>

typedef std::pair<size_t, size_t> Edge;  // A network edge is a pair of locus indices

// Accessory function
bool edgeCompare (const Edge &x, const Edge &y)
{
    if (x.first == y.first) {
        return (x.second < y.second);
    }
    else {
        return (x.first < y.first);
    }
}

// Constructor
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

// Generate a genetic architecture

// High-level function
void GeneticArchitecture::generateGeneticArchitecture(const ParameterSet& parameters)
{
    std::clog << "Generating a new genetic architecture\n";

    // Recombination map
    std::clog << "  Placing loci across the genome.";
    createRecombinationMap(parameters);
    std::clog << "..done\n";

    // Assign phenotypic traits
    std::clog << "  Assigning loci to traits.";
    assignPhenotypes(parameters);
    std::clog << "..done\n";

    // Assign additive and dominance effects
    std::clog << "  Sampling additive and dominance effects.";
    sampleEffectSizes(parameters);
    sampleDominanceCoeff(parameters);
    std::clog << "..done\n";

    // Make a regulatory network
    std::clog << "  Creating gene interaction network.";
    makeRegulatoryNetworks(parameters);
    std::clog << ".done\n";

}

// Low-level functions
void GeneticArchitecture::createRecombinationMap(const ParameterSet &parameters)
{


    setChromosomeSizes(parameters.nChromosomes);


    sampleGeneLocations(parameters);

}

void GeneticArchitecture::setChromosomeSizes(const size_t &nChromosomes)
{

    // Chromosomes all have the same size
    for (size_t i = 0u; i < nChromosomes; ++i) {
        chromosomeSizes.push_back((i + 1.0) / nChromosomes);
    }

}

void GeneticArchitecture::sampleGeneLocations(const ParameterSet &parameters)
{

    // Genes are distributed uniformly across the genome
    std::vector<double> geneLocations;
    for (size_t i = 0u; i < parameters.nLoci; ++i) {
        geneLocations[i] = rnd::uniform();
    }
    std::sort(geneLocations.begin(), geneLocations.end());
    for (size_t i = 0u, lg = 0u; i < parameters.nLoci; ++i) {
        while (lg < parameters.nChromosomes - 1u && geneLocations[i] > chromosomeSizes[lg]) {
            ++lg;
        }
        locusConstants[i].chromosome = lg;
        locusConstants[i].location = geneLocations[i];
    }

}

void GeneticArchitecture::assignPhenotypes(const ParameterSet &parameters)
{

    std::vector<size_t> nVertices {parameters.nEcoLoci, parameters.nMatLoci, parameters.nNtrLoci};

    // Randomized vector of locus indices
    for (size_t i = 0u; i < parameters.nLoci; ++i) {
        loci[i] = i;
    }
    std::shuffle(loci.begin(), loci.end(), rnd::rng);

    // For each phenotypic trait
    for (size_t crctr = 0u, k = 0u; crctr < parameters.nTraits; ++crctr) {

        // Assign loci to the trait
        for (size_t j = 0u; j < nVertices[crctr]; ++j, ++k) {
            networkVertices[crctr].push_back(loci[k]);
            locusConstants[loci[k]].trait = crctr;
        }
    }

}

void GeneticArchitecture::sampleEffectSizes(const ParameterSet &parameters)
{
    // For each phenotypic trait
    for (size_t crctr = 0u; crctr < parameters.nTraits; ++crctr) {
        double sumsqAdditive = 0.0;

        // For each vertex in the network
        for (size_t i : networkVertices[crctr]) {

            // Sample additive effect size from a bilateral Gamma distribution
            locusConstants[i].effectSize = std::gamma_distribution<double>(parameters.shapeEffectSizes, parameters.scaleEffectSizes)(rnd::rng);
            if (rnd::bernoulli(0.5)) {
                locusConstants[i].effectSize *= -1.0;
            }
            sumsqAdditive += sqr(locusConstants[i].effectSize);
        }

        // Normalize
        sumsqAdditive = sumsqAdditive > 0.0 ? sqrt(sumsqAdditive) : 1.0;
        for (size_t i : networkVertices[crctr]) {
            locusConstants[i].effectSize /= sumsqAdditive;
        }
    }
}

void GeneticArchitecture::sampleDominanceCoeff(const ParameterSet &parameters)
{
    // For each phenotypic trait
    for (size_t crctr = 0u; crctr < parameters.nTraits; ++crctr) {
        double sumsqDominance = 0.0;

        // For each vertex in the network
        for (size_t i : networkVertices[crctr]) {

            // Sample dominance coefficient
            locusConstants[i].dominanceCoeff = fabs(rnd::normal(0.0, 1.0));
            sumsqDominance += sqr(locusConstants[i].dominanceCoeff);
        }

        // Normalize
        sumsqDominance = sumsqDominance > 0.0 ? sqrt(sumsqDominance) : 1.0;
        for (size_t i : networkVertices[crctr]) {
            locusConstants[i].dominanceCoeff /= sumsqDominance;
        }
    }
}

void GeneticArchitecture::makeRegulatoryNetworks(const ParameterSet &parameters)
{

    std::vector<size_t> nEdges {parameters.nEcoInteractions, parameters.nMatInteractions, parameters.nNtrInteractions};
    std::vector<size_t> nVertices {parameters.nEcoLoci, parameters.nMatLoci, parameters.nNtrLoci};

    // For each phenotypic trait
    for (size_t crctr = 0u, offset = 0u; crctr < parameters.nTraits; ++crctr) {

        // Generate edges with the preferential attachment algorithm
        preferentialAttachmentNetwork(nVertices[crctr], nEdges[crctr], parameters.networkSkewness);

        // Sample interaction partners and weights for the current trait
        sampleInteractions(parameters, crctr, offset);

        // Update the offset
        offset += nVertices[crctr];
    }

}

void GeneticArchitecture::preferentialAttachmentNetwork(const size_t &nVertices, size_t &nEdges, const double &skewness)
{

    if (!(nVertices > 1u && skewness > 0.0)) {
        throw std::runtime_error("Invalid parameters in GeneticArchitecture::preferentialAttachmentNetwork()");
    }

    if (nEdges == 0u) {
        return;
    }

    // Start the network
    std::vector<size_t> degrees(nVertices, 0u);
    initializeNetwork(nVertices, nEdges, degrees);

    // Grow network by linking preferentially to well-connected nodes
    std::clog << ':';
    growNetwork(nVertices, nEdges, degrees, skewness);

    // Relabel node indices after sorting with respect to degree
    std::clog << ':';
    sortNetwork(nVertices, degrees);

    std::clog << ':';
}

void GeneticArchitecture::initializeNetwork(const size_t &nVertices, size_t &nEdges, std::vector<size_t> &degrees)
{

    if (!edges.empty()) {
        edges.clear();
    }

    // Create initial network
    edges.emplace_back(Edge {0u, 1u});
    degrees[0u] = degrees[1u] = 1u;
    --nEdges;

}

void GeneticArchitecture::growNetwork(const size_t &nVertices, size_t &nEdges, std::vector<size_t> &degrees, const double &skewness)
{

    for (size_t i = 2u; i < nVertices && nEdges > 0u; ++i) {
        std::vector<double> w(i);
        for (size_t j = 0u; j < i; ++j) {
            w[j] = pow(degrees[j], skewness);
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

}

void GeneticArchitecture::sortNetwork(const size_t &nVertices, const std::vector<size_t> &degrees)
{

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

}

void GeneticArchitecture::sampleInteractions(const ParameterSet &parameters, const size_t &crctr, const size_t &offset)
{

    double sumsqWeights = 0.0;  // For later normalization

    // For each edge in the network
    for (Edge &edge : edges) {

        // Select interacting partners
        size_t i = loci[offset + edge.first];
        size_t j = loci[offset + edge.second];

        // Make sure that both partner genes underlie the current phenotypic trait
        bool isSametrait = locusConstants[i].trait == crctr && locusConstants[j].trait == crctr;
        if (!isSametrait) {
            throw std::logic_error("Invalid epistatic interaction in GeneticArchitecture::generateGeneticArchitecture()");
        }

        // Sample interaction weight
        double interactionWeight = std::gamma_distribution<double>(parameters.shapeInteractionWeights, parameters.scaleInteractionWeights)(rnd::rng);
        if (rnd::bernoulli(0.5)) {
            interactionWeight *= -1.0;
        }

        // Record locus j as a partner of locus i
        locusConstants[i].neighbors.emplace_back(std::make_pair(j, interactionWeight));

        // Update sum of squares
        sumsqWeights += sqr(interactionWeight);
    }

    // Normalize interaction weights
    sumsqWeights = sumsqWeights > 0.0 ? sqrt(sumsqWeights) : 1.0;
    for (size_t i : networkVertices[crctr]) {
        for (std::pair<size_t, double> &edge : locusConstants[i].neighbors) {
            edge.second /= sumsqWeights;
        }
    }

}

// Load a genetic architecture

// High-level function
void GeneticArchitecture::loadGeneticArchitecture(const ParameterSet &parameters)
{
    std::clog << "Loading genetic architecture from file " << parameters.architectureFileName << '\n';

    // Open genetic architecture file
    std::ifstream ifs(parameters.architectureFileName);
    if (!ifs.is_open()) {
        throw std::runtime_error("Unable to open file in GeneticArchitecture::loadGeneticArchitecture()");
    }

    // Make sure the loaded architecture is valid
    std::clog << "  Validation.";
    bool isValidArchitecture = validateArchitecture(ifs, parameters);
    if (!isValidArchitecture) {
        throw std::logic_error("Genetic architecture in file is incompatible with current settings in GeneticArchitecture::loadGeneticArchitecture()");
    }
    std::clog << "..done\n";

    // The order of the following is important!

    // Chromosome sizes
    std::clog << "  Loading recombination map.";
    loadChromosomeSizes(ifs);
    std::clog << "..done\n";

    // Locus properties
    std::clog << "  Loading locus properties.";
    loadLocusConstants(ifs, parameters);
    std::clog << "..done\n";

    // Epistatic interactions
    std::clog << "  Loading epistatic interactions.";
    loadEpistaticInteractions(ifs);
    std::clog << "..done\n";
}

// Low-level functions
bool GeneticArchitecture::validateArchitecture(std::ifstream &ifs, const ParameterSet &parameters)
{

    // Read genetic architecture details from file
    size_t tmpEcoLoci;
    size_t tmpMatLoci;
    size_t tmpNtrLoci;
    size_t tmpEcoInteractions;
    size_t tmpMatInteractions;
    size_t tmpChromosomes;

    ifs >> tmpEcoLoci >> tmpMatLoci >> tmpNtrLoci >> tmpEcoInteractions >> tmpMatInteractions >> tmpChromosomes;

    bool isValid = tmpEcoLoci == parameters.nEcoLoci &&
                   tmpMatLoci == parameters.nMatLoci &&
                   tmpNtrLoci == parameters.nNtrLoci &&
                   tmpEcoInteractions == parameters.nEcoInteractions &&
                   tmpMatInteractions == parameters.nMatInteractions &&
                   tmpChromosomes == parameters.nChromosomes;

    return isValid;

}

void GeneticArchitecture::loadChromosomeSizes(std::ifstream &ifs)
{
    for (double &x : chromosomeSizes) {
        ifs >> x;
    }
}

void GeneticArchitecture::loadLocusConstants(std::ifstream &ifs, const ParameterSet &parameters)
{
    for (size_t i = 0u, lg = 0u; i < parameters.nLoci; ++i) {
        size_t j;
        ifs >> j;
        ifs >> locusConstants[j].trait
            >> locusConstants[j].location
            >> locusConstants[j].effectSize
            >> locusConstants[j].dominanceCoeff;
        networkVertices[locusConstants[j].trait].push_back(j);

        while (lg < parameters.nChromosomes - 1u && locusConstants[j].location > chromosomeSizes[lg]) {
            ++lg;
        }
        locusConstants[i].chromosome = lg;
    }
}

void GeneticArchitecture::loadEpistaticInteractions(std::ifstream &ifs)
{
    size_t i, j;
    double w;
    while (ifs >> i >> j >> w) {
        locusConstants[i].neighbors.emplace_back(std::make_pair(j, w));
    }
}

// Store a genetic architecture

// High-level function
void GeneticArchitecture::storeGeneticArchitecture(const ParameterSet& parameters)
{
    // Open genetic architecture file
    std::ofstream ofs(parameters.architectureFileName);
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
    writeChromosomeSizes(ofs);

    // Gene location, effect size and dominance coefficient
    writeLocusConstants(ofs, parameters);

    // Epistatic interactions
    writeEpistaticInteractions(ofs, parameters);
}

// Low-level functions
void GeneticArchitecture::writeChromosomeSizes(std::ofstream &ofs)
{
    for (double x : chromosomeSizes) {
        ofs << x << '\t';
    }
    ofs << '\n';
}

void GeneticArchitecture::writeLocusConstants(std::ofstream &ofs, const ParameterSet &parameters)
{
    for (size_t i = 0u; i < parameters.nLoci; ++i) {
        ofs << i << '\t'
            << locusConstants[i].trait << '\t'
            << locusConstants[i].location << '\t'
            << locusConstants[i].effectSize << '\t'
            << locusConstants[i].dominanceCoeff << '\n';
    }
}

void GeneticArchitecture::writeEpistaticInteractions(std::ofstream &ofs, const ParameterSet &parameters)
{
    for (size_t i = 0u; i < parameters.nLoci; ++i) {
        for (const std::pair<size_t, double> &edge : locusConstants[i].neighbors) {
            ofs << i << '\t' << edge.first << '\t' << edge.second << '\n';
        }
    }
}


























