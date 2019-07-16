#include "GeneticArchitecture.h"
#include "ParameterSet.h"
#include "random.h"
#include "utils.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>


/// A network edge is a pair of locus indices
typedef std::pair<size_t, size_t> Edge;


// Note:
// How about passing parameters to the genetic architecture constructor, so we don't have to pass each argument
// one by one?


/// Constructor of genetic architecture
GeneticArchitecture::GeneticArchitecture(const ParameterSet &pars) :
    nTraits(pars.getNTraits()),
    nChromosomes(pars.getNChromosomes()),
    nLociPerTrait(pars.getNLociPerTrait()),
    nEdgesPerTrait(pars.getNEdgesPerTrait()),
    skewnesses(pars.getSkewnesses()),
    chromosomeSizes(makeChromosomeSizes()),
    traitNetworkMaps(makeTraitNetworkMaps()),
    genome(makeGenome())
{}


/// Function to make a vector of chromosome sizes
std::vector<double> GeneticArchitecture::makeChromosomeSizes() const noexcept
{

    std::vector<double> chromsizes;

    // Chromosomes all have the same size
    for (size_t i = 0u; i < nChromosomes; ++i)
        chromsizes.push_back((i + 1.0) / nChromosomes);

    return chromsizes;

}


/// Function to make a vector of interacting partner loci for each trait
std::vector<Network> GeneticArchitecture::makeTraitNetworkMaps() const noexcept
{
    std::vector<Network> networks;

    // For each trait
    for (size_t trait = 0u; trait < nTraits; ++trait)
    {

        // Make a network map (a vector of edges) for the current trait using the preferential attachment algorithm
        Network network = Network(nLociPerTrait[trait], nEdgesPerTrait[trait], skewnesses[trait]);

        networks.push_back(network);
    }

    // The indices in these network maps are indices among the loci underlying a given trait,
    // not absolute loci indices across the genome

    assert(networks.size() == nTraits);
    for (size_t trait = 0u; trait < nTraits; ++trait)
        assert(networks[trait].map.size() == nEdgesPerTrait[trait]);

    return networks;
}


/// Network constructor
Network::Network(const size_t &nvertices, const size_t &nedges, const double &skew) :
        nVertices(nvertices),
        nEdges(nedges),
        skewness(skew),
        map(makeNetwork(nvertices, skew, nedges))
{}


/// Function to make a new interaction network based on the preferential attachment algorithm
std::vector<Edge> Network::makeNetwork(const size_t &nVertices, const double &skewness,
        size_t nEdges) const noexcept
{
    std::vector<Edge> network;

    assert(nVertices > 1u);
    assert(skewness > 0.0);

    if (nEdges == 0u) return network;

    // Start the network
    std::vector<size_t> degrees(nVertices, 0u);
    initializeNetwork(network, nEdges, degrees);

    // Grow network by linking preferentially to well-connected nodes
    growNetwork(network, nEdges, degrees, nVertices, skewness);

    // Relabel node indices after sorting with respect to degree
    sortNetwork(network, degrees, nVertices);

    return network;
}


/// Function to create a new network
void Network::initializeNetwork(std::vector<Edge> &network, size_t &nEdges, std::vector<size_t> &degrees)
const noexcept
{

    if (!network.empty()) network.clear();

    // Create initial network
    network.emplace_back(Edge {0u, 1u});
    degrees[0u] = degrees[1u] = 1u;
    --nEdges;

}


/// Function to iteratively grow a network based on the preferential attachment algorithm
void Network::growNetwork(std::vector<Edge> &network, size_t &nEdges,
        std::vector<size_t> &degrees, const size_t &nVertices, const double &skewness) const noexcept
{

    // For each vertex in the network...
    for (size_t vertex = 2u; vertex < nVertices && nEdges > 0u; ++vertex) {

        // Assign probability weights to all potential partners
        std::vector<double> weights(vertex);

        for (size_t partner = 0u; partner < vertex; ++partner) {
            weights[partner] = pow(degrees[partner], skewness);
        }

        // Sample the number of attachments of the current vertex
        size_t nAttachments = (vertex == nVertices - 1u ? nEdges : rnd::binomial(nEdges, 1.0 / (nVertices - vertex)));

        // For each new attachment...
        while (nAttachments) {

            // Compute the sum of weights, quit the loop if too small
            double sumWeights = 0.0;
            for (size_t partner = 0u; partner < vertex; ++partner)
                sumWeights += weights[partner];
            if (sumWeights < 1.0) break;

            // Sample the new partner without replacement
            std::discrete_distribution<size_t> attachmentProbs(weights.begin(), weights.end());
            size_t partner = attachmentProbs(rnd::rng);

            // Add the new partner to the network
            network.emplace_back(Edge {vertex, partner});
            weights[partner] = 0.0;
            ++degrees[partner];
            ++degrees[vertex];
            --nAttachments;
            --nEdges;
        }
    }
}


/// Function to compare edges in a network with respect to their degrees
bool edgeCompare(const Edge &x, const Edge &y) noexcept
{
    if (x.first == y.first) {
        return (x.second < y.second);
    }
    else {
        return (x.first < y.first);
    }
}


/// Function to sort the edges in the network by degree
void Network::sortNetwork(std::vector<Edge> &network, const std::vector<size_t> &degrees,
        const size_t &nVertices) const noexcept
{

    // Compute the ranks of all vertices with respect to their degrees
    std::vector<size_t> ranks(nVertices, 0u);
    for (size_t vertex = 0u; vertex < nVertices - 1u; ++vertex) {
        for (size_t othervertex = vertex + 1u; othervertex < nVertices; ++othervertex) {
            if (degrees[othervertex] > degrees[vertex]) {
                ++ranks[vertex];
            }
            else {
                ++ranks[othervertex];
            }
        }
    }

    // For each edge, swap partners to make sure the first partner always have the lower rank
    for (Edge &edge : network) {
        edge.first = ranks[edge.first];
        edge.second = ranks[edge.second];
        if (edge.first > edge.second) {
            std::swap(edge.first, edge.second);
        }
    }

    // Sort the edges
    std::sort(network.begin(), network.end(), edgeCompare);

}





/*
/// Function to make a set of layers of genetic features across all loci
Genome makeGenome()
{

}


/// Function to make a set of vectors of loci underlying each trait
std::vector<std::vector<size_t> > makeTraitUnderlyingLoci()
{

}

*/


//----------------------------------------------


/*
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
    initializeNetwork(nEdges, degrees);

    // Grow network by linking preferentially to well-connected nodes
    std::clog << '.';
    growNetwork(nVertices, nEdges, degrees, skewness);

    // Relabel node indices after sorting with respect to degree
    std::clog << '.';
    sortNetwork(nVertices, degrees);

    std::clog << '.';
}

void GeneticArchitecture::initializeNetwork(size_t &nEdges, std::vector<size_t> &degrees)
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

//----------------------------------------

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


/// Function to get a vector of chromosome lengths
std::vector<double> GeneticArchitecture::getChromosomeSizes(const size_t &nChromosomes)
{
    std::vector<double> chromosomeSizes;

    // Chromosomes all have the same size
    for (size_t i = 0u; i < nChromosomes; ++i) {
        chromosomeSizes.push_back((i + 1.0) / nChromosomes);
    }

    return chromosomeSizes;
}


/// Function to generate a vector of gene locations across the genome
std::vector<double> GeneticArchitecture::getGenomicLocations(const size_t &nLoci)
{
    std::vector<double> genomicLocations;

    // For every gene sample a uniform location
    for (size_t locus = 0u; locus < nLoci; ++locus)
        genomicLocations.push_back(rnd::uniform());

    // Sort gene locations by increasing order
    std::sort(genomicLocations.begin(), genomicLocations.end());

    return genomicLocations;
}


void GeneticArchitecture::sampleGeneLocations(const ParameterSet &parameters)
{

    // Declare a vector of gene locations
    std::vector<double> geneLocations;

    // For every gene sample a uniform location
    for (size_t i = 0u; i < parameters.nLoci; ++i) {
        geneLocations.push_back(rnd::uniform());
    }

    // Sort gene locations by increasing order
    std::sort(geneLocations.begin(), geneLocations.end());

    // Assign to each gene its chromosome and its location
    for (size_t i = 0u, lg = 0u; i < parameters.nLoci; ++i) {

        while (lg < parameters.nChromosomes - 1u && geneLocations[i] > chromosomeSizes[lg]) {
            ++lg;
        }

        LocusConstants gene;
        gene.chromosome = lg;
        gene.location = geneLocations[i];
        locusConstants.push_back(gene);

    }

}


/// Function to get a vector of effect sizes for all loci across the genome
std::vector<double> GeneticArchitecture::getEffectSizes(const size_t &nLoci, const double &shape, const double &scale,
        const size_t &nTraits)
{
    std::vector<double> effectSizes;
    std::vector<double> sumsqEffectSizes {0.0, 0.0, 0.0};

    // Loop through loci in the genome
    for (size_t locus = 0u; locus < nLoci; ++locus) {

        // For each one sample an effect size from a bilateral Gamma distribution
        double effectSize = std::gamma_distribution<double>(shape, scale)(rnd::rng);
        effectSize = rnd::bernoulli(0.5) ? effectSize * -1.0 : effectSize;
        effectSizes.push_back(effectSize);

        // Accumulate the sums of squared effect sizes for each trait separately
        sumsqEffectSizes[locusEncodedTraits[locus]] += sqr(effectSize);
    }

    // Take the square root of the sum of squares in order to normalize
    for (size_t trait = 0u; trait < nTraits; ++trait)
        sumsqEffectSizes[trait] > 0.0 ? sqrt(sumsqEffectSizes) : 1.0;

    // Normalize all locus effect sizes by the square rooted sum of squares of their respective trait
    for (size_t locus = 0u; locus < nLoci; ++locus)
        effectSizes[locus] /= sumsqEffectSizes[locusEncodedTraits[locus]];

    return effectSizes;
}


/// Function to get a vector of dominance coefficients across the genome
std::vector<double> GeneticArchitecture::getDominanceCoeffs(const size_t &nLoci, const size_t &nTraits)
{

    std::vector<double> dominanceCoeffs;
    std::vector<double> sumsqDominanceCoeffs {0.0, 0.0, 0.0};

    // Loop through loci across the genome
    for (size_t locus = 0u; locus < nLoci; ++locus) {

        // For each one sample a dominance coefficient from a unilateral normal distribution
        const double dominanceCoeff = fabs(rnd::normal(0.0, 1.0));
        dominanceCoeffs.push_back(dominanceCoeff);

        // Accumulate the sums of squared dominance coefficients for each trait
        sumsqDominanceCoeffs[locusEncodedTraits[locus]] += sqr(dominanceCoeff);
    }

    // Take the square root of the sum of squares in order to normalize
    for (size_t trait = 0u; trait < nTraits; ++trait)
        sumsqDominanceCoeffs[trait] > 0.0 ? sqrt(sumsqDominanceCoeffs) : 1.0;

    // Normalize the dominance coefficient of each locus by the square rooted sum of squares for its respective trait
    for (size_t locus = 0u; locus < nLoci; ++locus)
        dominanceCoeffs[locus] /= sumsqDominanceCoeffs[locusEncodedTraits[locus]];

    return dominanceCoeffs;
}


/// Function to get a vector of traits encoded by each locus across the genome
std::vector<size_t> GeneticArchitecture::getEncodedTraits(const size_t &nTraits, const std::vector<size_t> &nVertices)
{
    std::vector<size_t> encodedTraits;

    // Create an ordered vector of trait indices, one for each locus
    for (size_t trait = 0u; trait < nTraits; ++trait) {
        for (size_t vertex = 0u; vertex < nVertices[trait]; ++vertex) {
            encodedTraits.push_back(trait);
        }
    }

    // Shuffle the vector of encoded traits
    std::shuffle(encodedTraits.begin(), encodedTraits.end(), rnd::rng);

    return encodedTraits;
}


void GeneticArchitecture::assignPhenotypes(const ParameterSet &parameters)
{

    // Create a vector with the number of genes underlying each trait
    std::vector<size_t> nVertices {parameters.nEcoLoci, parameters.nMatLoci, parameters.nNtrLoci};

    // Randomize locus indices to make sure traits are assigned at random
    for (size_t i = 0u; i < parameters.nLoci; ++i) {
        loci.push_back(i);
    }
    std::shuffle(loci.begin(), loci.end(), rnd::rng);

    // For each phenotypic trait
    for (size_t trait = 0u, locus = 0u; trait < parameters.nTraits; ++trait) {

        std::vector<size_t> networkVerticesCurrTrait;

        // Assign loci to the trait
        for (size_t vertex = 0u; vertex < nVertices[trait]; ++vertex, ++locus) {

            networkVerticesCurrTrait.push_back(loci[locus]);
            locusConstants[loci[locus]].trait = trait;
        }

        networkVertices.push_back(networkVerticesCurrTrait);

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
        throw std::runtime_error("Unable to open genetic architecture file");
    }

    // Make sure the loaded architecture is
    std::clog << "  Validation.";
    bool isValidArchitecture = validateArchitecture(ifs, parameters);
    if (!isValidArchitecture) {
        throw std::logic_error("Genetic architecture in file is incompatible with current parameters");
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
        << parameters.nNtrInteractions << '\n'
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


*/























