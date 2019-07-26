#include "GeneticArchitecture.h"
#include "ParameterSet.h"
#include "Random.h"
#include "utils.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>


/// A network edge is a pair of locus indices
typedef std::pair<size_t, size_t> Edge;


/// Constructor of genetic architecture
GeneticArchitecture::GeneticArchitecture(const ParameterSet &pars,
 Random &rnd) :
    nTraits(pars.getNTraits()),
    nChromosomes(pars.getNChromosomes()),
    nLoci(pars.getNLoci()),
    nLociPerTrait(pars.getNLociPerTrait()),
    nEdgesPerTrait(pars.getNEdgesPerTrait()),
    skewnesses(pars.getSkewnesses()),
    effectSizeShape(pars.getEffectSizeShape()),
    effectSizeScale(pars.getEffectSizeScale()),
    interactionWeightShape(pars.getInteractionWeightShape()),
    interactionWeightScale(pars.getInteractionWeightScale()),
    chromosomeSizes(makeChromosomeSizes()),
    traitNetworks(makeTraitNetworks(rnd)),
    genome(makeGenome(rnd))
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
std::vector<Network> GeneticArchitecture::makeTraitNetworks(Random &rnd) const
 noexcept
{
    std::vector<Network> networks;

    // For each trait
    for (size_t trait = 0u; trait < nTraits; ++trait)
    {

        // Make a network map (a vector of edges) for the current trait using
        // the preferential
        // attachment algorithm
        Network network = Network(nLociPerTrait[trait], nEdgesPerTrait[trait],
         skewnesses[trait], interactionWeightShape, interactionWeightShape,
          rnd);

        networks.push_back(network);
    }

    // The indices in these network maps are indices among the loci underlying
    // a given trait,
    // not absolute loci indices across the genome

    assert(networks.size() == nTraits);
    for (size_t trait = 0u; trait < nTraits; ++trait)
        assert(networks[trait].map.size() == nEdgesPerTrait[trait]);

    return networks;
}


/// Network constructor
Network::Network(const size_t &nvertices, const size_t &nedges,
 const double &skew, const double &shape, const double &scale, Random &rnd) :
        nVertices(nvertices),
        nEdges(nedges),
        skewness(skew),
        map(makeNetwork(nedges, rnd)),
        weights(makeWeights(shape, scale, rnd))
{}


/// Function to make a new interaction network based on the preferential
/// attachment algorithm
std::vector<Edge> Network::makeNetwork(size_t nedges, Random &rnd) const
 noexcept
{
    std::vector<Edge> network;

    assert(nVertices > 1u);

    if (nedges == 0u) return network;

    // Start the network
    std::vector<size_t> degrees(nVertices, 0u);
    initializeNetwork(network, nedges, degrees);

    // Grow network by linking preferentially to well-connected nodes
    growNetwork(network, nedges, degrees, rnd);

    // Relabel node indices after sorting with respect to degree
    sortNetwork(network, degrees);

    return network;
}


/// Function to create a new network
void Network::initializeNetwork(std::vector<Edge> &network, size_t &nedges,
                                std::vector<size_t> &degrees) const noexcept
{

    if (!network.empty()) network.clear();

    // Create initial network
    network.emplace_back(Edge {0u, 1u});
    degrees[0u] = degrees[1u] = 1u;
    --nedges;

}


/// Function to iteratively grow a network based on the preferential attachment
/// algorithm
void Network::growNetwork(std::vector<Edge> &network, size_t &nedges,
 std::vector<size_t> &degrees, Random &rnd) const noexcept
{

    // For each vertex in the network...
    for (size_t vertex = 2u; vertex < nVertices && nedges > 0u; ++vertex) {

        // Assign probability weights to all potential partners
        std::vector<double> probs(vertex);

        for (size_t partner = 0u; partner < vertex; ++partner) {
            probs[partner] = pow(degrees[partner], skewness);
        }

        // Sample the number of attachments of the current vertex
        size_t nAttachments =
                (vertex == nVertices - 1u ?
                     nedges : rnd.binomial(nedges, 1.0 / (nVertices - vertex)));

        // For each new attachment...
        while (nAttachments) {

            // Compute the sum of weights, quit the loop if too small
            double sumProbs = 0.0;
            for (size_t partner = 0u; partner < vertex; ++partner)
                sumProbs += probs[partner];
            if (sumProbs < 1.0) break;

            // Sample the new partner without replacement
            std::discrete_distribution<size_t> attachmentProbs(probs.begin(),
             probs.end());
            size_t partner = attachmentProbs(rnd.rng);

            // Add the new partner to the network
            network.emplace_back(Edge {vertex, partner});
            probs[partner] = 0.0;
            ++degrees[partner];
            ++degrees[vertex];
            --nAttachments;
            --nedges;
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
void Network::sortNetwork(std::vector<Edge> &network,
 const std::vector<size_t> &degrees) const noexcept
{

    // Compute the ranks of all vertices with respect to their degrees
    std::vector<size_t> ranks(nVertices, 0u);
    for (size_t vertex = 0u; vertex < nVertices - 1u; ++vertex) {
        for (size_t othervertex = vertex + 1u; othervertex < nVertices;
             ++othervertex)
        {
            if (degrees[othervertex] > degrees[vertex])
            {
                ++ranks[vertex];
            }
            else
            {
                ++ranks[othervertex];
            }
        }
    }

    // For each edge, swap partners to make sure the first partner always have
    // the lower rank
    for (Edge &edge : network) {
        edge.first = ranks[edge.first];
        edge.second = ranks[edge.second];
        if (edge.first > edge.second)
        {
            std::swap(edge.first, edge.second);
        }
    }

    // Sort the edges
    std::sort(network.begin(), network.end(), edgeCompare);

}


/// Function to sample interaction weights across edges of a gene regulatory
/// network
std::vector<double> Network::makeWeights(const double &shape,
 const double &scale, Random &rnd) const noexcept
{
    std::vector<double> interweights;
    double sqrtsumsqWeights = 0.0;

    // For each edge in the network...
    for (size_t edge = 0u; edge < nEdges; ++edge) {

        // Sample the weight from a two-sided Gamma distribution
        double weight = std::gamma_distribution<double>(shape, scale)(rnd.rng);
        weight = rnd.bernoulli(0.5) ? weight * -1.0 : weight;
        interweights.push_back(weight);

        // Accumulate square rooted sum of squared interaction weights for
        // later normalization
        sqrtsumsqWeights += sqr(weight);
    }

    // Square root the normalizing factor
    sqrtsumsqWeights = sqrtsumsqWeights > 0.0 ? sqrt(sqrtsumsqWeights) : 1.0;

    // Normalize at the end
    for (size_t edge = 0u; edge < nEdges; ++edge)
        interweights[edge] /= sqrtsumsqWeights;

    return interweights;
}


/// Function from architecture to call the Genome constructor
Genome GeneticArchitecture::makeGenome(Random &rnd) const noexcept
{
    const Genome gen = Genome(nTraits, nLociPerTrait, nLoci, effectSizeShape,
     effectSizeScale, rnd);
    return gen;
}


/// Genome constructor
Genome::Genome(const size_t &nTraits, const std::vector<size_t> &nLociPerTrait,
 const size_t &nLoci, const double &shape, const double &scale, Random &rnd) :
    encodedTraits(makeEncodedTraits(nTraits, nLociPerTrait, rnd)),
    locations(std::vector<double> { 0.0 }),
    effectSizes(std::vector<double> { 0.0 }),
    dominanceCoeffs(std::vector<double> { 0.0 })
{

    locations.pop_back();
    effectSizes.pop_back();
    dominanceCoeffs.pop_back();

    // Sample locations, effect sizes and dominance coefficients across the
    // genome
    setLocationsEffectSizesAndDominance(nTraits, nLoci, shape, scale, rnd);

    assert(encodedTraits.size() == nLoci);
    assert(effectSizes.size() == nLoci);
    assert(dominanceCoeffs.size() == nLoci);
    assert(locations.size() == nLoci);
}


/// Function to randomly assign loci to their encoded traits
std::vector<size_t> Genome::makeEncodedTraits(const size_t &nTraits,
 const std::vector<size_t> &nLociPerTrait, Random &rnd) const noexcept
{

    std::vector<size_t> traits;

    // Make an ordered vector of trait indices
    for (size_t trait = 0u; trait < nTraits; ++trait) {
        for (size_t locus = 0u; locus < nLociPerTrait[trait]; ++locus) {
            traits.push_back(trait);
            assert(traits.back() <= 2u);
        }
    }

    // Shuffle encoded traits randomly
    std::shuffle(traits.begin(), traits.end(), rnd.rng);

    return traits;

}


/// Function to sample locations, effect sizes and dominance across the genome
void Genome::setLocationsEffectSizesAndDominance(const size_t &nTraits,
 const size_t &nLoci, const double &shape, const double &scale, Random &rnd)
{

    // Prepare squared roots of sums of squared effect sizes and dominance
    // coefficients
    std::vector<double> sqrtsumsqEffectSizes {0.0, 0.0, 0.0};
    std::vector<double> sqrtsumsqDominanceCoeffs {0.0, 0.0, 0.0};

    // Sample locations, effect sizes and dominance coefficients
    for (size_t locus = 0u; locus < nLoci; ++locus)
    {

        // Locations are sampled uniformly
        locations.push_back(rnd.uniform(1.0));

        assert(locations.back() > 0.0);
        assert(locations.back() < 1.0);

        // Effect sizes are sampled from a two-sided Gamma distribution
        double effectsize = std::gamma_distribution<double>(shape, scale)(
                    rnd.rng);
        effectsize = rnd.bernoulli(0.5) ? effectsize * -1.0 : effectsize;
        effectSizes.push_back(effectsize);

        // Squared effect sizes are accumulated for normalizing
        sqrtsumsqEffectSizes[encodedTraits[locus]] += sqr(effectsize);

        // Dominance coefficients are sampled from a one-sided normal
        // distribution
        const double dominance = fabs(rnd.normal(0.0, 1.0));
        assert(dominance > 0.0);
        dominanceCoeffs.push_back(dominance);

        // Squared dominance coefficients are accumulated for normalizing
        sqrtsumsqDominanceCoeffs[encodedTraits[locus]] += sqr(dominance);

    }

    // Sort gene locations by increasing order
    std::sort(locations.begin(), locations.end());

    for (size_t locus = 1u; locus < nLoci; ++locus)
        assert(locations[locus] > locations[locus - 1u]);

    // Take the square roots of the sums of squares for normalization
    for (size_t trait = 0u; trait < nTraits; ++trait)
    {
        sqrtsumsqEffectSizes[trait] = sqrt(sqrtsumsqEffectSizes[trait]);
        sqrtsumsqDominanceCoeffs[trait] = sqrt(sqrtsumsqDominanceCoeffs[trait]);

        assert(sqrtsumsqEffectSizes[trait] > 0.0);
        assert(sqrtsumsqDominanceCoeffs[trait] > 0.0);
    }

    // Normalize effect sizes and dominance coefficients across the genome
    for (size_t locus = 0u; locus < nLoci; ++locus)
    {
        effectSizes[locus] /= sqrtsumsqEffectSizes[encodedTraits[locus]];
        dominanceCoeffs[locus] /= sqrtsumsqDominanceCoeffs[
                encodedTraits[locus]];
    }
}
