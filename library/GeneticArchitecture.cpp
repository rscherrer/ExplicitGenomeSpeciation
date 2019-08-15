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
GeneticArchitecture::GeneticArchitecture(const ParameterSet &pars) :
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
    chromosomeSizes(makeChromosomes()),
    networks(makeNetworks()),
    genome(makeGenome())
{}


/// Function to make a vector of chromosome sizes
std::vector<double> GeneticArchitecture::makeChromosomes() const noexcept
{

    std::vector<double> chromsizes;

    // Chromosomes all have the same size
    for (size_t i = 0u; i < nChromosomes; ++i)
        chromsizes.push_back((i + 1.0) / nChromosomes);

    return chromsizes;

}


/// Function to make a vector of interacting partner loci for each trait
std::vector<Network> GeneticArchitecture::makeNetworks() const
 noexcept
{
    std::vector<Network> nets;

    // For each trait
    for (size_t trait = 0u; trait < nTraits; ++trait)
    {

        // Make a network map (a vector of edges) for the current trait using
        // the preferential
        // attachment algorithm
        Network network = Network(trait, nLociPerTrait[trait],
         nEdgesPerTrait[trait], skewnesses[trait], interactionWeightShape,
          interactionWeightShape, genome);

        nets.push_back(network);
    }

    assert(nets.size() == nTraits);

    for (size_t trait = 0u; trait < nTraits; ++trait)
        assert(nets[trait].map.size() == nEdgesPerTrait[trait]);

    return nets;
}


/// Network constructor
Network::Network(const size_t &character, const size_t &nvertices,
 const size_t &nedges, const double &skew, const double &shape,
  const double &scale, const Genome &genome) :
        trait(character),
        nVertices(nvertices),
        nEdges(nedges),
        skewness(skew),
        map(makeMap(nedges)),
        loci(makeLoci(genome)),
        edges(makeEdges()),
        weights(makeWeights(shape, scale))
{}


/// Function to make a new interaction network based on the preferential
/// attachment algorithm
std::vector<Edge> Network::makeMap(size_t nedges) const
 noexcept
{
    std::vector<Edge> network;

    assert(nVertices > 1u);

    if (nedges == 0u) return network;

    // Start the network
    std::vector<size_t> degrees(nVertices, 0u);
    initializeNetwork(network, nedges, degrees);

    // Grow network by linking preferentially to well-connected nodes
    growNetwork(network, nedges, degrees);

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
 std::vector<size_t> &degrees) const noexcept
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
                     nedges : rnd::binomial(nedges, 1.0 / (nVertices - vertex)));

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
            size_t partner = attachmentProbs(rnd::rng);

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


/// Detect the loci underlying a trait's network
std::vector<size_t> Network::makeLoci(const Genome &genome) const noexcept
{

    std::vector<size_t> underlying;

    for (size_t locus = 0u; locus < genome.nloci; ++locus)
        if (genome.traits[locus] == trait)
            underlying.push_back(locus);

    return underlying;

}


/// Map the network actual edges in the genome
std::vector<Edge> Network::makeEdges()
{
    std::vector<Edge> mapped;

    for (size_t edge = 0u; edge < map.size(); ++edge) {
        Edge partners;
        partners.first = loci[map[edge].first];
        partners.second = loci[map[edge].second];
        mapped.push_back(partners);
    }

    return mapped;
}


/// Function to sample interaction weights across edges of a gene regulatory
/// network
std::vector<double> Network::makeWeights(const double &shape,
 const double &scale) const noexcept
{
    std::vector<double> interweights;
    double sqrtsumsqWeights = 0.0;

    // For each edge in the network...
    for (size_t edge = 0u; edge < nEdges; ++edge) {

        // Sample the weight from a two-sided Gamma distribution
        double weight = std::gamma_distribution<double>(shape, scale)(rnd::rng);
        weight = rnd::bernoulli(0.5) ? weight * -1.0 : weight;
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
Genome GeneticArchitecture::makeGenome() const noexcept
{
    const Genome gen = Genome(nTraits, nLociPerTrait, nLoci, effectSizeShape,
     effectSizeScale);
    return gen;
}


/// Genome constructor
Genome::Genome(const size_t &nTraits, const std::vector<size_t> &nLociPerTrait,
 const size_t &nLoci, const double &shape, const double &scale) :
    traits(makeEncodedTraits(nTraits, nLociPerTrait)),

    locations(std::vector<double> { 0.0 }),
    effects(std::vector<double> { 0.0 }),
    dominances(std::vector<double> { 0.0 })
{

    locations.pop_back();
    effects.pop_back();
    dominances.pop_back();


    // Sample locations, effect sizes and dominance coefficients across the
    // genome
    setLocationsEffectSizesAndDominance(nTraits, nLoci, shape, scale);

    assert(traits.size() == nLoci);
    assert(effects.size() == nLoci);
    assert(dominances.size() == nLoci);
    assert(locations.size() == nLoci);
}


/// Function to randomly assign loci to their encoded traits
std::vector<size_t> Genome::makeEncodedTraits(const size_t &nTraits,
 const std::vector<size_t> &nLociPerTrait) const noexcept
{


    std::vector<size_t> phenotypes;


    // Make an ordered vector of trait indices
    for (size_t trait = 0u; trait < nTraits; ++trait) {
        for (size_t locus = 0u; locus < nLociPerTrait[trait]; ++locus) {

            phenotypes.push_back(trait);
            assert(phenotypes.back() <= 2u);

        }
    }

    // Shuffle encoded traits randomly

    std::shuffle(phenotypes.begin(), phenotypes.end(), rnd::rng);

    return phenotypes;


}


/// Function to sample locations, effect sizes and dominance across the genome
void Genome::setLocationsEffectSizesAndDominance(const size_t &nTraits,
 const size_t &nLoci, const double &shape, const double &scale)
{

    // Prepare squared roots of sums of squared effect sizes and dominance
    // coefficients
    std::vector<double> sqrtsumsqEffectSizes {0.0, 0.0, 0.0};
    std::vector<double> sqrtsumsqDominanceCoeffs {0.0, 0.0, 0.0};

    // Sample locations, effect sizes and dominance coefficients
    for (size_t locus = 0u; locus < nLoci; ++locus)
    {

        // Locations are sampled uniformly
        locations.push_back(rnd::uniform(1.0));

        assert(locations.back() > 0.0);
        assert(locations.back() < 1.0);

        // Effect sizes are sampled from a two-sided Gamma distribution
        double effectsize = std::gamma_distribution<double>(shape, scale)(
         rnd::rng);
        effectsize = rnd::bernoulli(0.5) ? effectsize * -1.0 : effectsize;
        effects.push_back(effectsize);

        // Squared effect sizes are accumulated for normalizing
        sqrtsumsqEffectSizes[traits[locus]] += sqr(effectsize);

        // Dominance coefficients are sampled from a one-sided normal
        // distribution
        const double dominance = fabs(rnd::normal(0.0, 1.0));
        assert(dominance > 0.0);
        dominances.push_back(dominance);

        // Squared dominance coefficients are accumulated for normalizing
        sqrtsumsqDominanceCoeffs[traits[locus]] += sqr(dominance);

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
        effects[locus] /= sqrtsumsqEffectSizes[traits[locus]];

        dominances[locus] /= sqrtsumsqDominanceCoeffs[traits[locus]];

    }
}
