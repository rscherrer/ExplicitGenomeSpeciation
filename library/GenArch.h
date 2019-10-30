#ifndef EXPLICITGENOMESPECIATION_GENARCH_H
#define EXPLICITGENOMESPECIATION_GENARCH_H

#include "Param.h"
#include "Random.h"
#include "Network.h"
#include "Types.h"
#include "Utilities.h"
#include <cassert>
#include <stddef.h>
#include <cstddef>

typedef std::pair<size_t, size_t> Edge;
typedef std::vector<Network> MultiNet;
typedef std::vector<std::shared_ptr<std::ofstream> > vecStreams;

// The genetic architecture contains locus-specific details about the
// genotype-phenotype map. It is created using the parameters, but contains
// large arrays of values across the whole genome, and is therefore
// larger than the Param class.

class GenArch {

    friend class Collector;

public:

    GenArch(const Param &pars) :
        isreset(resetseed(pars.archseed)),
        chromosomes(makeChromosomes(pars)),
        traits(makeEncodedTraits(pars)),
        locations(makeLocations(pars)),
        effects(makeEffects(pars)),
        dominances(makeDominances(pars)),
        networks(makeNetworks(pars)),
        archfilenames(whattosave()),
        archfiles({ })
    {

        // If we want to use a specific genetic architecture
        // We can supply a seed just for the architecture
        // At the end we need to reset the seed back to the original

        rnd::rng.seed(pars.archseed);

        assert(utl::sumu(pars.nvertices) == pars.nloci);
        assert(chromosomes.size() == pars.nchrom);
        assert(traits.size() == pars.nloci);
        assert(effects.size() == pars.nloci);
        assert(dominances.size() == pars.nloci);
        assert(locations.size() == pars.nloci);
        assert(networks.size() == 3u);
        assert(isreset);

        // Save the architecture if necessary
        if (pars.archsave) {

            archfiles.reserve(archfilenames.size());

            // Open files
            for (size_t f = 0u; f < archfilenames.size(); ++f) {

                const std::string filename = archfilenames[f] + ".dat";
                std::shared_ptr<std::ofstream> out(new std::ofstream);
                out->open(filename.c_str(), std::ios::binary);
                if (!out->is_open()) {
                    std::string msg = "Unable to open output file " + filename;
                    throw std::runtime_error(msg);
                }
                archfiles.push_back(out);
            }

            // Write in files
            for (size_t c = 0u; c < pars.nchrom; ++c) {
                stf::write(chromosomes[c], archfiles[0u]);
            }

            size_t f = 1u; // file id

            size_t off; // offset to write multiple loci to the same file

            for (size_t l = 0u; l < pars.nloci; ++l) {

                off = 0u; // reset the offset

                stf::write(utl::size2dbl(traits[l]), archfiles[f + off]); ++off;
                stf::write(locations[l], archfiles[f + off]); ++off;
                stf::write(effects[l], archfiles[f + off]); ++off;
                stf::write(dominances[l], archfiles[f + off]); ++off;

            }

            f += off; // move on to network files

            for (size_t t = 0u; t < 3u; ++t) {
                for (size_t e = 0u; e < getNetworkSize(t); ++e) {

                    stf::write(utl::size2dbl(networks[t].edges[e].first), archfiles[f + off]); ++off;
                    stf::write(utl::size2dbl(networks[t].edges[e].second), archfiles[f + off]); ++off;
                    stf::write(networks[t].weights[e], archfiles[f + off]); ++off;

                }
            }

            // Close files
            f = 0u;
            for (; f < archfiles.size(); ++f) archfiles[f]->close();
        }

        rnd::rng.seed(pars.seed);

    }

    bool isreset;

    vecDbl chromosomes;     // per chromosome
    vecUns traits;          // per locus
    vecDbl locations;       // per locus
    vecDbl effects;         // per locus
    vecDbl dominances;      // per locus
    MultiNet networks;      // per trait

    vecStrings archfilenames;
    vecStreams archfiles;

    // Getters called from tests
    size_t getNetworkSize() const
    {
        size_t nedges = 0u;
        for (size_t trait = 0u; trait < 3u; ++trait) {
            nedges += getNetworkSize(trait);
        }
        return nedges;
    }
    size_t getNetworkSize(const size_t &trait) const
    {
        return networks[trait].map.size();
    }
    Edge getEdge(const size_t &trait, const size_t &edge) const
    {
        return networks[trait].map[edge];
    }
    size_t getSumTraits() const
    {
        size_t sum = 0u;
        for (auto x : traits) sum += x;
        return sum;
    }
    double getSumEffects() const
    {
        double sum = 0.0;
        for (auto x : effects) sum += x;
        return sum;
    }
    double getSumDominances() const
    {
        double sum = 0.0;
        for (auto x : dominances) sum += x;
        return sum;
    }
    double getSumWeights(const size_t &trait) const
    {
        double sum = 0.0;
        for (size_t edge = 0u; edge < networks[trait].weights.size(); ++edge)
            sum += networks[trait].weights[edge];
        return sum;
    }

private:

    bool resetseed(const size_t&) const;
    MultiNet makeNetworks(const Param&) const;
    vecDbl makeChromosomes(const Param&) const;
    vecUns makeEncodedTraits(const Param&) const;
    vecDbl makeLocations(const Param&) const;
    vecDbl makeEffects(const Param&) const;
    vecDbl makeDominances(const Param&) const;

    vecStrings whattosave() const;

};

#endif
