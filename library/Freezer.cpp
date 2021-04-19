#include "Freezer.h"

// Constructor
Freezer::Freezer(const std::string &freezername, const std::string &lociname,
 const bool &gensave) :
    freezer(new std::ofstream),
    locivalues(new std::ofstream)
{

    if (gensave) {

        // Open the freezer
        freezer->open(freezername, std::ios::binary);
        if (!freezer->is_open()) {
            std::string msg = "Unable to open output freezer file";
            throw std::runtime_error(msg);
        }

        locivalues->open(lociname, std::ios::binary);
        if (!locivalues->is_open()) {
            std::string msg = "Unable to open loci values file";
            throw std::runtime_error(msg);
        }

    }

}

// Destructor
Freezer::~Freezer()
{
    shutdown();
}

void Freezer::shutdown() {
    freezer->close();
    locivalues->close();
}

// Member functions
//-----------------

// Save the genomes of all individuals in the population
void Freezer::freeze(const MetaPop &pop, const size_t &nloci) {

    for (size_t i = 0u; i < pop.getSize(); ++i) {
        saveIndivGenome(pop.population[i], nloci);
        for (size_t locus = 0u; locus < nloci; ++locus) {
            stf::write(pop.getLocusValue(i, locus), locivalues);
        }
    }

}

// Save the full genome of an individual in 64bit-chunks
void Freezer::saveIndivGenome(const Individual &ind, const size_t &n) {

    // ind: individual
    // n: number of loci

    // Number of 64bit-chunks that can contain the full genome
    const size_t nchunks = n * 2u / 64u + 1u;

    // Save each chunk as a long integer
    for (size_t i = 0u; i < nchunks; ++i) {

        // Integer representation of a 64bit-chunk of genome
        const size_t chunk = ind.getGenomeChunk(i).to_ulong();
        stf::write2(chunk, freezer);

    }

}
