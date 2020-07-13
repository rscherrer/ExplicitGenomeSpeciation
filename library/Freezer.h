#ifndef EXPLICITGENOMESPECIATION_FREEZER_H
#define EXPLICITGENOMESPECIATION_FREEZER_H

// This module is a special printer used for saving individual full genomes,
// directly from the MetaPop

#include "Utilities.h"
#include "MetaPop.h"
#include "Printer.h"
#include <cassert>

// Note: Using 64bit integers is a way to save space when saving full genomes
// The stf::write function saves in binary format, meaning that it will
// convert back the integers into their underlying bitsets
// The resulting binary output file should therefore be interpreted as
// a bit-wise array; each value (allele) is encoded by a single bit


class Freezer
{

public:

    Freezer(const std::string& = "freezer.dat",
     const std::string& = "locivalues.dat", const bool& = true);
    ~Freezer();

    void freeze(const MetaPop&, const size_t&);
    void saveIndivGenome(const Individual&, const size_t&);
    void shutdown();

private:

    Stream freezer;
    Stream locivalues;

};

#endif // FREEZER_H
