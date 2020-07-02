#ifndef EXPLICITGENOMESPECIATION_PRINTER_H
#define EXPLICITGENOMESPECIATION_PRINTER_H

// This module is used to save data to output files, either from a MetaPop
// (for raw data such as full genomes) or a Collector (for summary statistics),
// or both

#include "Utilities.h"
#include "Collector.h"
#include "MetaPop.h"
#include <cassert>

class Printer
{

public:

    Printer(const std::string& = "");
    ~Printer();

    void print(const size_t&, const Collector&, const MetaPop&);
    void freeze(const MetaPop&, const Param&);
    void shutdown();

private:

    std::vector<std::string> whattosave(const std::string&) const;

    std::vector<std::string> filenames;
    vecStreams files;
    Stream freezer;

};


#endif
