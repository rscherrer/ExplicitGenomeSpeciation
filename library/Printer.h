#ifndef EXPLICITGENOMESPECIATION_PRINTER_H
#define EXPLICITGENOMESPECIATION_PRINTER_H

// This module is used to save summary statistics from the Collector to output
// files

#include "Utilities.h"
#include "Collector.h"
#include "MetaPop.h"
#include <cassert>

class Printer
{

public:

    Printer(const std::string& = "", const bool& = true);
    ~Printer();

    void print(const size_t&, const Collector&, const MetaPop&);
    void shutdown();

private:

    std::vector<std::string> whattosave(const std::string&) const;

    std::vector<std::string> filenames;
    vecStreams files;

};


#endif
