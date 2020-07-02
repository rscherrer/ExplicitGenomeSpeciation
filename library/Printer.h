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

    Printer(const std::string &orderfile = "") :
        filenames(whattosave(orderfile)),
        files({ }),
        freezer(new std::ofstream)
    {

        files.reserve(filenames.size());

        // Open files
        for (size_t f = 0u; f < filenames.size(); ++f) {

            const std::string filename = filenames[f] + ".dat";
            Stream out(new std::ofstream);
            out->open(filename.c_str(), std::ios::binary);
            if (!out->is_open()) {
                std::string msg = "Unable to open output file " + filename;
                throw std::runtime_error(msg);
            }
            files.push_back(out);
        }

        // Open the freezer
        const std::string freezername = "freezer.dat";
        freezer->open(freezername.c_str(), std::ios::binary);
        if (!freezer->is_open()) {
            std::string msg = "Unable to open output freezer file";
            throw std::runtime_error(msg);
        }

    }

    ~Printer()
    {
        shutdown(); // close files
    }

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
