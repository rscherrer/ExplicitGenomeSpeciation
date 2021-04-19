#ifndef EXPLICITGENOMESPECIATION_PEDIGREE_H
#define EXPLICITGENOMESPECIATION_PEDIGREE_H

#include "Utilities.h"
#include "MetaPop.h"
#include "Printer.h"
#include <cassert>

class Pedigree
{

public:

    Pedigree(const std::string&, const bool&);
    ~Pedigree();

    void analyze(const MetaPop&, const Param&, const GenArch&);
    void shutdown();

private:

    Stream pedigree;

};

#endif // PEDIGREE_H
