//
// Created by p278834 on 15-5-2019.
//

#ifndef EXPLICITGENOMESPECIATION_OUTPUTFILE_H
#define EXPLICITGENOMESPECIATION_OUTPUTFILE_H

#include "ParameterSet.h"
#include <fstream>

class OutputFile {

public:

    std::ofstream file;

    void open(const size_t &seed, const std::string &extension);

    // For logFile
    void writeParameters(const ParameterSet &parameters);
    template <class T>
    void writeLine(std::string&, T&);
    void writeLineVector(std::string&, std::vector<double>&);

    // For datFile
    void writeHeader();
    void addColumn(std::string &name);

};


#endif //EXPLICITGENOMESPECIATION_OUTPUTFILE_H
