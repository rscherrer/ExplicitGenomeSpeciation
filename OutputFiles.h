//
// Created by p278834 on 15-5-2019.
//

#ifndef EXPLICITGENOMESPECIATION_OUTPUTFILES_H
#define EXPLICITGENOMESPECIATION_OUTPUTFILES_H

#include <fstream>
#include "ParameterSet.h"


class OutputFiles {

public:

    std::ofstream logFile;
    std::ofstream datFile;
    std::ofstream arcFile;

    void openAll(const size_t &seed);
    void writeParameters(const ParameterSet &parameters);
    template <class T>
    void writeLogLine(std::string&, T&);
    void writeLogLineVector(std::string&, std::vector<double>&);

};


#endif //EXPLICITGENOMESPECIATION_OUTPUTFILES_H
