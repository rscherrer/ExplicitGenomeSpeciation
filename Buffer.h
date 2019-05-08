//
// Created by p278834 on 8-5-2019.
//

#ifndef EXPLICITGENOMESPECIATION_BUFFER_H
#define EXPLICITGENOMESPECIATION_BUFFER_H

#include <string>
#include <vector>
#include <fstream>
#include "ParameterSet.h"
#include "Population.h"
#include "Genome.h"

class Buffer {

public:
    Buffer(const std::string&, const ParameterSet&, const Population&, const Genome&);


    ~Buffer() { ofs.close(); }
    double &operator[](size_t j) {return data[i][j];}
    void flush(const ParameterSet&);

private:
    size_t i, k;
    int t;
    const size_t n;
    const std::string label;
    std::ofstream ofs;
    std::vector< std::vector<double> > data;
    static const char sep = ',';
};



#endif //EXPLICITGENOMESPECIATION_BUFFER_H
