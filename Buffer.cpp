//
// Created by p278834 on 8-5-2019.
//

#include "Buffer.h"
#include <sstream>

Buffer::Buffer(const std::string &str,
               const ParameterSet& parameters,
               const Population& population,
               const Genome& genome) :
        i(0u), k(0u), t(-parameters.tBurnIn), n(static_cast<size_t>(parameters.tSavDat / parameters.tGetDat)), label(str)
{
    data = std::vector< std::vector<double> >(n, std::vector<double>(parameters.nLoci, 0.0));

    // open table data file
    std::ostringstream oss;
    oss << "simulation_" << parameters.seed << "_table_" << label << ".csv";
    ofs.open(oss.str());
    if(!ofs.is_open())
        throw std::runtime_error("unable to open output file in Buffer::Buffer()");

    // write header
    ofs << sep << "generation";
    for(size_t j = 0u; j < parameters.nLoci; ++j)
        ofs << sep << "loc." << j;
    ofs << '\n';
    ofs << "character" << sep << "NA";
    for(size_t j = 0u; j < parameters.nLoci; ++j)
        ofs << sep << genome.characterLocus[j].character;
    ofs << '\n';
    ofs << "linkage.group" << sep << "NA";
    for(size_t j = 0u; j < parameters.nLoci; ++j)
        ofs << sep << genome.characterLocus[j].linkageGroup;
    ofs << '\n';
    ofs << "degree" << sep << "NA";
    for(size_t j = 0u; j < parameters.nLoci; ++j)
        ofs << sep << genome.characterLocus[j].edges.size();
    ofs << '\n';
    ofs << "location" << sep << "NA";
    for(size_t j = 0u; j < parameters.nLoci; ++j)
        ofs << sep << genome.characterLocus[j].location;
    ofs << '\n';
    ofs << "effect.size" << sep << "NA";
    for(size_t j = 0u; j < parameters.nLoci; ++j)
        ofs << sep << genome.characterLocus[j].effectSize;
    ofs << '\n';
    ofs << "dominance.coeff" << sep << "NA";
    for(size_t j = 0u; j < parameters.nLoci; ++j)
        ofs << sep << genome.characterLocus[j].dominanceCoeff;
    ofs << '\n';
}

void Buffer::flush(const ParameterSet& parameters)
{
    t += parameters.tGetDat;
    ++i;
    if(t % parameters.tSavDat == 0u) {
        ++k;
        ofs << "point." << k << sep << t - (parameters.tSavDat >> 1u);
        for(size_t j = 0u; j < parameters.nLoci; ++j) {
            double sum = 0.0;
            for(i = 0u; i < n; ++i) sum += data[i][j];
            ofs << sep << sum / n;
        }
        ofs << '\n';
        i = 0u;
    }
    ofs.flush();
}