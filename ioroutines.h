//
// Created by p278834 on 9-5-2019.
//

#ifndef EXPLICITGENOMESPECIATION_IOROUTINES_H
#define EXPLICITGENOMESPECIATION_IOROUTINES_H

#include <string>
#include <iostream>
#include <fstream>
#include "ParameterSet.h"
#include "random.h"

template <class T>
bool read(const std::string &str, const std::string &name, T &par, std::ifstream &ifs)
{
    if(str == name) {
        ifs >> par;
        std::clog << "parameter " << str << " set to " << par << '\n';
        return true;
    }
    else return false;
}

void readParameters(const std::string& filename, ParameterSet& parameters)
{
    std::clog << "reading parameters from file " << filename << '\n';

    //open parameter file
    std::ifstream ifs(filename);
    if(!ifs.is_open())
        throw std::runtime_error("unable to open parameter file in readParameters()");

    std::string str;

    ifs >> str;
    if(str == "rng_seed_clock") parameters.seed = rnd::set_seed();
    else if(str == "rng_seed_user") {
        ifs >> parameters.seed;
        rnd::set_seed(parameters.seed);
    }
    else throw std::logic_error("\'rng_seed_clock\' or \'rng_seed_user <arg>\' expected at first line of parameterfile\n");

    ifs >> str;
    if(str == "architecture_generate") parameters.generateArchitecture = true;
    else if(str == "architecture_load") {
        parameters.generateArchitecture = false;
        ifs >> parameters.architecture;
    }
    else throw std::logic_error("\'architecture_generate\' or \'architecture_load <arg>\' expected at second line of parameterfile\n");

    while(ifs >> str) {
        if(read(str, "initial_sequence", parameters.sequenceString, ifs))
        {
            for(auto a : parameters.sequenceString)
            {
                parameters.sequence.push_back(a == '1');
            }
        }
        else if(str == "scale_A")
            for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr) {
                ifs >> parameters.scaleA[crctr];
                std::clog   << "parameter " << "scale_A[" << crctr
                            << "] set to " << parameters.scaleA[crctr] << '\n';
            }
        else if(str == "scale_D")
            for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr) {
                ifs >> parameters.scaleD[crctr];
                std::clog   << "parameter " << "scale_D[" << crctr
                            << "] set to " << parameters.scaleD[crctr] << '\n';
            }
        else if(str == "scale_I")
            for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr) {
                ifs >> parameters.scaleI[crctr];
                std::clog   << "parameter " << "scale_I[" << crctr
                            << "] set to " << parameters.scaleI[crctr] << '\n';
            }
        else if(str == "scale_E")
            for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr) {
                ifs >> parameters.scaleE[crctr];
                std::clog   << "parameter " << "scale_E[" << crctr
                            << "] set to " << parameters.scaleE[crctr] << '\n';
            }
        else if(read(str, "mutation_rate", parameters.mutationRate, ifs));
        else if(read(str, "genome_size_cm", parameters.mapLength, ifs));
        else if(read(str, "female_heterogamety", parameters.isFemaleHeteroGamety, ifs));
        else if(read(str, "dispersal_rate", parameters.dispersalRate, ifs));
        else if(read(str, "alpha", parameters.alpha, ifs));
        else if(read(str, "beta", parameters.beta, ifs));
        else if(read(str, "prob_survival", parameters.survivalProb, ifs));
        else if(read(str, "habitat_asymmetry", parameters.habitatAsymmetry, ifs));
        else if(read(str, "sel_coeff_ecol", parameters.ecoSelCoeff, ifs));
        else if(read(str, "preference_strength", parameters.matePreferenceStrength, ifs));
        else if(read(str, "preference_cost", parameters.mateEvaluationCost, ifs));
        else if(read(str, "incompatibility_cost", parameters.costIncompat, ifs));
        else if(read(str, "typeII_resource_utilisation", parameters.isTypeIIResourceUtilisation, ifs));
        else if(read(str, "typeII_mate_choice", parameters.isTypeIIMateChoice, ifs));
        else if(read(str, "network_skewness", parameters.networkSkewness, ifs));
        else if(str == "t_end") {
            ifs >> parameters.tBurnIn >> parameters.tEndSim;
            std::clog   << "burn-in period  " << parameters.tBurnIn << " generations \n";
            std::clog   << "simulation time " << parameters.tEndSim << " generations \n";
        }
        else if(str == "t_dat") {
            ifs >> parameters.tGetDat >> parameters.tSavDat;
            std::clog   << "data collected every " << parameters.tGetDat << " generations \n";
            std::clog   << "data stored    every " << parameters.tSavDat << " generations \n";
        }
        else throw std::runtime_error("unknown parameter " + str);
    }
    std::clog << "parameters were read in successfully\n";
}

void writeParameters(std::ofstream &ofs, const ParameterSet& parameters, const char sep = ' ')
{
    ofs << "rng_seed"       << sep  << parameters.seed         << '\n';
    ofs << "architecture"   << sep  << parameters.architecture << '\n';
    ofs << "scale_A";
    for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr)
        ofs << sep << parameters.scaleA[crctr];
    ofs << '\n';
    ofs << "scale_D";
    for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr)
        ofs << sep << parameters.scaleD[crctr];
    ofs << '\n';
    ofs << "scale_I";
    for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr)
        ofs << sep << parameters.scaleI[crctr];
    ofs << '\n';
    ofs << "scale_E";
    for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr)
        ofs << sep << parameters.scaleE[crctr];
    ofs << '\n';
    ofs << "mutation_rate"  << sep  << parameters.mutationRate << '\n';
    ofs << "genome_size_cm" << sep  << parameters.mapLength    << '\n';
    ofs << "female_heterogamety" << sep  << parameters.isFemaleHeteroGamety << '\n';
    ofs << "dispersal_rate" << sep  << parameters.dispersalRate << '\n';
    ofs << "alpha" << sep  << parameters.alpha << '\n';
    ofs << "beta" << sep  << parameters.beta << '\n';
    ofs << "prob_survival" << sep  << parameters.survivalProb << '\n';
    ofs << "habitat_asymmetry" << sep  << parameters.habitatAsymmetry << '\n';
    ofs << "sel_coeff_ecol" << sep  << parameters.ecoSelCoeff << '\n';
    ofs << "preference_strength" << sep  << parameters.matePreferenceStrength << '\n';
    ofs << "preference_cost" << sep  << parameters.mateEvaluationCost << '\n';
    ofs << "incompatibility_cost" << sep << parameters.costIncompat << '\n';
    ofs << "typeII_resource_utilisation" << sep << parameters.isTypeIIResourceUtilisation << '\n';
    ofs << "typeII_mate_choice" << sep << parameters.isTypeIIMateChoice << '\n';
    ofs << "network_skewness" << sep << parameters.networkSkewness << '\n';
    ofs << "t_end" << sep  << parameters.tBurnIn << sep << parameters.tEndSim << '\n';
    ofs << "t_dat" << sep  << parameters.tGetDat << sep << parameters.tSavDat << '\n';
    ofs << "initial_sequence" << sep << (parameters.sequence.size() == parameters.nBits ? to_string(parameters.sequence) : "random") << '\n';
}

#endif //EXPLICITGENOMESPECIATION_IOROUTINES_H
