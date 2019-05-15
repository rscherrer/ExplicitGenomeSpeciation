//
// Created by p278834 on 7-5-2019.
//

#include "ParameterSet.h"
#include <iostream>
#include <fstream>
#include "random.h"

ParameterSet::ParameterSet()
{

}

ParameterSet::ParameterSet(const std::string &parameterFileName)
{

}



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

// Getters

bool ParameterSet::getIsGenerateArchitecture() const
{
    return isGenerateArchitecture;
}

size_t ParameterSet::getSeed() const
{
    return seed;
}

std::string ParameterSet::getArchitectureFilename() const
{
    return architectureFilename;
}

// Setters

void ParameterSet::setDefaultSeed()
{
    seed = rnd::set_seed();
}

void ParameterSet::setIsGenerateArchitecture(const bool &yesOrNo)
{
    isGenerateArchitecture = yesOrNo;
}

void ParameterSet::readParameters(const std::string& filename)
{
    std::clog << "reading parameters from file " << filename << '\n';

    //open parameter file
    std::ifstream ifs(filename);
    if(!ifs.is_open())
        throw std::runtime_error("unable to open parameter file in readParameters()");

    std::string str;

    ifs >> str;
    if(str == "rng_seed_clock") seed = rnd::set_seed();
    else if(str == "rng_seed_user") {
        ifs >> seed;
        rnd::set_seed(seed);
    }
    else throw std::logic_error("\'rng_seed_clock\' or \'rng_seed_user <arg>\' expected at first line of parameterfile\n");

    ifs >> str;
    if(str == "architecture_generate") generateArchitecture = true;
    else if(str == "architecture_load") {
        generateArchitecture = false;
        ifs >> architecture;
    }
    else throw std::logic_error("\'architecture_generate\' or \'architecture_load <arg>\' expected at second line of parameterfile\n");

    while(ifs >> str) {
        if(read(str, "initial_sequence", sequenceString, ifs))
        {
            for(auto a : sequenceString)
            {
                sequence.push_back(a == '1');
            }
        }
        else if(str == "scale_A")
            for(size_t crctr = 0u; crctr < nCharacter; ++crctr) {
                ifs >> scaleA[crctr];
                std::clog   << "parameter " << "scale_A[" << crctr
                            << "] set to " << scaleA[crctr] << '\n';
            }
        else if(str == "scale_D")
            for(size_t crctr = 0u; crctr < nCharacter; ++crctr) {
                ifs >> scaleD[crctr];
                std::clog   << "parameter " << "scale_D[" << crctr
                            << "] set to " << scaleD[crctr] << '\n';
            }
        else if(str == "scale_I")
            for(size_t crctr = 0u; crctr < nCharacter; ++crctr) {
                ifs >> scaleI[crctr];
                std::clog   << "parameter " << "scale_I[" << crctr
                            << "] set to " << scaleI[crctr] << '\n';
            }
        else if(str == "scale_E")
            for(size_t crctr = 0u; crctr < nCharacter; ++crctr) {
                ifs >> scaleE[crctr];
                std::clog   << "parameter " << "scale_E[" << crctr
                            << "] set to " << scaleE[crctr] << '\n';
            }
        else if(read(str, "mutation_rate", mutationRate, ifs));
        else if(read(str, "genome_size_cm", mapLength, ifs));
        else if(read(str, "female_heterogamety", isFemaleHeteroGamety, ifs));
        else if(read(str, "dispersal_rate", dispersalRate, ifs));
        else if(read(str, "alpha", alpha, ifs));
        else if(read(str, "beta", beta, ifs));
        else if(read(str, "prob_survival", survivalProb, ifs));
        else if(read(str, "habitat_asymmetry", habitatAsymmetry, ifs));
        else if(read(str, "sel_coeff_ecol", ecoSelCoeff, ifs));
        else if(read(str, "preference_strength", matePreferenceStrength, ifs));
        else if(read(str, "preference_cost", mateEvaluationCost, ifs));
        else if(read(str, "incompatibility_cost", costIncompat, ifs));
        else if(read(str, "typeII_resource_utilisation", isTypeIIResourceUtilisation, ifs));
        else if(read(str, "typeII_mate_choice", isTypeIIMateChoice, ifs));
        else if(read(str, "network_skewness", networkSkewness, ifs));
        else if(str == "t_end") {
            ifs >> tBurnIn >> tEndSim;
            std::clog   << "burn-in period  " << tBurnIn << " generations \n";
            std::clog   << "simulation time " << tEndSim << " generations \n";
        }
        else if(str == "t_dat") {
            ifs >> tGetDat >> tSavDat;
            std::clog   << "data collected every " << tGetDat << " generations \n";
            std::clog   << "data stored    every " << tSavDat << " generations \n";
        }
        else throw std::runtime_error("unknown parameter " + str);
    }
    std::clog << "parameters were read in successfully\n";
}

void ParameterSet::setArchitectureFilename(const std::string &filename)
{
    architectureFilename = filename;
}

void ParameterSet::writeParameters(std::ofstream &ofs, const char sep = ' ')
{
    ofs << "rng_seed"       << sep  << seed         << '\n';
    ofs << "architecture"   << sep  << architecture << '\n';
    ofs << "scale_A";
    for(size_t crctr = 0u; crctr < nCharacter; ++crctr)
        ofs << sep << scaleA[crctr];
    ofs << '\n';
    ofs << "scale_D";
    for(size_t crctr = 0u; crctr < nCharacter; ++crctr)
        ofs << sep << scaleD[crctr];
    ofs << '\n';
    ofs << "scale_I";
    for(size_t crctr = 0u; crctr < nCharacter; ++crctr)
        ofs << sep << scaleI[crctr];
    ofs << '\n';
    ofs << "scale_E";
    for(size_t crctr = 0u; crctr < nCharacter; ++crctr)
        ofs << sep << scaleE[crctr];
    ofs << '\n';
    ofs << "mutation_rate"  << sep  << mutationRate << '\n';
    ofs << "genome_size_cm" << sep  << mapLength    << '\n';
    ofs << "female_heterogamety" << sep  << isFemaleHeteroGamety << '\n';
    ofs << "dispersal_rate" << sep  << dispersalRate << '\n';
    ofs << "alpha" << sep  << alpha << '\n';
    ofs << "beta" << sep  << beta << '\n';
    ofs << "prob_survival" << sep  << survivalProb << '\n';
    ofs << "habitat_asymmetry" << sep  << habitatAsymmetry << '\n';
    ofs << "sel_coeff_ecol" << sep  << ecoSelCoeff << '\n';
    ofs << "preference_strength" << sep  << matePreferenceStrength << '\n';
    ofs << "preference_cost" << sep  << mateEvaluationCost << '\n';
    ofs << "incompatibility_cost" << sep << costIncompat << '\n';
    ofs << "typeII_resource_utilisation" << sep << isTypeIIResourceUtilisation << '\n';
    ofs << "typeII_mate_choice" << sep << isTypeIIMateChoice << '\n';
    ofs << "network_skewness" << sep << networkSkewness << '\n';
    ofs << "t_end" << sep  << tBurnIn << sep << tEndSim << '\n';
    ofs << "t_dat" << sep  << tGetDat << sep << tSavDat << '\n';
    ofs << "initial_sequence" << sep << (sequence.size() == nBits ? to_string(sequence) : "random") << '\n';
}