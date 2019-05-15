//
// Created by p278834 on 7-5-2019.
//

#include "ParameterSet.h"
#include <iostream>
#include <fstream>
#include <sstream>
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



// Setters

void ParameterSet::setDefaultSeed()
{
    seed = rnd::set_seed();
}

void ParameterSet::readParameters(const std::string& filename)
{
    std::clog << "Reading parameters from file " << filename << '\n';

    // Open parameter file
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        throw std::runtime_error("Unable to open parameter file in readParameters()");
    }

    std::string input;

    inputFile >> input;
    if (input == "rng_seed_clock") {
        seed = rnd::set_seed();
    }
    else if (input == "rng_seed_user") {
        inputFile >> seed;
        rnd::set_seed(seed);
    }
    else {
        throw std::logic_error("\'rng_seed_clock\' or \'rng_seed_user <arg>\' expected at first line of parameterfile\n");
    }

    inputFile >> input;
    if (input == "architecture_generate") {
        isGenerateArchitecture = true;
    }
    else if (input == "architecture_load") {
        isGenerateArchitecture = false;
        inputFile >> architectureFileName;
    }
    else {
        throw std::logic_error("\'architecture_generate\' or \'architecture_load <arg>\' expected at second line of parameterfile\n");
    }

    while (inputFile >> input) {
        if (read(input, "initial_sequence", strProvidedSequence, inputFile)) {
            for (auto a : strProvidedSequence) {
                providedSequence.push_back(a == '1');
            }
        }
        else if (input == "scale_A") {
            for (size_t crctr = 0u; crctr < nTraits; ++crctr) {
                inputFile >> scaleA[crctr];
                std::clog   << "parameter " << "scale_A[" << crctr
                << "] set to " << scaleA[crctr] << '\n';
            }
        }
        else if (input == "scale_D") {
            for (size_t crctr = 0u; crctr < nTraits; ++crctr) {
                inputFile >> scaleD[crctr];
                std::clog   << "parameter " << "scale_D[" << crctr
                << "] set to " << scaleD[crctr] << '\n';
            }
        }
        else if (input == "scale_I") {
            for (size_t crctr = 0u; crctr < nTraits; ++crctr) {
                inputFile >> scaleI[crctr];
                std::clog   << "parameter " << "scale_I[" << crctr
                << "] set to " << scaleI[crctr] << '\n';
            }
        }
        else if (input == "scale_E") {
            for(size_t crctr = 0u; crctr < nTraits; ++crctr) {
                inputFile >> scaleE[crctr];
                std::clog   << "parameter " << "scale_E[" << crctr
                << "] set to " << scaleE[crctr] << '\n';
            }
        }
        else if (read(input, "mutation_rate", mutationRate, inputFile));
        else if (read(input, "genome_size_cm", genomeLength, inputFile));
        else if (read(input, "female_heterogamety", isFemaleHeterogamy, inputFile));
        else if (read(input, "dispersal_rate", dispersalRate, inputFile));
        else if (read(input, "beta", birthRate, inputFile));
        else if (read(input, "prob_survival", survivalProb, inputFile));
        else if (read(input, "habitat_asymmetry", habitatAsymmetry, inputFile));
        else if (read(input, "sel_coeff_ecol", ecoSelCoeff, inputFile));
        else if (read(input, "preference_strength", matePreferenceStrength, inputFile));
        else if (read(input, "preference_cost", mateEvaluationCost, inputFile));
        else if (read(input, "network_skewness", networkSkewness, inputFile));
        else if (input == "t_end") {
            inputFile >> tBurnIn >> tEndSim;
            std::clog   << "Burn-in period  " << tBurnIn << " generations \n";
            std::clog   << "Simulation time " << tEndSim << " generations \n";
        }
        else if (input == "t_dat") {
            inputFile >> tGetDat >> tSavDat;
            std::clog   << "Data collected every " << tGetDat << " generations \n";
            std::clog   << "Data stored every " << tSavDat << " generations \n";
        }
        else {
            throw std::runtime_error("Unknown parameter " + input);
        }
    }
    std::clog << "Parameters were read in successfully\n";
}

void ParameterSet::newArchitectureFileName()
{
    std::ostringstream oss;
    oss << "architecture_" << seed << ".txt";
    setArchitectureFileName(oss.str());
}


void ParameterSet::setArchitectureFileName(const std::string &filename)
{
    architectureFileName = filename;
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