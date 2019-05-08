/*==================================================================================================================================
                                                     main.cpp
====================================================================================================================================

C++-code accompanying:	
		 
		(ms. in prep).

Written by:
        G. Sander van Doorn
       	Centre for Ecological and Evolutionary Studies - Theoretical Biology Group
        University of Groningen
        the Netherlands

Program version
		xx/xx/2018	:	

Instructions for compiling and running the program
		
	Versions of this program were compiled and run on Windows and Mac, using Microsoft Visual C++
	2010 and XCode, respectively. The code is written in standard C++ and should be compatible with 
	other compilers. 

=================================================================================================================================*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <chrono>
#include <string>
#include <vector>
#include <list>
#include <queue>
#include <cassert>
#include <gtest/gtest.h>
#include "random.h"
#include "Individual.h"
#include "analysis.h"
#include "Genome.h"
#include "Buffer.h"

/*=======================================================================================================
                                Parameter definitions and default values
========================================================================================================*/




/*=======================================================================================================
                                            I/O routines
========================================================================================================*/

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




/*=======================================================================================================
                                    biological model implementation
========================================================================================================*/



bool tradeOffCompare (const Individual::TradeOffPt &x, const Individual::TradeOffPt &y) {
    bool yOnLeft = y.first < y.second;
    if(x.first < x.second) {
        if(yOnLeft) return (x.first < y.first);
        else return true;
    }
    else {
        if(yOnLeft) return false;
        else return (x.second < y.second);
    }
}






/*=======================================================================================================
                                            main()
========================================================================================================*/

int main(int argc, char * argv[])
{

    #if !defined(NDEBUG)

        // In debug mode run tests only
        std::cout << "Running tests in debug mode...";
        testing::InitGoogleTest(&argc, argv);
        int testResult = RUN_ALL_TESTS();

    #else

    std::ofstream logFile; // must be accessible from catch

    try {

        // *** preliminaries ***

        // Parameters
        ParameterSet parameters;

        std::ofstream datFile, arcFile;

        // Declare pointers to addresses where buffers are stored
        // The values (the buffers) are accessed with (*pointerName)
        Buffer *bufferFreq, *bufferF_it, *bufferF_is, *bufferF_st,
                *bufferP_st, *bufferG_st, *bufferQ_st, *bufferC_st,
                *bufferVarP, *bufferVarG, *bufferVarA, *bufferVarD, *bufferVarI;

        // set parameter values
        if(argc == 1) {
            parameters.seed = rnd::set_seed(); // use default parameters values and use clock to set random seed
            parameters.generateArchitecture = true;
        }
        else if(argc == 2)
            readParameters(argv[1], parameters);
        else throw std::runtime_error("invalid number of program arguments in main()");

        Population population;
        Genome genome;

        // initialise genetic architecture
        if(parameters.generateArchitecture) {
            std::ostringstream oss;
            oss << "architecture_" << parameters.seed << ".txt";
            parameters.architecture = oss.str();
            genome.generateGeneticArchitecture(parameters);
            genome.storeGeneticArchitecture(parameters.architecture, parameters);
        }
        else genome.loadGeneticArchitecture(parameters.architecture, parameters, population);



        // open files and data buffers
        std::clog << "opening files and data buffers.";
        std::ostringstream oss;
        oss << "simulation_" << parameters.seed;
        logFile.open(oss.str() + ".log");
        datFile.open(oss.str() + ".dat");
        arcFile.open(oss.str() + "_fossil_record.txt");
        if(!(logFile.is_open() && datFile.is_open() && arcFile.is_open()))
            throw std::runtime_error("unable to open output files in main()");
        BufferBox bufferPointers;
        bufferPointers.bufferFreq = new Buffer("freq", parameters, population, genome);
        bufferPointers.bufferF_it = new Buffer("Fit", parameters, population, genome);
        bufferPointers.bufferF_is = new Buffer("Fis", parameters, population, genome);
        bufferPointers.bufferF_st = new Buffer("Fst", parameters, population, genome);
        bufferPointers.bufferP_st = new Buffer("Pst", parameters, population, genome);
        bufferPointers.bufferG_st = new Buffer("Gst", parameters, population, genome);
        bufferPointers.bufferQ_st = new Buffer("Qst", parameters, population, genome);
        bufferPointers.bufferC_st = new Buffer("Cst", parameters, population, genome);
        bufferPointers.bufferVarP = new Buffer("varP", parameters, population, genome);
        bufferPointers.bufferVarG = new Buffer("varG", parameters, population, genome);
        bufferPointers.bufferVarA = new Buffer("varA", parameters, population, genome);
        bufferPointers.bufferVarD = new Buffer("varD", parameters, population, genome);
        bufferPointers.bufferVarI = new Buffer("varI", parameters, population, genome);
        std::clog << "..done\n";

        // store parameter values
        std::clog << "storing parameter values..";
        logFile << "parameters: ";
        if(argc == 1) logFile << "default values\n";
        else logFile << "imported from file " << argv[1] << '\n';
        writeParameters(logFile, parameters);
        std::clog << "..done\n";

        // write data file header
        std::clog << "writing data file header.";
        datFile << '\t' << "pop.size"
                << '\t' << "females"
                << '\t' << "males";
        for(size_t hab = 0u; hab < parameters.nHabitat; ++hab)
            datFile << '\t' << "pop.size." << hab;
        for(size_t hab = 0u; hab < parameters.nHabitat; ++hab)
            datFile << '\t' << "attack.rate." << hab << ".1" << '\t' << "attack.rate." << hab << ".2";
        for(size_t hab = 0u; hab < parameters.nHabitat; ++hab)
            datFile << '\t' << "resource." << hab << ".1" << '\t' << "resource." << hab << ".2";
        for(size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr)
            datFile << '\t' << "pop.1.size" << crctr
                    << '\t' << "pop.2.size" << crctr
                    << '\t' << "mean.all." << crctr
                    << '\t' << "mean.1." << crctr
                    << '\t' << "mean.2." << crctr
                    << '\t' << "varP." << crctr
                    << '\t' << "varG." << crctr
                    << '\t' << "varA." << crctr
                    << '\t' << "varD." << crctr
                    << '\t' << "varI." << crctr
                    << '\t' << "Fst."  << crctr
                    << '\t' << "Pst."  << crctr
                    << '\t' << "Gst."  << crctr
                    << '\t' << "Qst."  << crctr
                    << '\t' << "Cst."  << crctr;
        datFile << '\t' << "speciation.cube.spatial.isolation"
                << '\t' << "speciation.cube.ecological.isolation"
                << '\t' << "speciation.cube.mating.isolation"
                << '\t' << "post.zygotic.isolation"<< '\n';
        std::clog << "..done\n";





        // record start of simulation
        auto tStart = std::chrono::system_clock::now();

        // *** simulation ***
        // create initial population

        std::vector<std::pair<double, double> > resourceConsumption, resourceEql;
        std::vector<std::pair<size_t, size_t> > genderCounts;
        Individual::TradeOffPt breakEvenPoint;

        std::clog << "creating initial population.";
        if(parameters.sequence.size() == parameters.nBits)
            for(size_t i = 0u; i < parameters.nIndividualInit ; ++i)
                population.individuals.push_back(new Individual(parameters.sequence, parameters, genome));
        else for(size_t i = 0u; i < parameters.nIndividualInit ; ++i)
                population.individuals.push_back(new Individual(parameters, genome));
        std::clog << "..done\n";

        // enter simulation loop
        std::clog << "entering simulation loop:\n";
        for(int t = 1 - parameters.tBurnIn; t <= parameters.tEndSim; ++t) {
            if(t > 0) {
                // default
                population.dispersal(parameters);
                population.competitionAndReproduction(0u, parameters, population, resourceConsumption, breakEvenPoint, resourceEql, genderCounts);
                population.competitionAndReproduction(1u, parameters, population, resourceConsumption, breakEvenPoint, resourceEql, genderCounts);

            }
            else population.competitionAndReproduction(0u, parameters, population, resourceConsumption, breakEvenPoint, resourceEql, genderCounts, 1u); // burn-in period
            if(population.individuals.size() < 2u) {
                std::clog << "population size underflow at t = " << t << '\n';
                break;
            }
            if(t % parameters.tGetDat == 0u) decomposeVariance(t, parameters, bufferPointers, breakEvenPoint, resourceConsumption, resourceEql, arcFile, datFile, genderCounts, population, genome);
            if(t % parameters.tSavDat == 0u) analyseNetwork(t, parameters, population, genome);
        }

        // *** finalisation ***
        // record end of simulation
        auto tEnd = std::chrono::system_clock::now();
        std::chrono::duration<double> diff = tEnd-tStart;
        logFile << "Time to complete simulation: " << diff.count() << " s\n";

        // close output file
        logFile.close();
        datFile.close();
        arcFile.close();

        // free allocated memory
        while(!population.individuals.empty()) {
            delete population.individuals.back();
            population.individuals.pop_back();
        }
    }
    catch(const std::exception &err) {
        std::cerr << "exception: " << err.what() << '\n';
        logFile << "exception: " << err.what() << '\n';
        exit(EXIT_FAILURE);
    }

    #endif


    return EXIT_SUCCESS;
}