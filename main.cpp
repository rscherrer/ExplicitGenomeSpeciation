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
#include "GeneticArchitecture.h"
#include "Buffer.h"
#include "Population.h"
#include "ioroutines.h"


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
    std::ofstream datFile;
    std::ofstream arcFile;

    try {

        // *** Preliminaries ***

        ParameterSet parameters;

        // Set default parameters
        if (argc == 1) {
            parameters.seed = rnd::set_seed();
            parameters.isGenerateArchitecture = true;
        }

        // Or read parameters from a file
        else if (argc == 2) {
            parameters.readParameters(argv[1]);
        }
        else {
            throw std::runtime_error("invalid number of program arguments in main()");
        }

        // Set the name of the architecture file if needed
        if (parameters.isGenerateArchitecture) {
            std::ostringstream oss;
            oss << "architecture_" << parameters.seed << ".txt";
            parameters.architectureFilename = oss.str();
        }

        // Initalize genetic architecture
        GeneticArchitecture geneticArchitecture = GeneticArchitecture(parameters);

        // Open output files
        std::clog << "Opening files and data buffers.";
        std::ostringstream oss;
        oss << "simulation_" << parameters.seed;
        logFile.open(oss.str() + ".log");
        datFile.open(oss.str() + ".dat");
        arcFile.open(oss.str() + "_fossil_record.txt");
        if (!(logFile.is_open() && datFile.is_open() && arcFile.is_open())) {
            throw std::runtime_error("Unable to open output files in main()");
        }

        // Open a buffer box (work on this later)
        BufferBox bufferPointers = BufferBox(parameters, population, genome);

        std::clog << "..done\n";

        // Store parameter values
        std::clog << "Storing parameter values..";
        logFile << "Parameters: ";
        if (argc == 1) {
            logFile << "Default values\n";
        }
        else {
            logFile << "Imported from file " << argv[1] << '\n';
        }
        parameters.writeParameters(logFile);
        std::clog << "..done\n";

        // Write data file header
        std::clog << "Writing data file header.";
        datFile << '\t' << "pop.size"
                << '\t' << "females"
                << '\t' << "males";

        for (size_t hab = 0u; hab < parameters.nHabitat; ++hab) {
            datFile << '\t' << "pop.size." << hab;
        }

        for (size_t hab = 0u; hab < parameters.nHabitat; ++hab) {
            datFile << '\t' << "attack.rate." << hab << ".1" << '\t' << "attack.rate." << hab << ".2";
        }

        for (size_t hab = 0u; hab < parameters.nHabitat; ++hab) {
            datFile << '\t' << "resource." << hab << ".1" << '\t' << "resource." << hab << ".2";
        }

        for (size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr) {
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
        }

        datFile << '\t' << "speciation.cube.spatial.isolation"
                << '\t' << "speciation.cube.ecological.isolation"
                << '\t' << "speciation.cube.mating.isolation"
                << '\t' << "post.zygotic.isolation"<< '\n';

        std::clog << "..done\n";


        // *** Simulation ***

        // Record start of simulation
        auto tStart = std::chrono::system_clock::now();

        std::clog << "Creating initial population.";
        Population population = Population(parameters);
        std::clog << "..done\n";

        // Enter simulation loop
        std::clog << "Entering simulation loop:\n";

        // Loop through time
        for (int t = 1 - parameters.tBurnIn; t <= parameters.tEndSim; ++t) {
            if (t > 0) {
                population.dispersal(parameters);
                population.sortByHabitat();
                population.resourceDynamics(0u);
                population.reproduction(0u);
                population.resourceDynamics(1u);
                population.reproduction(1u);
            }
            else {

                // Burnin period
                population.resourceDynamics(0u);
                population.reproduction(0u);
            }

            // Check for extinction
            if (population.individuals.size() < 2u) {
                std::clog << "population size underflow at t = " << t << '\n';
                break;
            }

            // Time to analyze?
            if(t % parameters.tGetDat == 0u) decomposeVariance(t, parameters, bufferPointers, arcFile, datFile, population, genome);
            if(t % parameters.tSavDat == 0u) analyseNetwork(t, parameters, population, genome);
        }

        // *** Finalisation ***

        // Record end and duration of simulation
        auto tEnd = std::chrono::system_clock::now();
        std::chrono::duration<double> diff = tEnd-tStart;
        logFile << "Time to complete simulation: " << diff.count() << " s\n";

        // Close output files
        logFile.close();
        datFile.close();
        arcFile.close();

        // Free allocated memory (what needs to be deleted?)
        while(!population.individuals.empty()) {
            delete population.individuals.back();
            population.individuals.pop_back();
        }


    }
    catch(const std::exception &err) {
        std::cerr << "Exception: " << err.what() << '\n';
        logFile << "Exception: " << err.what() << '\n';
        exit(EXIT_FAILURE);
    }

    #endif


    return EXIT_SUCCESS;
}