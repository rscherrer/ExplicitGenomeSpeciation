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

#include "ParameterSet.h"
#include "GeneticArchitecture.h"
#include "Individual.h"
#include "Population.h"
#include "OutputFile.h"
#include "random.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <chrono>
#include <string>
#include <vector>
#include <list>
#include <cassert>


/*=======================================================================================================
                                            main()
========================================================================================================*/

int doMain(const std::vector<std::string> &argv)
{
    OutputFile logFile;
    OutputFile datFile;
    OutputFile arcFile;

    try {

        // *** Preliminaries ***

        const int argc = argv.size();

        // Set default parameters
        ParameterSet parameters;
        parameters.setDefaultSeed();

        // Update parameters if a parameter file is provided
        if (argc == 2) {
            std::string parameterFileName = argv[1];
            parameters.readParameters(parameterFileName);
        }

        // Exception handling
        if (argc < 1 || argc > 2) {
            throw std::runtime_error("Invalid number of program arguments in main()");
        }

        // Seed the random number generator
        rnd::rng.seed(parameters.seed);

        // Make an architecture file if architecture is to be generated
        if (parameters.isGenerateArchitecture) {
            parameters.newArchitectureFileName();
        }
        assert(parameters.architectureFileName.size() > 1u);

        // Initialize genetic architecture
        GeneticArchitecture geneticArchitecture = GeneticArchitecture(parameters);

        // Output files
        std::clog << "Opening output files.";
        logFile.open(parameters.seed, ".log");
        datFile.open(parameters.seed, ".dat");
        arcFile.open(parameters.seed, "_fossil_record.txt");
        if (!(logFile.file.is_open() && datFile.file.is_open() && arcFile.file.is_open())) {
            throw std::runtime_error("Unable to open output files in main()");
        }
        std::clog << "..done\n";

        // Buffers
        // std::clog << "Opening data buffers.";
        // BufferBox bufferPointers = BufferBox(parameters, population, genome);
        // std::clog << "..done\n";

        // Store parameter values
        std::clog << "Storing parameter values..";
        if (argc == 1) {
            logFile.file << "Default parameters values\n";
        }
        else {
            logFile.file << "Parameters imported from file " << argv[1] << '\n';
        }
        logFile.writeParameters(parameters);
        std::clog << "..done\n";

        // Write data file header
        std::clog << "Writing data file header.";
        datFile.writeHeader();
        std::clog << "..done\n";

        // *** Simulation ***

        // Record start of simulation
        auto tStart = std::chrono::system_clock::now();

        std::clog << "Creating initial population.";
        auto population = new Population(parameters, geneticArchitecture);
        std::clog << "..done\n";

        // Enter simulation loop
        std::clog << "Entering simulation loop:\n";

        // Loop through time
        for (int t = 1 - parameters.tBurnIn; t <= parameters.tEndSim; ++t) {

            bool isBurnin = t < 0;
            if (isBurnin) {

                // Burnin period
                if (population->getNResources() > 1u) {
                    population->setBurnin();
                }
                population->resourceDynamics(0u, parameters.ecoSelCoeff);
                population->reproduction(0u, parameters, geneticArchitecture);

            }
            else {

                if (population->getNResources() < 2u) {
                    population->endBurnin();
                }

                population->dispersal(parameters);
                population->sortByHabitat();

                // First habitat
                population->resourceDynamics(0u, parameters.ecoSelCoeff);
                population->reproduction(0u, parameters, geneticArchitecture);

                // Second habitat
                population->resourceDynamics(1u, parameters.ecoSelCoeff);
                population->reproduction(1u, parameters, geneticArchitecture);

            }

            population->survival(parameters.survivalProb);

            // Check for extinction
            bool isExtinct = population->getPopSize() < 2u;
            if (isExtinct) {
                std::clog << "Population size underflow at t = " << t << '\n';
                break;
            }

            // Time to analyze?

            if (t % parameters.tGetDat == 0u) {

                population->sortByHabitat();
                population->assignEcotypes();

                // Genome-wide variance decomposition
                population->decomposeVarianceAlongGenome(parameters.tiny);

                // Overall variance decomposition
                population->decomposeVariance(parameters.tiny);

                // Record data
                // population->writePopulationData(datFile);
                // population->screenshotIndividuals();
                // population->screenshotLoci();
                // population->screenshotNetworkEdges();

            }

            // if(t % parameters.tGetDat == 0u) decomposeVariance(t, parameters, bufferPointers, arcFile, datFile, population, genome);
            // if(t % parameters.tSavDat == 0u) analyseNetwork(t, parameters, population, genome);
        }

        // *** Finalization ***

        // Record end and duration of simulation
        auto tEnd = std::chrono::system_clock::now();
        std::chrono::duration<double> diff = tEnd-tStart;
        logFile.file << "Time to complete simulation: " << diff.count() << " s\n";

        // Close output files
        logFile.file.close();
        datFile.file.close();
        arcFile.file.close();

        // Free allocated memory (what needs to be deleted?)
        population->massExtinction();

    }
    catch(const std::exception &err) {
        std::cerr << "Exception: " << err.what() << '\n';
        logFile.file << "Exception: " << err.what() << '\n';
        return 1;
    }

    std::cout << "The simulation completed. Hurray!\n";

    return 0;
}

/// Program to simulate an individual-based speciation event with explicit genomic features
int main(int argc, char * argv[])
{

    // Preparation
    const std::vector<std::string> args(argv, argv + argc);

    // Run the program
    return doMain(args);

}