#include "doMain.h"
#include "ParameterSet.h"
#include "GeneticArchitecture.h"
#include "Population.h"
#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <cassert>


/// Program to run the main function
int doMain(const std::vector<std::string> &args)
{

    try
    {
        // Return an error if there are more than one argument
        if (args.size() > 2u)
        {
            throw std::runtime_error("More than one argument was supplied");
        }

        // Create a default parameter set
        ParameterSet parameters;
        parameters.setDefaultSeed();

        /*
        // Update parameters from a parameter file if needed
        if (args.size() == 2u)
        {

            // Open the parameter file
            std::ifstream inputFile(args[1u]);
            if (!inputFile.is_open()) {
                throw std::runtime_error("Unable to open parameter file");
            }

            std::clog << "Reading parameters from file " << args[1u] << '\n';

            // Update parameters based on what is read from the parameter file
            parameters.readParameters(inputFile);

            std::clog << "Parameters were read in successfully\n";
        }
        else
        {
            std::clog << "Using default parameters\n";
        }
         */

        GeneticArchitecture geneticArchitecture = GeneticArchitecture(parameters.nChromosomes);





        /*
        // Declare a genetic architecture
        GeneticArchitecture geneticArchitecture(parameters.nChromosomes);

        // Set the genetic architecture
        if (parameters.isGenerateArchitecture)
        {

            // Create a new genetic architecture if needed
            geneticArchitecture.generateGeneticArchitecture(parameters);

            // Create a file name for the newly created architecture
            parameters.newArchitectureFileName();

            // Store the newly created genetic architecture
            geneticArchitecture.storeGeneticArchitecture(parameters);

            std::clog << "New genetic architecture saved as " << parameters.architectureFileName << '\n';

        }

        else
        {

            std::clog << "Loading a genetic architecture from file " << parameters.architectureFileName << '\n';

            // Otherwise load the genetic architecture from a genetic architecture file
            geneticArchitecture.loadGeneticArchitecture(parameters);

        }

        // Record start of simulation
        auto tStart = std::chrono::system_clock::now();

        std::clog << "Creating initial population.";
        auto population = new Population(parameters, geneticArchitecture);
        assert(population->getPopSize() == parameters.initialPopSize);
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

                // Open question: is is better to keep all individuals in a single vector and order them by habitat,
                // Or have two separate vectors, subpopulations 1 and 2, members of class Population?

                population->dispersal(parameters);
                std::clog << " Dispersal\n";

                population->sortByHabitat();
                std::clog << " Individuals sorted by habitat\n";

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
        */





    }
    catch (const std::runtime_error &err)
    {
        std::cerr << "Exception: " << err.what() << '\n';
        return 1;
    }

    return 0;
}
