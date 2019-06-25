#include "doMain.h"
#include "ParameterSet.h"
#include <iostream>
#include <vector>
#include <string>


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

    }
    catch (const std::runtime_error &err)
    {
        std::cout << "Exception: " << err.what() << '\n';
        return 1;
    }

    return 0;
}
