#ifndef EXPLICITGENOMESPECIATION_TESTFIXTURECREATEARCHITECTUREFILE_H
#define EXPLICITGENOMESPECIATION_TESTFIXTURECREATEARCHITECTUREFILE_H

#include <fstream>
#include <string>


/// Test fixture for creating a dummy valid genetic architecture
struct createValidArchitectureFile
{

    createValidArchitectureFile() :
    validArchitectureFileName("valid_architecture.txt"),
    expectedChromosomeSizes( {0.5, 1.0} )
    {

        // Write a dummy architecture file
        validArchitectureFile.open(validArchitectureFileName);
        validArchitectureFile << "2\n"
                              << "2\n"
                              << "3\n"
                              << "0\n"
                              << "1\n"
                              << "2\n"
                              << "2\n"
                              << "0.5\t1.0\n\n"
                              << "0\t0\t0.1\t0.01\t0.0\n"
                              << "1\t0\t0.2\t-0.01\t0.0\n"
                              << "2\t1\t0.3\t0.01\t0.0\n"
                              << "3\t1\t0.4\t-0.01\t0.0\n"
                              << "4\t2\t0.5\t0.01\t0.0\n"
                              << "5\t2\t0.6\t-0.01\t0.0\n"
                              << "6\t2\t0.7\t0.01\t0.0\n"
                              << "2\t3\t0.01\n"
                              << "4\t5\t-0.01\n"
                              << "4\t6\t0.01\n";
        validArchitectureFile.close();

        // Make sure that the model parameters match
        parameters.nEcoLoci = 2u;
        parameters.nMatLoci = 2u;
        parameters.nNtrLoci = 3u;
        parameters.nEcoInteractions = 0u;
        parameters.nMatInteractions = 1u;
        parameters.nNtrInteractions = 2u;
        parameters.nChromosomes = 2u;
    }
    ~createValidArchitectureFile() {}

    std::ofstream validArchitectureFile;
    std::string validArchitectureFileName;

    ParameterSet parameters;

    std::vector<double> expectedChromosomeSizes;

};

#endif //EXPLICITGENOMESPECIATION_TESTFIXTURECREATEARCHITECTUREFILE_H
