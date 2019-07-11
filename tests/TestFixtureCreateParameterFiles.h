#ifndef EXPLICITGENOMESPECIATION_TESTFIXTURECREATEPARAMETERFILES_H
#define EXPLICITGENOMESPECIATION_TESTFIXTURECREATEPARAMETERFILES_H

#include <fstream>
#include <string>


/// Test fixture for creating a dummy but valid parameter file
struct createValidParameterFile
{

    createValidParameterFile() : validParameterFileName("valid_parameters.txt")
    {
        validParameterFile.open(validParameterFileName);
        validParameterFile << "initialPopSize 1000\nbirthRate 4.0\nscaleA 0.0 0.5 1.0\n";
        validParameterFile.close();
    }

    ~createValidParameterFile() {}

    std::ofstream validParameterFile;
    std::string validParameterFileName;

};


/// Test fixture for creating a dummy parameter file with invalid parameters in it
struct createInvalidParameterFile
{

    createInvalidParameterFile() : invalidParameterFileName("invalid_parameters.txt")
    {
        invalidParameterFile.open(invalidParameterFileName);
        invalidParameterFile << "height 1000\nweight 4.0\n";
        invalidParameterFile.close();
    }

    ~createInvalidParameterFile() {}

    std::ofstream invalidParameterFile;
    std::string invalidParameterFileName;

};

#endif //EXPLICITGENOMESPECIATION_TESTFIXTURECREATEPARAMETERFILES_H
