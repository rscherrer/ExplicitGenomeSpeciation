#ifndef EXPLICITGENOMESPECIATION_DOMAIN_H
#define EXPLICITGENOMESPECIATION_DOMAIN_H

#include "Population.h"
#include <vector>
#include <string>


/// Function to run the simulation
void runSimulation(size_t&, Population&, const size_t&, const double&,
 const double&);

/// Function to run the program
int doMain(const std::vector<std::string>&);


#endif //EXPLICITGENOMESPECIATION_DOMAIN_H
