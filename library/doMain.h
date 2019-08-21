#ifndef EXPLICITGENOMESPECIATION_DOMAIN_H
#define EXPLICITGENOMESPECIATION_DOMAIN_H

#include "Population.h"
#include <vector>
#include <string>

typedef std::vector<Network> MultiNet;
typedef std::vector<std::string> sVector;

/// Function to run the simulation
void runSimulation(size_t&, Population&, const size_t&, const double&,
 const double&, const double&, const Genome&, const MultiNet&);

/// Function to run the program
int doMain(const sVector&);


#endif //EXPLICITGENOMESPECIATION_DOMAIN_H
