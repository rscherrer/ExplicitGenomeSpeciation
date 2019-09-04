#ifndef EXPLICITGENOMESPECIATION_DOMAIN_H
#define EXPLICITGENOMESPECIATION_DOMAIN_H

#include "Population.h"
#include <vector>
#include <string>
#include <stddef.h>


typedef std::vector<Population> vecPop;
typedef std::vector<Network> MultiNet;
typedef std::vector<std::string> vecStr;

/// Function to run the simulation
size_t runSimulation(vecPop&, const size_t&, const size_t&, const double&,
 const double&, const double&, const double&, const Genome&, const MultiNet&,
  const bool& = false);

/// Function to run the program
int doMain(const vecStr&);


#endif
