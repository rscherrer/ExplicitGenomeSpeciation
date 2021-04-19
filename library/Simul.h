#ifndef EXPLICITGENOMESPECIATION_SIMUL_H
#define EXPLICITGENOMESPECIATION_SIMUL_H

#include "Param.h"
#include "GenArch.h"
#include "Random.h"

#include "MetaPop.h"
#include "Collector.h"
#include "Printer.h"
#include "Freezer.h"
#include "Pedigree.h"
#include <stddef.h>
#include <iostream>

int simulate(const std::vector<std::string>&);

bool timetosave(const int &t, const Param &p);
bool timetofreeze(const int &t, const Param &p);

#endif
