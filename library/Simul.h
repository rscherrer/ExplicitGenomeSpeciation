#ifndef EXPLICITGENOMESPECIATION_SIMUL_H
#define EXPLICITGENOMESPECIATION_SIMUL_H

#include "Param.h"
#include "GenArch.h"
#include "Random.h"
#include "Types.h"
#include "MetaPop.h"
#include "Collector.h"
#include <stddef.h>
#include <iostream>

int simulate(const vecStrings&);

bool timetosave(const int &t, const Param &p);

#endif
