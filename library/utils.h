#ifndef EXPLICITGENOMESPECIATION_UTILS_H
#define EXPLICITGENOMESPECIATION_UTILS_H

#include "types.h"
#include <vector>
#include <cstddef>
#include <string>
#include <stddef.h>

double sqr(const double&);

vecDbl ones(const size_t &n);

vecDbl zeros(const size_t&);

vecUns uzeros(const size_t&); // unsigned zeros

vecBool falses(const size_t&);

double sum(vecDbl&);

size_t argmin(vecDbl&);

size_t sumbool(vecBool&);

size_t sumu(vecUns&);

double size2dbl(const size_t&);

#endif
