#include <vector>
#include <cstddef>
#include <string>
#include <stddef.h>

typedef std::vector<double> vecDbl;
typedef std::vector<size_t> vecUns;
typedef std::vector<bool> Haplotype;

double sqr(const double&);

vecDbl ones(const size_t &n);

vecDbl zeros(const size_t&);

vecUns uzeros(const size_t&); // unsigned zeros

Haplotype falses(const size_t&);

double sum(vecDbl&);

size_t argmin(vecDbl&);

size_t sumbool(Haplotype&);

size_t sumu(vecUns&);
