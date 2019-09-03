#include <vector>
#include <cstddef>
#include <string>
#include <stddef.h>

typedef std::vector<double> dVector;
typedef std::vector<size_t> uVector;
typedef std::vector<bool> Haplotype;

double sqr(const double&);

dVector ones(const size_t &n);

dVector zeros(const size_t&);

uVector uzeros(const size_t&); // unsigned zeros

Haplotype falses(const size_t&);

double sum(dVector&);

size_t argmin(dVector&);

size_t sumbool(Haplotype&);

size_t sumu(uVector&);
