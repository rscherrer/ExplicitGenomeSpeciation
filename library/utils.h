#include <vector>
#include <cstddef>
#include <string>

typedef std::vector<double> Vector;
typedef std::vector<size_t> Uector;
typedef std::vector<bool> Haplotype;

double sqr(const double&);

std::vector<double> zeros(const size_t&);

std::vector<size_t> uzeros(const size_t&); // unsigned zeros

double sum(std::vector<double>&);

size_t argmin(std::vector<double>&);

size_t sumbool(std::vector<bool>&);

size_t sumu(Uector&);

void printvec(Vector, std::string = "");

void printuec(Uector, std::string = "");

void printbvec(Haplotype, std::string = "");
