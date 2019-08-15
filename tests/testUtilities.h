#ifndef EXPLICITGENOMESPECIATION_TESTUTILITIES_H
#define EXPLICITGENOMESPECIATION_TESTUTILITIES_H

#include <boost/test/included/unit_test.hpp>
#include <vector>


/// Function to compare vectors of doubles to their expectedectations
void checkVectorOfDoubles(std::vector<double> expected,
 std::vector<double> realized, const double &factor = 1.0e6)
{
    for (size_t i = 0u; i < realized.size(); ++i) {
        realized[i] = round(realized[i] * factor);
        expected[i] = round(expected[i] * factor);
    }

    BOOST_CHECK_EQUAL_COLLECTIONS(realized.begin(), realized.end(),
     expected.begin(), expected.end());
}





struct GenFixture {

    GenFixture() :
        pars(ParameterSet()),
        arch(GeneticArchitecture(pars)),
        genome(arch.getGenome())
    {}
    ~GenFixture() {}

    ParameterSet pars;
    GeneticArchitecture arch;
    Genome genome;

};
#endif
