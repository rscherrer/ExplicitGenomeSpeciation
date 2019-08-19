#ifndef EXPLICITGENOMESPECIATION_TESTUTILITIES_H
#define EXPLICITGENOMESPECIATION_TESTUTILITIES_H

#include <boost/test/included/unit_test.hpp>
#include <vector>


/// Function to compare vectors of doubles to their expectedectations
void compareVectors(std::vector<double> expected,
 std::vector<double> realized, const double &factor = 1.0e6)
{
    for (size_t i = 0u; i < realized.size(); ++i) {
        realized[i] = round(realized[i] * factor);
        expected[i] = round(expected[i] * factor);
    }

    BOOST_CHECK_EQUAL_COLLECTIONS(realized.begin(), realized.end(),
     expected.begin(), expected.end());
}


#endif
