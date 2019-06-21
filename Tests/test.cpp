#include <gtest/gtest.h>
#include <iostream>

#include "../Individual.h"
#include "../GeneticArchitecture.h"

using testing::Eq;

// Dummy test
TEST(DummyTests, OneEqualsOne)
{
    EXPECT_EQ(1, 1);
}

// Test fixture to create and clean up a parameter set
struct ParameterSetTest : testing::Test{

    ParameterSet* parameters;

    ParameterSetTest()
    {
        parameters = new ParameterSet;
    }

    virtual ~ParameterSetTest()
    {
        delete parameters;
    }

};

// Test fixture to create and clean up a genetic architecture
struct GeneticArchitectureTest : testing::Test
{

    GeneticArchitecture geneticArchitecture;

    GeneticArchitectureTest(const ParameterSet &parameters)
    {
        geneticArchitecture = GeneticArchitecture(parameters);
    }

};

TEST(geneticArchitectureTests, threeChromosomesOfEqualSizes)
{
    setChromosomeSizes(nChromosomes);
}
