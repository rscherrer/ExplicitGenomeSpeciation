#ifndef EXPLICITGENOMESPECIATION_TESTFIXTUREMAKEAPOPULATION_H
#define EXPLICITGENOMESPECIATION_TESTFIXTUREMAKEAPOPULATION_H

#include "Population.h"
#include "Individual.h"
#include "TestFixtureDefaultPopulationSartingKit.h"
#include <cassert>

// This fixture is to create a population with specific features.
// The population is not initialized upon construction (that is why it is in a vector)
// There are several ways to initialize the population and so the fixture can be used for different tests
// Different functions are provided in the structure to set-up a population that meets the needs of a given test
// Member variables are also there that can be used for this setup
// The common feature is that all populations will be constructed on the same parameters and architecture as loaded with the sarting kit

struct makeADefaultPopulation
{

    makeADefaultPopulation() {}
    ~makeADefaultPopulation() {}

    createPopulationStartingKit startingKit;
    std::vector<Population> populations;  // This is just a hack not to have to initialize the population upon construction of the structure, there is just one population!

    // Fixtures to allow testing habitat sorting
    std::vector<size_t> expectedInitialHabitats {0u, 0u, 0u, 0u};
    std::vector<size_t> expectedSortedHabitats {0u, 0u, 0u, 0u};
    std::vector<size_t> sortedHabitats;

    
    // Function to create a habitat-sorted population of nIndiv individuals based on explicitely provided initial habitats
    void setPopulationWithSortedHabitats()
    {
        
        // Check that there are indeed a list of initial habitats to work with
        assert(expectedInitialHabitats.size() > 0u);
        assert(expectedInitialHabitats.size() == expectedSortedHabitats.size());
        
        // Set the number of individuals in the population to be created
        const size_t nIndiv = expectedInitialHabitats.size();
        
        // Make a vector of individuals with default parameters and genetic attributes
        std::vector<PInd> individuals;
        for (size_t i = 0u; i != nIndiv; ++i) individuals.push_back(new Individual(startingKit.parameters, startingKit.geneticArchitecture));

        // Set their habitats to what we want
        for (size_t i = 0u; i != nIndiv; ++i) {
            individuals[i]->resetHabitat(expectedInitialHabitats[i]);
        }

        // Put all these individuals into a Population so we can use the sortByHabitat member function
        auto population = new Population(individuals, startingKit.parameters);

        // Sort them using the sortByHabitat function from class Population
        population->sortByHabitat();

        // Record the resulting habitats
        sortedHabitats = population->getHabitatVector();
        
        // Check the size compatibility between expected and observed sorted habitats
        assert(sortedHabitats.size() == expectedSortedHabitats.size());

        // Add the population to the unit length vector of populations (again, just a hack to be able to construct the structure without init the population)
        populations.push_back(*population);
    }

};


#endif //EXPLICITGENOMESPECIATION_TESTFIXTUREMAKEAPOPULATION_H

