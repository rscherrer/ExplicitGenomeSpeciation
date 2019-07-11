set (Boost_USE_STATIC_LIBS OFF)
set(BOOST_ROOT D:/boost_1_68_0)
find_package (Boost REQUIRED COMPONENTS unit_test_framework)
include_directories (${Boost_INCLUDE_DIRS})

add_executable (ExplicitGenomeSpeciationTestRun utilsTest.cpp  doMainTest.cpp ParameterSetTest.cpp TestFixtureCreateParameterFiles.h TestFixtureCreateArchitectureFile.h GeneticArchitectureTest.cpp PopulationTest.cpp TestFixtureDefaultPopulationSartingKit.h TestFixtureMakeAPopulation.h)
target_link_libraries (ExplicitGenomeSpeciationTestRun ${Boost_LIBRARIES})
target_link_libraries (ExplicitGenomeSpeciationTestRun lib)