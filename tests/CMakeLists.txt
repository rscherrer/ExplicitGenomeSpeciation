set(BOOST_ROOT D:/boost_1_68_0)
find_package(Boost REQUIRED COMPONENTS unit_test_framework)
include_directories(${Boost_INCLUDE_DIR})

add_executable(ExplicitGenomeSpeciation_test tests.cpp doMainTest.cpp ParameterSetTest.cpp PopulationTest.cpp GeneticArchitectureTest.cpp utilsTest.cpp)

target_link_libraries(ExplicitGenomeSpeciation_test library)