set(BOOST_ROOT D:/boost_1_68_0)
find_package(Boost REQUIRED COMPONENTS unit_test_framework)
include_directories(${Boost_INCLUDE_DIR})

set(
        TEST_FILES
        tests.cpp
        doMainTest.cpp
        GeneticArchitectureTest.cpp
        ParameterSetTest.cpp
        PopulationTest.cpp
        utilsTest.cpp
)

add_executable(ExplicitGenomeSpeciation_test ${TEST_FILES})

target_link_libraries(ExplicitGenomeSpeciation_test library)