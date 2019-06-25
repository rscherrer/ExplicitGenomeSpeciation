cmake_minimum_required(VERSION 3.13)
project(ExplicitGenomeSpeciation)
set(CMAKE_CXX_STANDARD 14)

message(STATUS "Start running cmake...")

set(Boost_DEBUG 1)
set(BOOST_ROOT D:/boost_1_68_0)
set(BOOST_LIBRARYDIR D:/boost_1_68_0/stage/lib)

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

if (Boost_FOUND)

    include_directories(${Boost_INCLUDE_DIR})

endif()

add_executable(
        ExplicitGenomeSpeciation
        main.cpp
        ParameterSet.cpp
        GeneticArchitecture.cpp
        Individual.cpp
        Population.cpp
        OutputFile.cpp
        random.cpp
        utils.cpp
        utilsTest.cpp
)

if(Boost_FOUND)

    target_link_libraries(ExplicitGenomeSpeciation ${Boost_LIBRARIES})

endif()