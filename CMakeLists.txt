cmake_minimum_required(VERSION 3.13)
project(ExplicitGenomeSpeciation)

set(CMAKE_CXX_STANDARD 11)

add_subdirectory(lib/googletest-master)
include_directories(lib/googletest-master/googletest/include)
include_directories(lib/googletest-master/googlemock/include)

add_executable(ExplicitGenomeSpeciation main.cpp random.cpp Individual.cpp analysis.cpp Tests/test.cpp ParameterSet.cpp ParameterSet.h BufferBox.cpp BufferBox.h Population.cpp Population.h Genome.cpp Genome.h Buffer.cpp Buffer.h ioroutines.h accessories.h)
target_link_libraries(ExplicitGenomeSpeciation gtest gtest_main)

