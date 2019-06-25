cmake_minimum_required(VERSION 3.13)
project(ExplicitGenomeSpeciation)
set(CMAKE_CXX_STANDARD 14)

add_executable(ExplicitGenomeSpeciation main.cpp ParameterSet.cpp GeneticArchitecture.cpp Individual.cpp Population.cpp OutputFile.cpp random.cpp)

