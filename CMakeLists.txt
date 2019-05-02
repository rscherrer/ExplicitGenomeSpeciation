cmake_minimum_required(VERSION 3.10)
project(ExplicitGenomeSpeciation)

set(CMAKE_CXX_STANDARD 11)
set(SOURCE_FILES ../ExplicitGenomeSpeciation/main.cpp ../ExplicitGenomeSpeciation/random.cpp ../ExplicitGenomeSpeciation/individual.cpp ../ExplicitGenomeSpeciation/architecture.cpp ../ExplicitGenomeSpeciation/analysis.cpp)

add_executable(ExplicitGenomeSpeciation ../ExplicitGenomeSpeciation/main.cpp ../ExplicitGenomeSpeciation/random.cpp ../ExplicitGenomeSpeciation/random.h ../ExplicitGenomeSpeciation/individual.cpp ../ExplicitGenomeSpeciation/individual.h ../ExplicitGenomeSpeciation/analysis.cpp ../ExplicitGenomeSpeciation/analysis.h ../ExplicitGenomeSpeciation/architecture.cpp)