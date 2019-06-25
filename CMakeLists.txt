cmake_minimum_required(VERSION 3.13)
project(ExplicitGenomeSpeciation)
set(CMAKE_CXX_STANDARD 14)

message(STATUS "Start running cmake...")

set(Boost_DEBUG on)

set(BOOST_ROOT C:/boost_1_68_0/stage)
set(BOOST_INCLUDEDIR D:/boost_1_68_0/stage/lib)

find_package(Boost 1.68.0 COMPONENTS system filesystem REQUIRED)

if(Boost_FOUND)

    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")

    include_directories(${Boost_INCLUDE_DIRS})

endif()

add_executable(ExplicitGenomeSpeciation main.cpp ParameterSet.cpp GeneticArchitecture.cpp Individual.cpp Population.cpp OutputFile.cpp random.cpp utils.cpp)

if(Boost_FOUND)

    target_link_libraries(ExplicitGenomeSpeciation ${Boost_LIBRARIES})

endif()