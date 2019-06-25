cmake_minimum_required(VERSION 3.13)
project(ExplicitGenomeSpeciation)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES main.cpp)
add_executable(ExplicitGenomeSpeciation ${SOURCE_FILES})

include_directories(lib)
add_subdirectory(lib)

target_link_libraries(ExplicitGenomeSpeciation lib)

add_subdirectory(test)
