project(lib)

set(
        HEADER_FILES
        GeneticArchitecture.h
        Individual.h
        OutputFile.h
        ParameterSet.h
        Population.h
        random.h
        utils.h
        doMain.h
)

set(
        SOURCE_FILES
        GeneticArchitecture.cpp
        Individual.cpp
        OutputFile.cpp
        ParameterSet.cpp
        Population.cpp
        random.cpp
        utils.cpp
        doMain.cpp
)

add_library(lib STATIC ${SOURCE_FILES} ${HEADER_FILES})