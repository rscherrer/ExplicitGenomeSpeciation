project(library)

set(
        HEADER_FILES
        doMain.h
        GeneticArchitecture.h
        Individual.h
        OutputFile.h
        ParameterSet.h
        Population.h
        random.h
        utils.h
)

set(
        SOURCE_FILES
        doMain.cpp
        GeneticArchitecture.cpp
        Individual.cpp
        OutputFile.cpp
        ParameterSet.cpp
        Population.cpp
        random.cpp
        utils.cpp
)


add_library(library STATIC ${SOURCE_FILES} ${HEADER_FILES})