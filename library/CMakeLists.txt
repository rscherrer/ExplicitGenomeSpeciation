project(library)

set(
        HEADER_FILES
        function1.h
        function2.h
        function3.h
        function4.h
)

set(
        SOURCE_FILES
        function1.cpp
        function2.cpp
        function3.cpp
        function4.cpp
)

add_library(library STATIC ${SOURCE_FILES} ${HEADER_FILES})