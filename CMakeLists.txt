cmake_minimum_required(VERSION 3.0)

## Set our project name
project(hakmer-ng-redesign CXX C)

## Use all the *.cpp files we found under this folder for the project

file(GLOB_RECURSE SRCS src/*.cpp)
file(GLOB_RECURSE HDRS src/*.hpp)

INCLUDE_DIRECTORIES(include)
LINK_DIRECTORIES(lib)

add_compile_options(-std=c++14)

if (NOT (CMAKE_BUILD_TYPE MATCHES Debug))
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
endif()

if (CMAKE_BUILD_TYPE MATCHES Debug)
# Enfore debug info for valgrind
set( CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        -g" )
endif()

include_directories( "src" )

add_executable( hakmer-ng-redesign ${SRCS} ${HDRS} )
target_link_libraries ( hakmer-ng-redesign sdsl divsufsort divsufsort64 )
