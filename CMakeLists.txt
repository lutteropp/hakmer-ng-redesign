cmake_minimum_required(VERSION 3.0)

## Set our project name
project(hakmer-ng-redesign CXX C)

#find_library(SDSL_LIB sdsl)
#find_library(DIVSUFSORT_LIB divsufsort)
#find_library(DIVSUFSORT64_LIB divsufsort64)

## Use all the *.cpp files we found under this folder for the project

file(GLOB_RECURSE SRCS src/*.cpp)
file(GLOB_RECURSE HDRS src/*.hpp)

add_compile_options(-std=c++14)

if (NOT (CMAKE_BUILD_TYPE MATCHES Debug))
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Add Genesis as dependency.
add_subdirectory(${PROJECT_SOURCE_DIR}/genesis)
# Use everything that Genesis exports, just to be sure that we use the same setup.
add_definitions( ${GENESIS_DEFINITIONS} )
include_directories( ${GENESIS_INCLUDE_DIR} )
set( CMAKE_C_FLAGS          "${CMAKE_C_FLAGS}          ${GENESIS_C_FLAGS}")
set( CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        ${GENESIS_CXX_FLAGS}" )
set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GENESIS_EXE_LINKER_FLAGS}" )
add_definitions(-DWITH_GENESIS -DWITH_OPENMP)
endif()

if (CMAKE_BUILD_TYPE MATCHES Debug)
# Enfore debug info for valgrind
set( CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        -g" )
endif()

include_directories( "src" )

add_executable( hakmer-ng-redesign ${SRCS} ${HDRS} )
#target_link_libraries( hakmer-ng-redesign "${SDSL_LIB}" )
#target_link_libraries( hakmer-ng-redesign "${DIVSUFSORT_LIB}" )
#target_link_libraries( hakmer-ng-redesign "${DIVSUFSORT64_LIB}" )

if (NOT (CMAKE_BUILD_TYPE MATCHES Debug))
# Link it against Genesis, and against all dependencies of Genesis.
target_link_libraries ( hakmer-ng-redesign ${GENESIS_LINK_LIBRARIES} )
endif()