cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
cmake_policy(SET CMP0022 NEW)

if(NOT CMAKE_BUILD_TYPE)
    message("Using default build type: Release")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

project(ALPSCoreCTINT CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_SOURCE_DIR})

#To include alps/fastupdate
#set(CMAKE_INCLUDE_CURRENT_DIR on)

#SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Disable in-source builds
if (${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    message(FATAL_ERROR "In source builds are disabled. Please use a separate build directory")
endif()

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

# RPATH fix
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
else()
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

#policy update CMP0042
if(APPLE)
    set(CMAKE_MACOSX_RPATH ON)
endif()


# The project relies on the ALPSCore package. If ALPSCore is not found
# automatically, specify its location using:
# export ALPSCore_DIR=/location/to/ALPSCORE/
find_package(ALPSCore REQUIRED)

#Find MPI
include(EnableMPI)


# Eigen is provided by ALPSCore
if (NOT ALPSCore_HAS_EIGEN_VERSION)
  find_package(Eigen3 3.2.8 REQUIRED)
endif()

#ALPSCore disable debug for gf library
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DBOOST_DISABLE_ASSERTS -DNDEBUG")

#include directories
set(CTINT_LIBRARY_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR} ${MPI_CXX_INCLUDE_PATH} ${MPI_C_INCLUDE_PATH} ${CMAKE_SOURCE_DIR}/include) #rest taken care of by libraries dependencies
include_directories(${CTINT_LIBRARY_INCLUDE_DIRS})

#source files
set (LIB_FILES
        src/fouriertransform.C
        src/observables.ipp
        src/io.ipp
        src/interaction_expansion.ipp
        src/util.cpp
        src/legendre.cpp
        src/legendre.h
        src/fastupdate_formula.h
        src/U_matrix.cpp
        src/update_statistics.h
        src/wang_landau.h
        src/measurements.hpp
        src/selfenergy.ipp
        src/operator.cpp)

#static library
add_library(alpscore_ctint ${LIB_FILES})

#Compiler dependent libraries
set(EXTRA_LIBS "")

#Set link libraries
list(APPEND CTINT_LIBRARIES ${ALPSCore_LIBRARIES} ${MPI_CXX_LIBRARIES} ${EXTRA_LIBS})
if(DEFINED GOOGLE_PROFILER)
  list(APPEND CTINT_LIBRARIES profiler)
endif(DEFINED GOOGLE_PROFILER)
target_link_libraries(alpscore_ctint ${CTINT_LIBRARIES})

#executable
ADD_EXECUTABLE(interaction_real src/main_real.cpp)
target_link_libraries(interaction_real alpscore_ctint ${ALPSCore_LIBRARIES} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES} ${EXTRA_LIBS})

#testing setup
option(Testing "Enable testing" ON)
#include(EnableGtests) #defined in ./cmake
include_directories("test")

function(add_gtest test link_lib)
    if (TestXMLOutput)
        set (test_xml_output --gtest_output=xml:${test}.xml)
    endif(TestXMLOutput)

    set(source "${ARGV1}/${test}.cpp")
    set(gtest_src "${ARGV1}/gtest_main.cc;${ARGV1}/gtest-all.cc")

    add_executable(${test} ${source} ${gtest_src})
        message("DEBUG " ${ARGV2})
    target_link_libraries(${test} ${ARGV2})
    add_test(NAME ${test} COMMAND ${test} ${test_xml_output})
endfunction(add_gtest)

#set (GTEST_SOURCE src/legendre.cpp src/util.cpp src/operator.cpp src/U_matrix.cpp)
#add_executable(gtest ./test/gtest.cpp ./test/gtest-all.cc ./test/gtest_main.cc ${GTEST_SOURCE})
#target_link_libraries(gtest ${CTINT_LIBRARIES})

list(APPEND LINK_ALL alpscore_ctint ${ALPSCore_LIBRARIES} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES} ${EXTRA_LIBS})
list(APPEND LINK_REST ${ALPSCore_LIBRARIES} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES} ${EXTRA_LIBS})
#set(unittest_src unittests)
enable_testing(test)

add_gtest(fastupdate test ${LINK_REST})

#foreach(test ${unittest_src})
    #add_gtest(${test} test)
#endforeach(test)


