cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
cmake_policy(SET CMP0022 NEW)

if(NOT CMAKE_BUILD_TYPE)
    message("Using default build type: Release")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

project(ALPSCoreCTINT CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_SOURCE_DIR})

#To include alps/fastupdate
set(CMAKE_INCLUDE_CURRENT_DIR on)

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
        fouriertransform.C
        observables.ipp
        io.ipp
        interaction_expansion.ipp
        util.cpp
        legendre.cpp legendre.h fastupdate_formula.h U_matrix.cpp update_statistics.h wang_landau.h measurements.hpp selfenergy.ipp operator.cpp)

#static library
add_library(alpscore_ctint ${LIB_FILES})
#set_target_properties(alpscore_ctint PROPERTIES PUBLIC_HEADER src/solver.hpp)

#Compiler dependent libraries
set(EXTRA_LIBS "")

#Set link libraries
list(APPEND CTINT_LIBRARIES ${ALPSCore_LIBRARIES} ${MPI_CXX_LIBRARIES} ${EXTRA_LIBS})
if(DEFINED GOOGLE_PROFILER)
  list(APPEND CTINT_LIBRARIES profiler)
endif(DEFINED GOOGLE_PROFILER)
target_link_libraries(alpscore_ctint ${CTINT_LIBRARIES})

#executable
ADD_EXECUTABLE(interaction_real main_real.cpp)
target_link_libraries(interaction_real alpscore_ctint ${ALPSCore_LIBRARIES} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES} ${EXTRA_LIBS})

#Gtest
set (GTEST_SOURCE legendre.cpp util.cpp operator.cpp U_matrix.cpp)
add_executable(gtest gtest.cpp ./gtest/gtest-all.cc ./gtest/gtest_main.cc ${GTEST_SOURCE})
target_link_libraries(gtest ${CTINT_LIBRARIES})
