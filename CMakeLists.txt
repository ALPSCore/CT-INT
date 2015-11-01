cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(impsolver)

# find ALPS Library
if(NOT ALPS_DIR AND ALPS_ROOT_DIR)
  set(ALPS_DIR ${ALPS_ROOT_DIR}/share/alps)
endif(NOT ALPS_DIR AND ALPS_ROOT_DIR)
if(NOT ALPS_DIR AND DEFINED ENV{ALPS_ROOT})
  set(ALPS_DIR $ENV{ALPS_ROOT}/share/alps)
endif(NOT ALPS_DIR AND DEFINED ENV{ALPS_ROOT})

#set(ALPS_DMFT_INCLUDE_DIR "" ${ALPS_DMFT_INCLUDE_DIR})

set(PREVENT_ALPS_COMPILERS ON)
find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)

set(ALPS_LIBRARY_DIRS ${ALPS_LIBRARY_DIRS} CACHE STRINGS "ALPS library dirs")
set(ALPS_INCLUDE_DIRS ${ALPS_INCLUDE_DIRS} CACHE STRINGS "ALPS include dirs")
if(ALPS_HDF5_INCLUDE_DIR) # support for old ALPS revisions
  set(HDF5_INCLUDE_DIR ${ALPS_HDF5_INCLUDE_DIR} CACHE STRINGS "")
endif(ALPS_HDF5_INCLUDE_DIR)
if(ALPS_HDF5_LIBRARIES) # support for old ALPS revisions
  set(HDF5_LIBRARIES ${ALPS_HDF5_LIBRARIES} CACHE STRINGS "")
endif(ALPS_HDF5_LIBRARIES)

mark_as_advanced(ALPS_LIBRARY_DIRS ALPS_INCLUDE_DIRS HDF5_INCLUDE_DIR HDF5_LIBRARIES)

if (ALPS_PYTHON_INTERPRETER)
  set(ALPS_INCLUDE_DIRS ${ALPS_INCLUDE_DIRS} ${ALPS_PYTHON_INCLUDE_DIRS} ${ALPS_PYTHON_NUMPY_INCLUDE_DIR})
  list(APPEND CHYB_LIBRARIES ${ALPS_PYTHON_LIBRARY})
endif (ALPS_PYTHON_INTERPRETER)

# enable C and C++ compilers
enable_language(CXX)

#include directories
set(CHYB_LIBRARY_INCLUDE_DIRS ${PROJECT_SOURCE_DIR} ${ALPS_INCLUDE_DIRS} $ENV{BLAS_LAPACK_INCLUDE_DIRS})
message(CHYB_LIBRARY_INCLUDE_DIRS "${CHYB_LIBRARY_INCLUDE_DIRS}")
include_directories(${CHYB_LIBRARY_INCLUDE_DIRS})

#library directories
set(CHYB_LIBRARY_DIRS ${ALPS_LIBRARY_DIRS} $ENV{BLAS_LAPACK_DIRS} CACHE STRINGS "CHYB library dirs")
message(CHYB_LIBRARY_DIRS "${CHYB_LIBRARY_DIRS}")
link_directories(${CHYB_LIBRARY_DIRS})

#executable
set(ALPS_ROOT ${ALPS_LIBRARY_DIRS} $ENV{BLAS_LAPACK_DIRS} CACHE STRINGS "CHYB library dirs")
set (DMFT_QMC_INTERACTION_EXPANSION3_SOURCE
       fouriertransform.C
       auxiliary.cpp
       observables.cpp
       fastupdate.cpp
       selfenergy.cpp
       solver.cpp
       io.cpp
       splines.cpp
       interaction_expansion.cpp
       measurements.cpp
       model.cpp
       util.cpp
       legendre.cpp legendre.h multi_vertex_fastupdate.cpp multi_vertex_fastupdate.h fastupdate_formula.h U_matrix.cpp)
add_executable(interaction main.cpp ${DMFT_QMC_INTERACTION_EXPANSION3_SOURCE})

set (GTEST_SOURCE legendre.cpp util.cpp)
add_executable(gtest gtest.cpp ./gtest/gtest-all.cc ./gtest/gtest_main.cc ${GTEST_SOURCE})

#libraries
list(APPEND CHYB_LIBRARIES alps boost $ENV{BLAS_LAPACK_LIBRARIES})
message(CHYB_LIBRARIES "${CHYB_LIBRARIES}")
target_link_libraries(interaction ${CHYB_LIBRARIES})
target_link_libraries(gtest ${CHYB_LIBRARIES})

message("CXX, CC=")
message(STATUS "${CMAKE_CXX_COMPILER}")
message(STATUS "${CMAKE_C_COMPILER}")
