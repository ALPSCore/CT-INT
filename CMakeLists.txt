cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
cmake_policy(SET CMP0022 NEW)

if(NOT CMAKE_BUILD_TYPE)
    message("Using default build type: Release")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

project(ALPSCoreCTINT C CXX)

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

find_package(Threads)

# Eigen is provided by ALPSCore
if (NOT ALPSCore_HAS_EIGEN_VERSION)
  if (${USE_BLAS_LAPACK})
    # To use Lapack with Eigen3, version >= 3.3 is required.
    find_package(Eigen3 3.3 REQUIRED)
  else()
    find_package(Eigen3 3.2.8 REQUIRED)
  endif()
endif()

if (${USE_BLAS_LAPACK})
  add_definitions(-DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE)
  set(BLA_VENDOR Generic)
  find_package(BLAS)
  find_package(LAPACK)
  set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  set(EXTRA_LIBS ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

#ALPSCore disable debug for gf library
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DBOOST_DISABLE_ASSERTS -DNDEBUG")

if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    message( STATUS "Using Intel Compiler. Adding -fp-model strict to prevent aggressive optimization.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fp-model strict")
endif()

#include directories
set(CTINT_LIBRARY_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR} ${MPI_CXX_INCLUDE_PATH} ${MPI_C_INCLUDE_PATH} ${CMAKE_SOURCE_DIR}/include) #rest taken care of by libraries dependencies
include_directories(${CTINT_LIBRARY_INCLUDE_DIRS})

#source files
set (LIB_FILES
        src/measurement.ipp
        src/interaction_expansion.ipp
        src/legendre.h
        src/fastupdate_formula.h
        src/update_statistics.h
        src/postprocess.hpp
        src/measurement.ipp
        )

#executable
add_executable(ctint_real src/main_real.cpp ${LIB_FILES})
target_link_libraries(ctint_real ${ALPSCore_LIBRARIES} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${EXTRA_LIBS})

add_executable(ctint_complex src/main_complex.cpp ${LIB_FILES})
target_link_libraries(ctint_complex ${ALPSCore_LIBRARIES} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${EXTRA_LIBS})

install (TARGETS ctint_real RUNTIME DESTINATION bin)
install (TARGETS ctint_complex RUNTIME DESTINATION bin)


#testing setup
option(Testing "Enable testing" ON)
include_directories("test")

function(add_gtest test link_lib)
    if (TestXMLOutput)
        set (test_xml_output --gtest_output=xml:${test}.xml)
    endif(TestXMLOutput)

    set(source "test/${test}.cpp")
    set(gtest_src "test/gtest_main.cc;test/gtest-all.cc")

    add_executable(${test} ${source} ${gtest_src})
    target_link_libraries(${test} ${link_lib})
    add_test(NAME ${test} COMMAND ${test} ${test_xml_output})
endfunction(add_gtest)

set(LINK_ALL alpscore_ctint ${ALPSCore_LIBRARIES} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${EXTRA_LIBS})
set(LINK_REST ${ALPSCore_LIBRARIES} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${EXTRA_LIBS})
enable_testing(test)

add_gtest(fastupdate "${LINK_REST}")
add_gtest(misc "${LINK_REST}")
