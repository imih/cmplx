cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/usr/share/cmake/")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")

find_package(MPI REQUIRED)

set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS}  ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(${MPI_INCLUDE_PATH})

add_library(mpi_detection mpi_paral.cc directmc_paral.cc softmc_paral.cc seqis_paral.cc)
target_include_directories(mpi_detection PUBLIC mpi)

link_directories(source_detector/)
target_link_libraries(mpi_detection mpi source_detector)

