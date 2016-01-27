cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/usr/share/cmake/")

find_package(MPI REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")

set(Igraph_DIR "~/dipl/cmplx/common/")
find_package(Igraph REQUIRED)
include_directories(${IGRAPH_INCLUDE_DIR})
include_directories(${IGRAPH_LIBRARY})

project (cmplx)
add_subdirectory(common)
add_subdirectory(simul)
add_library(source_detection source_detector.cc)
target_include_directories(source_detection PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${MPI_LIBRARIES}) 
target_link_libraries(source_detection simulator igraph)

set(CMAKE_CXX_COMPILER mpicxx)
add_executable(ss_para single_source_paral.cc)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(ss_para common simulator source_detection mpi igraph)
