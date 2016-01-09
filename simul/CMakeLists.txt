cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/usr/share/cmake/")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")

set(Igraph_DIR "~/dipl/cmplx/common/")
find_package(Igraph REQUIRED)
include_directories(${IGRAPH_INCLUDE_DIR})
include_directories(${IGRAPH_LIBRARY})

add_library(simulator simulator.cc)
target_include_directories(simulator PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${IGRAPH_INCLUDE_DIR})

