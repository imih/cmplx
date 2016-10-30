cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/usr/share/cmake/")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")

set(Igraph_DIR "~/dipl/cmplx/common/")
find_package(Igraph REQUIRED)

include_directories(${IGRAPH_INCLUDE_DIR})

add_library(common igraph.cc ivector.cc bit_array.cc idqueue.cc realization.cc)
target_include_directories(common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${IGRAPH_INCLUDE_DIR})

