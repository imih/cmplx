cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/usr/share/cmake/")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")

set(Igraph_DIR "~/dipl/cmplx/common/")
find_package(Igraph REQUIRED)
include_directories(${IGRAPH_INCLUDE_DIR})
include_directories(${IGRAPH_LIBRARY})

project (cmplx)
add_subdirectory(common)
add_subdirectory(simul)
add_library(source_detection common simul source_detector.cc)
target_include_directories(source_detection PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}) 

