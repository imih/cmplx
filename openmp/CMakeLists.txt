cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/usr/share/cmake/")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")


find_package(OpenMP REQUIRED) 
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g ${OpenMP_CXX_FLAGS}")
endif()

add_library(openmp_detection omp_directmc.cc omp_softmc.cc omp_seqis.cc)
target_include_directories(openmp_detection PUBLIC)

link_directories(source_detector/)
target_link_libraries(openmp_detection source_detector)

