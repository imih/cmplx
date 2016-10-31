cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/usr/share/cmake/")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")

add_library(source_detector source_detector.cc seq_sample.cc)
target_include_directories(source_detector PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

