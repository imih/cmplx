cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/usr/share/cmake/")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")

add_library(source_detector seq_sample.cc source_detector.cc source_detection_params.cc common_paral.cc paral_softmc.cc)
link_directories(common/ simul/)
target_link_libraries(source_detector common simulator)

