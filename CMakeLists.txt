cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE Release)

project (cmplx)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/usr/share/cmake/")

find_package(MPI REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")

set(Igraph_DIR "~/dipl/cmplx/common/")
find_package(Igraph REQUIRED)
include_directories(${IGRAPH_INCLUDE_DIR})

add_subdirectory(common)
add_subdirectory(simul)
add_library(source_detection source_detector.cc seq_sample.cc)
target_include_directories(source_detection PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${MPI_LIBRARIES})

target_link_libraries(source_detection common simulator ${IGRAPH_LIBRARY})

add_executable(ss_seq source_detection_seq.cc source_detection_params.cc)
target_link_libraries(ss_seq common simulator source_detection ${IGRAPH_LIBRARY})

set(CMAKE_CXX_COMPILER mpicxx)
add_executable(main main.cc source_detection_paral.cc source_detection_params.cc)
add_executable(main_iss main_iss.cc source_detection_paral.cc source_detection_params.cc)

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS}  ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(main common simulator source_detection mpi ${IGRAPH_LIBRARY})
target_link_libraries(main_iss common simulator source_detection mpi ${IGRAPH_LIBRARY})

