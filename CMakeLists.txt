cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE Release)

project (cmplx)

option(test "Build all tests." OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/usr/share/cmake/")

find_package(MPI REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")

set(Igraph_DIR "~/dipl/cmplx/common/")
find_package(Igraph REQUIRED)
include_directories(${IGRAPH_INCLUDE_DIR})

add_subdirectory(common)
add_subdirectory(simul)
add_subdirectory(source_detector)

#add_library(source_detection source_detector.cc seq_sample.cc)
#target_include_directories(source_detection PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${MPI_LIBRARIES})

#target_link_libraries(source_detection common simulator ${IGRAPH_LIBRARY})

add_executable(ss_seq source_detection_seq.cc source_detection_params.cc)
target_link_libraries(ss_seq common simulator source_detector ${IGRAPH_LIBRARY})

set(CMAKE_CXX_COMPILER mpicxx)
add_executable(main main.cc source_detection_paral.cc source_detection_params.cc)
add_executable(main_iss main_iss.cc source_detection_paral.cc source_detection_params.cc)

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS}  ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(main common simulator source_detector mpi ${IGRAPH_LIBRARY})
target_link_libraries(main_iss common simulator source_detector mpi ${IGRAPH_LIBRARY})

if(test)
  enable_testing()
  find_package(GTest REQUIRED)
  include_directories(${GTEST_BOTH_LIBRARIES} )
  set(EXECUTABLE_OUTPUT_PATH "test")

  #BitArray tests and comp
  add_executable(bitarraytest common/bit_array_test.cc common/bit_array.cc)
  target_link_libraries(bitarraytest ${GTEST_BOTH_LIBRARIES})
  add_test(AllTestsInBitarraytest test/bitarraytest)
  #add_executable(ba_vs_bitset bit_array_bitset.cc bit_array.cc)
  
  #IVector test
  add_executable(ivectortest common/ivector.cc common/ivector_test.cc)
  target_link_libraries(ivectortest ${GTEST_BOTH_LIBRARIES} ${IGRAPH_LIBRARY})
  add_test(AllTestInIvectortest test/ivectortest)

  #IGraph test
  add_executable(igraphtest common/igraph.cc common/igraph_test.cc common/ivector.cc)
  target_link_libraries(igraphtest ${GTEST_BOTH_LIBRARIES} ${IGRAPH_LIBRARY})
  add_test(AllTestInIgraphtest test/igraphtest)

  #Simulator test
  add_executable(simulatortest simul/simulator_test.cc
   simul/simulator.cc common/bit_array.cc common/igraph.cc common/realization.cc common/idqueue.cc)
  target_link_libraries(simulatortest ${GTEST_BOTH_LIBRARIES} ${IGRAPH_LIBRARY})
  add_test(AllTestInSimulatortest test/simulatortest)
endif()

