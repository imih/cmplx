cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/usr/share/cmake/")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")

set(Igraph_DIR "~/dipl/cmplx/common/")
find_package(Igraph REQUIRED)
include_directories(${IGRAPH_INCLUDE_DIR})
include_directories(${IGRAPH_LIBRARY})

add_library(common igraph.cc ivector.cc bit_array.cc idqueue.cc sir_params.cc)
target_include_directories(common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${IGRAPH_INCLUDE_DIR})

enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_BOTH_LIBRARIES} )

#BitArray tests and comp
add_executable(bitarraytest bit_array_test.cc bit_array.cc)
target_link_libraries(bitarraytest ${GTEST_BOTH_LIBRARIES})
add_test(AllTestsInBitarraytest bitarraytest)

#IVector test
add_executable(ivectortest ivector.cc ivector_test.cc)
target_link_libraries(ivectortest ${GTEST_BOTH_LIBRARIES} ${IGRAPH_LIBRARY})
add_test(AllTestInIvectortest ivectortest)

#IGraph test
add_executable(igraphtest igraph.cc igraph_test.cc ivector.cc)
target_link_libraries(igraphtest ${GTEST_BOTH_LIBRARIES} ${IGRAPH_LIBRARY})
add_test(AllTestInIgraphtest igraphtest)


#Random test
add_executable(randomtest random.cc random_test.cc)
target_link_libraries(randomtest ${GTEST_BOTH_LIBRARIES})
add_test(AllTestInIgraphtest randomtest)


add_executable(ba_vs_bitset bit_array_bitset.cc bit_array.cc)
