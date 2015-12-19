cmake_minimum_required(VERSION 2.8) 

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")
add_library(common bit_array.cc)

target_include_directories (common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_BOTH_LIBRARIES})
add_test(AllTestsInBitarraytest bitarraytest)
add_executable(bitarraytest bit_array_test.cc bit_array.cc)
target_link_libraries(bitarraytest ${GTEST_BOTH_LIBRARIES})

add_executable(ba_vs_bitset bit_array_bitset.cc bit_array.cc)

