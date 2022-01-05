cmake_minimum_required(VERSION 3.0)
project(quadprog VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)

file(GLOB_RECURSE SRC_FILES src/*.cpp)

add_executable (quadprog ${SRC_FILES})

target_include_directories(quadprog PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(quadprog PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen)
