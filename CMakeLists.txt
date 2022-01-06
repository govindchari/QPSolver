cmake_minimum_required(VERSION 3.0)
project(QPSolver VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)

file(GLOB_RECURSE SRC_FILES src/*.cpp)

add_executable (QPSolver ${SRC_FILES})

target_include_directories(QPSolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(QPSolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen)
