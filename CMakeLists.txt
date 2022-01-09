cmake_minimum_required(VERSION 3.0)
project(QPSolver VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)

file(GLOB_RECURSE SRC_FILES src/*.cpp)

add_library (QPSolver ${SRC_FILES})

target_include_directories(QPSolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(QPSolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen)

include_directories(tests)
add_executable(run_tests tests/run_tests.cpp)
target_include_directories(run_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(run_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen)
target_link_libraries(run_tests QPSolver)

