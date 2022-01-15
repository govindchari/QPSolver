cmake_minimum_required(VERSION 3.0)
project(QPSolver VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_FLAGS "-O3")

add_library (QPSolver src/QP.cpp)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries (QPSolver Eigen3::Eigen)

target_include_directories(QPSolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(QPSolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen)

#Testing
if (ENABLE_TESTING)
    include_directories(tests)
    add_executable(run_tests tests/run_tests.cpp)
    target_include_directories(run_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_include_directories(run_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen)
    target_include_directories(run_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/tests)

    set(ENABLE_OSQP TRUE)
    set(ENABLE_ECOS TRUE)
    add_subdirectory(tests/Epigraph)
    target_link_libraries(run_tests QPSolver epigraph)
endif()
