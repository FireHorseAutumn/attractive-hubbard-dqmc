cmake_minimum_required(VERSION 3.20)
project(attractive_hubbard_dqmc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Eigen3 REQUIRED NO_MODULE)
enable_testing()

add_library(dqmc_core src/dqmc.cpp src/current.cpp)
target_include_directories(dqmc_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(dqmc_core PUBLIC Eigen3::Eigen)
target_compile_options(dqmc_core PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_core tests/test_core.cpp)
target_link_libraries(test_core PRIVATE dqmc_core)
target_compile_options(test_core PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME core COMMAND test_core)

add_executable(test_current tests/test_current.cpp)
target_link_libraries(test_current PRIVATE dqmc_core)
target_compile_options(test_current PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME current COMMAND test_current)
