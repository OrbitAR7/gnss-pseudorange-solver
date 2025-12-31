cmake_minimum_required(VERSION 3.14)
project(gnss_pseudorange_solver 
    VERSION 1.0.0
    DESCRIPTION "C++ GNSS Pseudorange Position Solver"
    LANGUAGES CXX)

# C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
elseif(MSVC)
    add_compile_options(/W4)
endif()

# Header-only library interface
add_library(gnss_solver INTERFACE)
target_include_directories(gnss_solver INTERFACE 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Examples
option(BUILD_EXAMPLES "Build example programs" ON)
if(BUILD_EXAMPLES)
    add_executable(gps_position_fix examples/gps_position_fix.cpp)
    target_link_libraries(gps_position_fix PRIVATE gnss_solver)

    add_executable(dop_analysis examples/dop_analysis.cpp)
    target_link_libraries(dop_analysis PRIVATE gnss_solver)
endif()

# Tests
option(BUILD_TESTS "Build tests" ON)
if(BUILD_TESTS)
    enable_testing()
    add_executable(test_gnss_solver tests/test_gnss_solver.cpp)
    target_link_libraries(test_gnss_solver PRIVATE gnss_solver)
    add_test(NAME GNSSSolverTests COMMAND test_gnss_solver)
endif()

# Installation
include(GNUInstallDirs)
install(TARGETS gnss_solver
    EXPORT gnss_solver_targets
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(FILES include/gnss_solver.hpp 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
