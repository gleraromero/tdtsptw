cmake_minimum_required(VERSION 2.8.4)
project(tdtsptw)

# Set compiler flags.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++14 -m64")
if ("${CMAKE_BUILD_TYPE}" MATCHES release OR "${CMAKE_BUILD_TYPE}" MATCHES Release)
    message("-- Compiling in Release Mode")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -DNDEBUG")
elseif ("${CMAKE_BUILD_TYPE}" MATCHES debug OR "${CMAKE_BUILD_TYPE}" MATCHES Debug)
    message("-- Compiling in Debug Mode")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()

if ("${RUNNER}" MATCHES ON)
    message("-- Runner compilation")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DRUNNER")
endif()

# Add include directory.
include_directories("include/")

# Including GOC.
add_subdirectory(goc)
include_directories(goc/include)

# Create binaries.
add_executable(main src/main.cpp src/vrp_instance.cpp src/preprocess_travel_times.cpp src/preprocess_time_windows.cpp src/preprocess_waiting_times.cpp src/spf.cpp src/pricing_problem.cpp src/heuristic.cpp src/preprocess_remove_time_windows.cpp src/preprocess_remove_time_dependency.cpp src/ngl_info.cpp src/core.cpp src/label_sequence_ti.cpp src/label_sequence_td.cpp src/dynamic_neighbour_augmentation.cpp src/exact_solver.cpp src/bounding_tree.cpp src/preprocess_time_precedence.cpp src/relaxation_solver.cpp src/column_generation.cpp)
target_link_libraries(main goc)
target_link_libraries(main)