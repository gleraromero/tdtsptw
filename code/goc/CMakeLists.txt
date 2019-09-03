cmake_minimum_required(VERSION 2.8.4)
project(goc)

# Check if CPLEX environment variables are set.
if("$ENV{CPLEX_INCLUDE}" STREQUAL "")
    message(SEND_ERROR "CPLEX_INCLUDE environment variable is not set. Make sure to point the variable to the CPLEX include directory before compiling.\nTo do this run on the terminal (for example):\nexport CPLEX_INCLUDE=/opt/ibm/ILOG/CPLEX_Studio129/cplex/include\nThis can be added permanently by adding the export command to the file /etc/environment and rebooting (on Linux).")
endif()
if("$ENV{CPLEX_BIN}" STREQUAL "")
    message(SEND_ERROR "CPLEX_BIN environment variable is not set. Make sure to point the variable to the CPLEX binary (.a) before compiling.\nTo do this run on the terminal (for example):\nexport CPLEX_BIN=/opt/ibm/ILOG/CPLEX_Studio129/cplex/lib/x86-64_linux/static_pic/libcplex.a\nThis can be added permanently by adding the export command to the file /etc/environment and rebooting (on Linux).")
endif()

# Check if BOOST environment variables are set.
if("$ENV{BOOST_INCLUDE}" STREQUAL "")
    message(SEND_ERROR "BOOST_INCLUDE environment variable is not set. Make sure to point the variable to the BOOST include directory before compiling.\nTo do this run on the terminal (for example):\nexport BOOST_INCLUDE=/usr/include\nThis can be added permanently by adding the export command to the file /etc/environment and rebooting (on Linux).")
endif()
if("$ENV{BOOST_BIN}" STREQUAL "")
    message(SEND_ERROR "BOOST_BIN environment variable is not set. Make sure to point the variable to the BOOST binary (.a) before compiling.\nTo do this run on the terminal (for example):\nexport BOOST_BIN=/opt/ibm/ILOG/CPLEX_Studio129/cplex/lib/x86-64_linux/static_pic/libcplex.a\nThis can be added permanently by adding the export command to the file /etc/environment and rebooting (on Linux).")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++14")

include_directories(include)
add_library(goc src/collection/collection_utils.cpp src/graph/arc.cpp src/graph/digraph.cpp src/math/interval.cpp src/math/linear_function.cpp src/linear_programming/model/variable.cpp src/linear_programming/model/expression.cpp src/linear_programming/model/constraint.cpp src/linear_programming/cplex/cplex_formulation.cpp src/linear_programming/model/valuation.cpp src/time/duration.cpp src/time/stopwatch.cpp src/time/watch.cpp src/time/date.cpp src/time/point_in_time.cpp src/print/string_utils.cpp src/runner/runner_utils.cpp src/json/json_utils.cpp src/print/printable.cpp src/linear_programming/cplex/cplex_solver.cpp src/log/lp_execution_log.cpp src/log/bcp_execution_log.cpp src/linear_programming/cplex/cplex_wrapper.cpp src/linear_programming/solver/lp_solver.cpp src/linear_programming/solver/bc_solver.cpp src/linear_programming/cuts/separation_algorithm.cpp src/log/mlb_execution_log.cpp src/log/blb_execution_log.cpp src/linear_programming/colgen/colgen.cpp src/log/cg_execution_log.cpp src/linear_programming/solver/cg_solver.cpp src/graph/path_finding.cpp src/graph/graph_path.cpp src/print/table_stream.cpp src/graph/maxflow_mincut.cpp src/linear_programming/cuts/separation_strategy.cpp src/math/pwl_function.cpp src/log/log.cpp src/log/bc_execution_log.cpp src/math/point_2d.cpp src/graph/edge.cpp src/graph/graph.cpp src/vrp/route.cpp src/vrp/vrp_solution.cpp)

include_directories($ENV{CPLEX_INCLUDE})
include_directories($ENV{BOOST_INCLUDE})
target_link_libraries(goc "$ENV{CPLEX_BIN}" -ldl -lm)
target_link_libraries(goc "$ENV{BOOST_BIN}")