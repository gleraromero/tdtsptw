//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_LOG_BC_EXECUTION_LOG_H
#define GOC_LOG_BC_EXECUTION_LOG_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "goc/base/maybe.h"
#include "goc/lib/json.hpp"
#include "goc/linear_programming/model/valuation.h"
#include "goc/log/log.h"
#include "goc/time/duration.h"

namespace goc
{
// All the status that can result from a branch and cut execution.
enum class BCStatus {
	DidNotStart, Infeasible, Unbounded, TimeLimitReached, MemoryLimitReached, Optimum, NodeLimitReached
};

// This class stores information about the execution of a branch-and-cut solver.
// - It is JSON serializable and compatible with the Kaleidoscope kd_type "bc".
class BCExecutionLog : public Log
{
public:
	Maybe<std::string> screen_output; // output of the algorithm in the screen.
	Maybe<Duration> time; // total time spent solving the problem.
	Maybe<BCStatus> status; // the status of the execution.
	Maybe<int> constraint_count; // number of constraints of the initial formulation.
	Maybe<int> variable_count; // number of variables of the initial formulation.
	Maybe<int> nodes_open; // number of open nodes in the BB tree at the end of the execution.
	Maybe<int> nodes_closed; // number of nodes closed during the BB.
	Maybe<double> root_lp_value; // objective value of the root node relaxation (after cuts).
	Maybe<double> root_int_value; // value of the integer solution found at the root node.
	Maybe<Valuation> root_int_solution; // integer solution found at the root node.
	Maybe<double> best_bound; // best dual bound found for the problem.
	Maybe<Valuation> best_int_solution; // best integer solution found.
	Maybe<double> best_int_value; // value of the best integer solution found.
	Maybe<int> cut_count; // total number of cuts added.
	Maybe<Duration> cut_time; // total time spent separating cuts.
	Maybe<std::vector<std::string>> cut_families; // vector of families considered in the BB algorithm.
	Maybe<std::unordered_map<std::string, int>> cut_family_cut_count; // number of cuts added per family.
	Maybe<std::unordered_map<std::string, int>> cut_family_iteration_count; // number of cut iterations done per family.
	Maybe<std::unordered_map<std::string, Duration>> cut_family_cut_time; // time spent separating cuts per family.
	
	BCExecutionLog() = default;
	
	virtual nlohmann::json ToJSON() const;
};

std::ostream& operator<<(std::ostream& os, BCStatus status);
} // namespace goc

#endif //GOC_LOG_BC_EXECUTION_LOG_H
