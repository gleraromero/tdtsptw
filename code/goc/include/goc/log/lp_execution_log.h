//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_LOG_LP_EXECUTION_LOG_H
#define GOC_LOG_LP_EXECUTION_LOG_H

#include <iostream>
#include <string>
#include <vector>

#include "goc/base/maybe.h"
#include "goc/linear_programming/model/valuation.h"
#include "goc/log/log.h"
#include "goc/time/duration.h"

namespace goc
{
// All the status that can result from a simplex execution.
enum class LPStatus { DidNotStart, Infeasible, Unbounded, TimeLimitReached, MemoryLimitReached, Optimum };

// This class stores information about the execution of a LP relaxation solver.
// It is compatible with the Kaleidoscope kd_type "lp".
class LPExecutionLog : public Log
{
public:
	Maybe<std::string> screen_output; // output of the algorithm in the screen.
	Maybe<Duration> time; // total time spent solving the relaxation.
	Maybe<LPStatus> status; // the status of the execution
	Maybe<int> simplex_iterations; // number of iterations the simplex algorithm did.
	Maybe<Valuation> incumbent; // best solution found.
	Maybe<double> incumbent_value; // value of the best solution found.
	Maybe<int> variable_count; // number of variables in the lp.
	Maybe<int> constraint_count; // number of constraints in the lp.
	Maybe<std::vector<double>> duals; // vector of the dual variables associated to the rows in the solution.
	
	LPExecutionLog() = default;
	
	// Serialize log.
	virtual nlohmann::json ToJSON() const;
};

std::ostream& operator<<(std::ostream& os, LPStatus status);
} // namespace goc

#endif //GOC_LOG_LP_EXECUTION_LOG_H
