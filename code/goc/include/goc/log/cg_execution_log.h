//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_LOG_CG_EXECUTION_LOG_H
#define GOC_LOG_CG_EXECUTION_LOG_H

#include <iostream>
#include <string>
#include <vector>

#include "goc/base/maybe.h"
#include "goc/lib/json.hpp"
#include "goc/linear_programming/model/valuation.h"
#include "goc/log/log.h"
#include "goc/time/duration.h"

namespace goc
{
// All the status that can result from a column generation execution.
enum class CGStatus { DidNotStart, Infeasible, Unbounded, TimeLimitReached, MemoryLimitReached, Optimum };

// This class stores information about the execution of a column generation algorithm.
// It is compatible with the Kaleidoscope kd_type "cg".
class CGExecutionLog : public Log
{
public:
	Maybe<std::string> screen_output; // output of the algorithm in the screen.
	Maybe<Duration> time; // total time spent solving.
	Maybe<CGStatus> status; // the status of the execution
	Maybe<Valuation> incumbent; // best solution found.
	Maybe<double> incumbent_value; // value of the best solution found.
	Maybe<int> columns_added; // total number of columns added in the colgen.
	Maybe<int> iteration_count; // number of pricing iterations solved.
	Maybe<Duration> pricing_time; // time spent solving the pricing problem.
	Maybe<Duration> lp_time; // time spent solving the lp relaxation.
	Maybe<std::vector<nlohmann::json>> iterations; // logs of the pricing iterations.
	
	CGExecutionLog() = default;
	
	// Serialize log.
	virtual nlohmann::json ToJSON() const;
};

std::ostream& operator<<(std::ostream& os, CGStatus status);
} // namespace goc

#endif //GOC_LOG_CG_EXECUTION_LOG_H
