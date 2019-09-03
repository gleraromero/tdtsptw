//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_LOG_BLB_EXECUTION_LOG_H
#define GOC_LOG_BLB_EXECUTION_LOG_H

#include <iostream>
#include <string>

#include "goc/base/maybe.h"
#include "goc/lib/json.hpp"
#include "goc/log/mlb_execution_log.h"
#include "goc/log/log.h"
#include "goc/time/duration.h"

namespace goc
{
// All the status that can result from a monodirectional labeling algorithm.
enum class BLBStatus { DidNotStart, TimeLimitReached, SolutionLimitReached, Finished };

// This class stores information about the execution of a BIDIRECTIONAL labeling algorithm.
// It is compatible with the Kaleidoscope kd_type "blb".
class BLBExecutionLog : public Log
{
public:
	Maybe<BLBStatus> status; // status at the end of the execution.
	Maybe<Duration> time; // total execution time.
	Maybe<std::string> screen_output; // output of the algorithm in the screen.
	Maybe<MLBExecutionLog> forward_log; // log of the forward labeling.
	Maybe<MLBExecutionLog> backward_log; // log of the backward labeling.
	Maybe<Duration> merge_time; // time spent merging labels.
	
	// init_defaults: if true, then all properties are initialized with their default constructor.
	BLBExecutionLog(bool init_defaults=false);
	
	// Serialize log.
	virtual nlohmann::json ToJSON() const;
};

std::ostream& operator<<(std::ostream& os, BLBStatus status);
} // namespace goc

#endif //GOC_LOG_BLB_EXECUTION_LOG_H
