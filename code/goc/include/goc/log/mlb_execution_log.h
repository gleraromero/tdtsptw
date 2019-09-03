//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_LOG_MLB_EXECUTION_LOG_H
#define GOC_LOG_MLB_EXECUTION_LOG_H

#include <iostream>
#include <string>
#include <vector>

#include "goc/base/maybe.h"
#include "goc/lib/json.hpp"
#include "goc/log/log.h"
#include "goc/time/duration.h"

namespace goc
{
// All the status that can result from a monodirectional labeling algorithm.
enum class MLBStatus { DidNotStart, TimeLimitReached, ProcessLimitReached, Finished };

// This class stores information about the execution of a MONODIRECTIONAL labeling algorithm.
// We consider that a label in the labeling algorithm has to pass through the following stages:
//	1 - enumeration: label is enumerated and put in the queue.
//	2 - extension: label was actually extended (notice with lazy extension this is different than 1).
//	3 - domination: check if a label is dominated.
//	4 - process: a non dominated label should be added to the structure.
// It is compatible with the Kaleidoscope kd_type "mlb".
class MLBExecutionLog : public Log
{
public:
	Maybe<std::string> screen_output; // output of the algorithm in the screen.
	Maybe<Duration> time; // total time spent solving the problem.
	Maybe<MLBStatus> status; // the status of the execution.
	Maybe<int> enumerated_count; // number of labels enumerated.
	Maybe<int> extended_count; // number of labels extended.
	Maybe<int> dominated_count; // number of labels dominated.
	Maybe<int> corrected_count; // number of labels corrected.
	Maybe<int> processed_count; // number of labels processed.
	Maybe<std::vector<int>> count_by_length; // count_by_length[i] indicates how many labels of length i were processed.
	Maybe<Duration> queuing_time; // time pushing and popping from the queue.
	Maybe<Duration> enumeration_time; // time spent in the enumeration phase.
	Maybe<Duration> extension_time; // time spent in the extension phase.
	Maybe<Duration> domination_time; // time spent in the domination phase.
	Maybe<Duration> correction_time; // time spent in the correction phase.
	Maybe<Duration> process_time; // time spent in the process phase.
	Maybe<Duration> positive_domination_time; // time spent in the domination phase (when the result was DOMINATED).
	Maybe<Duration> negative_domination_time; // time spent in the domination phase (when the result was NOT DOMINATED).
	
	// init_defaults: if true, then all properties are initialized with their default constructor.
	MLBExecutionLog(bool init_defaults=false);
	
	// Serialize log.
	virtual nlohmann::json ToJSON() const;
};

std::ostream& operator<<(std::ostream& os, MLBStatus status);
} // namespace goc

#endif //GOC_LOG_MLB_EXECUTION_LOG_H
