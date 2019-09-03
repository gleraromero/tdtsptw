//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/log/mlb_execution_log.h"

#include <unordered_map>

#include "goc/string/string_utils.h"

using namespace std;
using namespace nlohmann;

namespace goc
{
MLBExecutionLog::MLBExecutionLog(bool init_defaults)
{
	if (init_defaults)
	{
		screen_output = "";
		time = 0.0_hr;
		status = MLBStatus::DidNotStart;
		enumerated_count = extended_count = dominated_count = corrected_count = processed_count = 0;
		count_by_length = vector<int>();
		queuing_time = enumeration_time = extension_time = domination_time = correction_time = process_time = negative_domination_time = positive_domination_time = 0.0_hr;
	}
}

json MLBExecutionLog::ToJSON() const
{
	json j;
	j["kd_type"] = "mlb"; // ID of the log type.
	if (screen_output.IsSet()) j["screen_output"] = screen_output.Value();
	if (time.IsSet()) j["time"] = time.Value();
	if (status.IsSet()) j["status"] = STR(status.Value());
	if (enumerated_count.IsSet()) j["enumerated_count"] = enumerated_count.Value();
	if (extended_count.IsSet()) j["extended_count"] = extended_count.Value();
	if (dominated_count.IsSet()) j["dominated_count"] = dominated_count.Value();
	if (corrected_count.IsSet()) j["corrected_count"] = corrected_count.Value();
	if (processed_count.IsSet()) j["processed_count"] = processed_count.Value();
	if (count_by_length.IsSet()) j["count_by_length"] = count_by_length.Value();
	if (queuing_time.IsSet()) j["queuing_time"] = queuing_time.Value();
	if (enumeration_time.IsSet()) j["enumeration_time"] = enumeration_time.Value();
	if (extension_time.IsSet()) j["extension_time"] = extension_time.Value();
	if (domination_time.IsSet()) j["domination_time"] = domination_time.Value();
	if (correction_time.IsSet()) j["correction_time"] = correction_time.Value();
	if (process_time.IsSet()) j["process_time"] = process_time.Value();
	if (positive_domination_time.IsSet()) j["positive_domination_time"] = positive_domination_time.Value();
	if (negative_domination_time.IsSet()) j["negative_domination_time"] = negative_domination_time.Value();
	
	return j;
}

ostream& operator<<(ostream& os, MLBStatus status)
{
	unordered_map<MLBStatus, string> mapper = {{MLBStatus::DidNotStart, "DidNotStart"},
											  {MLBStatus::TimeLimitReached, "TimeLimitReached"},
											  {MLBStatus::ProcessLimitReached, "ProcessLimitReached"},
											  {MLBStatus::Finished, "Finished"}};
	return os << mapper[status];
}
} // namespace goc