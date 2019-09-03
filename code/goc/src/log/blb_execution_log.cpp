//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/log/blb_execution_log.h"

#include <unordered_map>

#include "goc/string/string_utils.h"

using namespace std;
using namespace nlohmann;

namespace goc
{
BLBExecutionLog::BLBExecutionLog(bool init_defaults)
{
	if (init_defaults)
	{
		status = BLBStatus::DidNotStart;
		time = merge_time = Duration();
		screen_output = "";
		forward_log = backward_log = MLBExecutionLog(true);
	}
}

json BLBExecutionLog::ToJSON() const
{
	json j;
	j["kd_type"] = "blb"; // ID of the log type.
	if (screen_output.IsSet()) j["screen_output"] = screen_output.Value();
	if (time.IsSet()) j["time"] = time.Value();
	if (status.IsSet()) j["status"] = STR(status.Value());
	if (forward_log.IsSet()) j["forward"] = forward_log.Value();
	if (backward_log.IsSet()) j["backward"] = backward_log.Value();
	if (merge_time.IsSet()) j["merge_time"] = merge_time.Value();
	
	return j;
}

ostream& operator<<(ostream& os, BLBStatus status)
{
	unordered_map<BLBStatus, string> mapper = {{BLBStatus::DidNotStart, "DidNotStart"},
											   {BLBStatus::TimeLimitReached, "TimeLimitReached"},
											   {BLBStatus::SolutionLimitReached, "SolutionLimitReached"},
											   {BLBStatus::Finished, "Finished"}};
	return os << mapper[status];
}
} // namespace goc