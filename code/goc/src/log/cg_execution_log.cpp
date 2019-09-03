//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/log/cg_execution_log.h"

#include "goc/string/string_utils.h"

using namespace std;
using namespace nlohmann;

namespace goc
{
json CGExecutionLog::ToJSON() const
{
	json j;
	j["log_type"] = "cg"; // ID of the log type.
	if (screen_output.IsSet()) j["screen_output"] = screen_output.Value();
	if (time.IsSet()) j["time"] = time.Value().Amount(DurationUnit::Seconds);
	if (status.IsSet()) j["status"] = STR(status.Value());
	if (incumbent.IsSet()) j["incumbent"] = incumbent.Value();
	if (incumbent_value.IsSet()) j["incumbent_value"] = incumbent_value.Value();
	if (columns_added.IsSet()) j["columns_added"] = columns_added.Value();
	if (iteration_count.IsSet()) j["iteration_count"] = iteration_count.Value();
	if (pricing_time.IsSet()) j["pricing_time"] = pricing_time.Value();
	if (lp_time.IsSet()) j["lp_time"] = lp_time.Value();
	if (iterations.IsSet())
	{
		j["iterations"] = vector<json>();
		for (auto& iteration: *iterations) j["iterations"].push_back(iteration);
	}
	
	return j;
}

ostream& operator<<(ostream& os, CGStatus status)
{
	unordered_map<CGStatus, string> mapper = {{CGStatus::DidNotStart, "DidNotStart"},
											  {CGStatus::Infeasible, "Infeasible"},
											  {CGStatus::Unbounded, "Unbounded"},
											  {CGStatus::TimeLimitReached, "TimeLimitReached"},
											  {CGStatus::MemoryLimitReached, "MemoryLimitReached"},
											  {CGStatus::Optimum, "Optimum"}};
	return os << mapper[status];
}
} // namespace goc