//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/log/lp_execution_log.h"

#include <unordered_map>

#include "goc/string/string_utils.h"

using namespace std;
using namespace nlohmann;

namespace goc
{

json LPExecutionLog::ToJSON() const
{
	json j;
	j["kd_type"] = "lp"; // ID of the log type.
	if (screen_output.IsSet()) j["screen_output"] = screen_output.Value();
	if (time.IsSet()) j["time"] = time.Value().Amount(DurationUnit::Seconds);
	if (status.IsSet()) j["status"] = STR(status.Value());
	if (simplex_iterations.IsSet()) j["simplex_iterations"] = simplex_iterations.Value();
	if (incumbent.IsSet()) j["incumbent"] = incumbent.Value();
	if (incumbent_value.IsSet()) j["incumbent_value"] = incumbent_value.Value();
	if (constraint_count.IsSet()) j["constraint_count"] = constraint_count.Value();
	if (variable_count.IsSet()) j["variable_count"] = variable_count.Value();
	if (duals.IsSet()) j["duals"] = duals.Value();
	return j;
}

ostream& operator<<(ostream& os, LPStatus status)
{
	unordered_map<LPStatus, string> mapper = {{LPStatus::DidNotStart, "DidNotStart"},
											  {LPStatus::Infeasible, "Infeasible"},
											  {LPStatus::Unbounded, "Unbounded"},
											  {LPStatus::TimeLimitReached, "TimeLimitReached"},
											  {LPStatus::MemoryLimitReached, "MemoryLimitReached"},
											  {LPStatus::Optimum, "Optimum"}};
	return os << mapper[status];
}
} // namespace goc