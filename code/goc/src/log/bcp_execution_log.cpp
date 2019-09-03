//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/log/bcp_execution_log.h"

using namespace std;
using namespace nlohmann;

namespace goc
{
json BCPExecutionLog::ToJSON() const
{
	json j = BCExecutionLog::ToJSON();
	j["kd_type"] = "bcp";
	if (root_constraint_count.IsSet()) j["root_constraint_count"] = root_constraint_count.Value();
	if (root_variable_count.IsSet()) j["root_variable_count"] = root_variable_count.Value();
	if (final_constraint_count.IsSet()) j["final_constraint_count"] = final_constraint_count.Value();
	if (final_variable_count.IsSet()) j["final_variable_count"] = final_variable_count.Value();
	if (root_log.IsSet()) j["root_log"] = root_log.Value();
	if (lp_time.IsSet()) j["lp_time"] = lp_time.Value();
	if (pricing_time.IsSet()) j["pricing_time"] = pricing_time.Value();
	if (branching_time.IsSet()) j["branching_time"] = branching_time.Value();
	return j;
}
} // namespace goc