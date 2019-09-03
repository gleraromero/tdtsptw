//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/log/bc_execution_log.h"

#include "goc/collection/collection_utils.h"
#include "goc/string/string_utils.h"

using namespace std;
using namespace nlohmann;

namespace goc
{
json BCExecutionLog::ToJSON() const
{
	json j;
	j["kd_type"] = "bc";
	if (screen_output.IsSet()) j["screen_output"] = screen_output.Value();
	if (time.IsSet()) j["time"] = time.Value().Amount(DurationUnit::Seconds);
	if (status.IsSet()) j["status"] = STR(status.Value());
	if (constraint_count.IsSet()) j["constraint_count"] = constraint_count.Value();
	if (variable_count.IsSet()) j["variable_count"] = variable_count.Value();
	if (nodes_open.IsSet()) j["nodes_open"] = nodes_open.Value();
	if (nodes_closed.IsSet()) j["nodes_closed"] = nodes_closed.Value();
	if (root_lp_value.IsSet()) j["root_lp_value"] = root_lp_value.Value();
	if (root_int_solution.IsSet()) j["root_int_solution"] = root_int_solution.Value();
	if (root_int_value.IsSet()) j["root_int_value"] = root_int_value.Value();
	if (best_bound.IsSet()) j["best_bound"] = best_bound.Value();
	if (best_int_solution.IsSet()) j["best_int_solution"] = best_int_solution.Value();
	if (best_int_value.IsSet()) j["best_int_value"] = best_int_value.Value();
	if (cut_count.IsSet()) j["cut_count"] = cut_count.Value();
	if (cut_time.IsSet()) j["cut_time"] = cut_time.Value();
	if (cut_families.IsSet() && !cut_families->empty())
	{
		j["cut_families"] = vector<json>();
		for (auto& family: cut_families.Value())
		{
			json cut_family_json;
			cut_family_json["name"] = family;
			if (cut_family_cut_count.IsSet() && includes_key(cut_family_cut_count.Value(), family))
				cut_family_json["cut_count"] = cut_family_cut_count.Value().at(family);
			if (cut_family_iteration_count.IsSet() && includes_key(cut_family_iteration_count.Value(), family))
				cut_family_json["cut_iterations"] = cut_family_iteration_count.Value().at(family);
			if (cut_family_cut_time.IsSet() && includes_key(cut_family_cut_time.Value(), family))
				cut_family_json["cut_time"] = cut_family_cut_time.Value().at(family);
			j["cut_families"].push_back(cut_family_json);
		}
	}
	return j;
}

ostream& operator<<(ostream& os, BCStatus status)
{
	unordered_map<BCStatus, string> mapper = {{BCStatus::DidNotStart, "DidNotStart"},
											   {BCStatus::Infeasible, "Infeasible"},
											   {BCStatus::Unbounded, "Unbounded"},
											   {BCStatus::TimeLimitReached, "TimeLimitReached"},
											   {BCStatus::MemoryLimitReached, "MemoryLimitReached"},
											   {BCStatus::Optimum, "Optimum"},
											   {BCStatus::NodeLimitReached, "NodeLimitReached"}};
	return os << mapper[status];
}
} // namespace goc