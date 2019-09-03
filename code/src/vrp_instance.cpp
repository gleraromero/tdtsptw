//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "vrp_instance.h"

using namespace std;
using namespace goc;
using namespace nlohmann;

namespace tdtsptw
{
TimeUnit VRPInstance::ReadyTime(const GraphPath& p, TimeUnit t0) const
{
	TimeUnit t = t0;
	for (int k = 0; k < (int)p.size()-1; ++k)
	{
		Vertex i = p[k], j = p[k+1];
		if (!tau[i][j].Domain().Includes(t)) return INFTY;
		t += tau[i][j](t);
	}
	return t;
}

void VRPInstance::Print(ostream& os) const
{
	os << json(*this);
}

void to_json(json& j, const VRPInstance& instance)
{
	j["digraph"] = instance.D;
	j["start_depot"] = instance.o;
	j["end_depot"] = instance.d;
	j["horizon"] = vector<TimeUnit>({0, instance.T});
	j["time_windows"] = instance.tw;
	j["travel_times"] = instance.tau;
}

void from_json(const json& j, VRPInstance& instance)
{
	int n = j["digraph"]["vertex_count"];
	instance.D = j["digraph"];
	instance.o = j["start_depot"];
	instance.d = j["end_depot"];
	instance.T = j["horizon"][1];
	for (int i = 0; i < n; ++i) instance.tw.push_back(j["time_windows"][i]);
	instance.tau = j["travel_times"];
}
} // namespace tdtsptw