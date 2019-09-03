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

TimeUnit VRPInstance::ArrivalTime(Arc e, TimeUnit t0) const
{
	auto& tau_e = tau[e.tail][e.head];
	if (epsilon_smaller(tau_e.Domain().right, t0)) return INFTY;
	else if (epsilon_bigger(tau_e.Domain().left, t0)) return tau_e.Domain().left+tau_e(tau_e.Domain().left);
	return t0 + tau_e(t0);
}

TimeUnit VRPInstance::TravelTime(Arc e, TimeUnit t0) const
{
	double t = ArrivalTime(e, t0);
	if (t == INFTY) return INFTY;
	return t - t0;
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
	instance.tw = vector<Interval>(j["time_windows"].begin(), j["time_windows"].end());
	instance.EAT = j["EAT"];
	instance.LDT = j["LDT"];
	instance.tau = instance.arr = instance.dep = instance.pretau = Matrix<PWLFunction>(n, n);
	for (Vertex u: instance.D.Vertices())
	{
		for (Vertex v: instance.D.Successors(u))
		{
			instance.tau[u][v] = j["travel_times"][u][v];
			instance.arr[u][v] = instance.tau[u][v] + PWLFunction::IdentityFunction(instance.tau[u][v].Domain());
			instance.dep[u][v] = instance.arr[u][v].Inverse();
			instance.pretau[u][v] = PWLFunction::IdentityFunction(instance.dep[u][v].Domain()) - instance.dep[u][v];
		}
	}
	// Add travel functions for (i, i) (for boundary reasons).
	for (Vertex u: instance.D.Vertices())
	{
		instance.tau[u][u] = instance.pretau[u][u] = PWLFunction::ConstantFunction(0.0, instance.tw[u]);
		instance.dep[u][u] = instance.arr[u][u] = PWLFunction::IdentityFunction(instance.tw[u]);
	}
}
} // namespace tdtsptw