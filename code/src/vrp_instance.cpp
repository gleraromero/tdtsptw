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
TimeUnit VRPInstance::TravelTime(Arc e, TimeUnit t0) const
{
	double t = ArrivalTime(e, t0);
	if (t == INFTY) return INFTY;
	return t - t0;
}

TimeUnit VRPInstance::PreTravelTime(Arc e, TimeUnit tf) const
{
	double t = DepartureTime(e, tf);
	if (t == INFTY) return INFTY;
	return tf - t;
}

TimeUnit VRPInstance::ArrivalTime(Arc e, TimeUnit t0) const
{
	auto& tau_e = tau[e.tail][e.head];
	if (epsilon_smaller(tau_e.Domain().right, t0)) return INFTY;
	else if (epsilon_bigger(tau_e.Domain().left, t0)) return tau_e.Domain().left+tau_e(tau_e.Domain().left);
	return t0 + tau_e(t0);
}

TimeUnit VRPInstance::DepartureTime(Arc e, TimeUnit tf) const
{
	auto& dep_e = dep[e.tail][e.head];
	if (epsilon_smaller(tf, dep_e.Domain().left)) return INFTY;
	if (epsilon_bigger(tf, dep_e.Domain().right)) return max(img(dep_e));
	return dep[e.tail][e.head](tf);
}

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

TimeUnit VRPInstance::PathDepartureTime(const GraphPath& p, TimeUnit tf) const
{
	TimeUnit t = tf;
	for (int k = (int)p.size()-1; k > 0; --k)
	{
		Vertex i = p[k-1], j = p[k];
		t = DepartureTime({i,j}, t);
		if (t == INFTY) return INFTY;
	}
	return t;
}

TimeUnit VRPInstance::MinimumTravelTime(Arc e, TimeUnit t0, TimeUnit tf) const
{
	TimeUnit tmin = INFTY;
	int j = 0;
	int v = e.tail, w = e.head;
	while (j < tau[v][w].PieceCount())
	{
		if (epsilon_bigger(tau[v][w][j].domain.left, tf)) break;
		if (tau[v][w][j].domain.Intersects({t0, tf}))
			tmin = min(tmin, tau[v][w][j].image.left);
		if (epsilon_bigger_equal(tau[v][w][j].domain.right,tf)) break;
		++j;
	}
	return tmin;
}

Route VRPInstance::BestDurationRoute(const GraphPath& p) const
{
	PWLFunction Delta = arr[p[0]][p[0]];
	if (Delta.Empty()) return {{}, 0.0, INFTY};
	for (int k = 0; k < (int)p.size()-1; ++k)
	{
		Vertex i = p[k], j = p[k+1];
		Delta = arr[i][j].Compose(Delta);
		if (Delta.Empty()) return {{}, 0.0, INFTY};
	}
	Delta = Delta - PWLFunction::IdentityFunction(dom(Delta));
	return Route(p, Delta.PreValue(min(img(Delta))), min(img(Delta)));
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
	for (int i = 0; i < n; ++i) instance.a.push_back(instance.tw[i].left);
	for (int i = 0; i < n; ++i) instance.b.push_back(instance.tw[i].right);
	if (has_key(j, "EAT")) instance.EAT = j["EAT"];
	if (has_key(j, "LDT")) instance.LDT = j["LDT"];
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
	if (has_key(j, "precedence_matrix")) {
		instance.prec = Matrix<bool>(n, n);
		for (Vertex i: instance.D.Vertices())
			for (Vertex k: instance.D.Vertices())
				instance.prec[i][k] = (bool) j["precedence_matrix"][i][k];

		instance.prec_count = vector<int>(n, 0);
		instance.suc_count = vector<int>(n, 0);
		for (Vertex i: instance.D.Vertices()) {
			for (Vertex k: instance.D.Vertices()) {
				if (instance.prec[i][k]) {
					instance.prec_count[k]++;
					instance.suc_count[i]++;
				}
			}
		}

		instance.time_prec = Matrix<pair<double, Vertex> >(n, n+1);
		instance.time_succ = Matrix<pair<double, Vertex> >(n, n+1);

		for (int k = 0; k <= n; k++)
		{
			instance.time_prec[instance.o][k] = pair<double, Vertex>(instance.tw[instance.o].left, -1);
			instance.time_prec[instance.d][k] = pair<double, Vertex>(instance.tw[instance.d].right, -1);

			instance.time_succ[instance.o][k] = pair<double, Vertex>(instance.tw[instance.o].left, -1);
			instance.time_succ[instance.d][k] = pair<double, Vertex>(instance.tw[instance.d].right, -1);
		}

		// Compute LDT times for the predecesors.
		for (Vertex i : instance.D.Vertices())
		{
		    if (i == instance.o || i == instance.d) continue;

			vector<pair<double, Vertex> > limit_times;
			for (Vertex k : instance.D.Vertices())
			{
				if (k == i || k == instance.o || k == instance.d) continue;

				limit_times.push_back(pair<double, Vertex>(instance.LDT(i,k), k));
			}
			sort(limit_times.begin(), limit_times.end());

			instance.time_prec[i][0] = pair<double, Vertex>(instance.tw[i].left, -1);
			instance.time_prec[i][1] = pair<double, Vertex>(instance.tw[i].left, -1);
			instance.time_prec[i][2] = pair<double, Vertex>(instance.tw[i].left, instance.o);

			for (int k = 3; k < n; k++)
			{
				instance.time_prec[i][k] = pair<double, Vertex>(max(limit_times[k-3].first, instance.tw[i].left), limit_times[k-3].second);
			}

			instance.time_prec[i][n] = pair<double, Vertex>(instance.tw[i].right, instance.d);
		}

		// Compute EAT times for the successors
		for (Vertex i : instance.D.Vertices())
		{
			if (i == instance.o || i == instance.d) continue;
			vector<pair<double, Vertex> > limit_times;
			for (Vertex k : instance.D.Vertices())
			{
				if (k == i || k == instance.o || k == instance.d) continue;

				limit_times.push_back(pair<double, Vertex>(instance.EAT(k,i), k));
			}

			sort(limit_times.begin(), limit_times.end());

			instance.time_succ[i][0] = pair<double, Vertex>(instance.tw[i].right, -1);
			instance.time_succ[i][1] = pair<double, Vertex>(instance.tw[i].right, -1);
			instance.time_succ[i][2] = pair<double, Vertex>(instance.tw[i].right, instance.d);

            for (int k = 3; k < n; k++)
			{
				instance.time_succ[i][n - 1 - k + 3] = pair<double, Vertex>(min(limit_times[k-3].first, instance.tw[i].right), limit_times[k-3].second);
			}

			instance.time_succ[i][n] = pair<double, Vertex>(instance.tw[i].left, instance.o);
		}

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j <= n; j++)
			{
				clog << instance.time_succ[i][j] << " ";
			}
			clog << endl;
		}
		clog << n << endl;
		clog << instance.tw	<< endl;
		exit(1);
	}
}
} // namespace tdtsptw