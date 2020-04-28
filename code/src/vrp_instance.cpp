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
	if (has_key(j, "precedence_matrix"))
	{
		instance.prec = Matrix<bool>(n, n);
		for (Vertex i: instance.D.Vertices())
			for (Vertex k: instance.D.Vertices())
				instance.prec[i][k] = (bool) j["precedence_matrix"][i][k];

		instance.prec_count = vector<int>(n, 0);
		instance.suc_count = vector<int>(n, 0);
		for (Vertex i: instance.D.Vertices())
		{
			for (Vertex k: instance.D.Vertices())
			{
				if (instance.prec[i][k])
				{
					instance.prec_count[k]++;
					instance.suc_count[i]++;
				}
			}
		}
	}
}

VRPInstance reverse_instance(const VRPInstance& vrp)
{
	int n = vrp.D.VertexCount();

	VRPInstance rev;
	rev.D = vrp.D.Reverse();
	rev.o = vrp.d, rev.d = vrp.o;
	rev.T = vrp.T;
	for (Vertex v: vrp.D.Vertices()) rev.tw.push_back({-vrp.tw[v].right, -vrp.tw[v].left});
	for (Vertex v: vrp.D.Vertices()) rev.a.push_back(rev.tw[v].left);
	for (Vertex v: vrp.D.Vertices()) rev.b.push_back(rev.tw[v].right);
	rev.prec = Matrix<bool>(n, n, false);
	rev.prec_count = vector<int>(n, 0);
	rev.suc_count = vector<int>(n, 0);
	for (Vertex v: vrp.D.Vertices())
	{
		for (Vertex w: vrp.D.Vertices())
		{
			if (vrp.prec[w][v])
			{
				rev.prec[v][w] = true;
				rev.prec_count[w]++;
				rev.suc_count[v]++;
			}
		}
	}
	rev.LDT = rev.EAT = Matrix<double>(n, n);
	for (Vertex v: vrp.D.Vertices())
	{
		for (Vertex w: vrp.D.Vertices())
		{
			rev.EAT[v][w] = -vrp.LDT[w][v];
			rev.LDT[v][w] = -vrp.EAT[w][v];
		}
	}
	rev.arr = rev.tau = rev.dep = rev.pretau = Matrix<PWLFunction>(n, n);
	for (Vertex u: vrp.D.Vertices())
	{
		for (Vertex v: vrp.D.Successors(u))
		{
	// Compute reverse travel functions.
			rev.tau[v][u] = vrp.pretau[u][v].Compose(PWLFunction::IdentityFunction({-vrp.T, 0.0}) * -1);
			auto init_piece = LinearFunction({-vrp.T, rev.tau[v][u](min(dom(rev.tau[v][u]))) + min(dom(rev.tau[v][u])) + vrp.T}, {min(dom(rev.tau[v][u])), rev.tau[v][u](min(dom(rev.tau[v][u])))});
			rev.tau[v][u] = Min(rev.tau[v][u], PWLFunction({init_piece}));
			rev.arr[v][u] = rev.tau[v][u] + PWLFunction::IdentityFunction(rev.tau[v][u].Domain());
			rev.dep[v][u] = rev.arr[v][u].Inverse();
			rev.pretau[v][u] = PWLFunction::IdentityFunction(rev.dep[v][u].Domain()) - rev.dep[v][u];
		}
	}
	// Add travel functions for (i, i) (for boundary reasons).
	for (Vertex u: rev.D.Vertices())
	{
		rev.tau[u][u] = rev.pretau[u][u] = PWLFunction::ConstantFunction(0.0, rev.tw[u]);
		rev.dep[u][u] = rev.arr[u][u] = PWLFunction::IdentityFunction(rev.tw[u]);
	}
	return rev;
}
} // namespace tdtsptw