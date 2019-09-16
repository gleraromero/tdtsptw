//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "preprocess_time_windows.h"

#include "vrp_instance.h"
#include <vector>
#include <queue>

using namespace std;
using namespace goc;
using namespace nlohmann;

namespace tdtsptw
{
namespace {
// Returns: the latest we can arrive to k if departing from i (and traversing arc (i, k)) without waiting.
double latest_arrival(const VRPInstance& vrp, Vertex i, Vertex k)
{
	double t = vrp.DepartureTime({i,k}, vrp.b[k]);
	if (t == INFTY) return vrp.ArrivalTime({i,k}, vrp.b[i]);
	return vrp.b[k];
}

// Returns: the earliest we can depart from i, to reach k inside its time window without waiting.
double earliest_departure(const VRPInstance& vrp, Vertex i, Vertex k)
{
	double t = vrp.DepartureTime({i,k}, vrp.a[k]);
	if (t == INFTY) return vrp.a[i];
	return t;
}

// Earliest arrival time from i to all vertices if departing at a_i.
// EAT[i][j] = INFTY if j < i.
vector<double> compute_EAT(const VRPInstance& vrp, Vertex i)
{
	return compute_earliest_arrival_time(vrp.D, i, vrp.a[i], [&] (Vertex u, Vertex v, double t0) {
		return vrp.TravelTime({u,v}, t0);
	});
}

// Latest departure time from all vertices to j if arriving to j at tf.
// LDT[i][j] = -INFTY if j < i.
vector<double> compute_LDT(const VRPInstance& vrp, Vertex j)
{
	auto LDT = compute_latest_departure_time(vrp.D, j, vrp.b[j], [&] (Vertex u, Vertex v, double t0) {
		return vrp.DepartureTime({u,v}, t0);
	});
	for (double& t: LDT) if (t == INFTY) t = -INFTY;
	return LDT;
}

// Removes the arc ij from the instance.
void remove_arc(json& instance, Vertex i, Vertex j)
{
	if (instance["digraph"]["arcs"][i][j] == 0) return;
	instance["digraph"]["arcs"][i][j] = 0;
	int arc_count = instance["digraph"]["arc_count"];
	instance["digraph"]["arc_count"] = arc_count - 1;
	if (has_key(instance, "travel_times")) instance["travel_times"][i][j] = vector<json>({});
}

// Returns: if the instance includes the arc.
bool includes_arc(json& instance, Arc ij)
{
	return instance["digraph"]["arcs"][ij.tail][ij.head] == 1;
}
}

void preprocess_time_windows(json& instance)
{
	for (int iter = 0; iter < 10; ++iter)
	{
		VRPInstance vrp = instance;
		int n = vrp.D.VertexCount();
		auto& V = vrp.D.Vertices();
		auto a = [&](Vertex i) -> double { return instance["time_windows"][i][0]; };
		auto b = [&](Vertex i) -> double { return instance["time_windows"][i][1]; };
		Vertex o = vrp.o;
		Vertex d = vrp.d;
		auto set_a = [&](Vertex i, double t) { instance["time_windows"][i][0] = t; };
		auto set_b = [&](Vertex i, double t) { instance["time_windows"][i][1] = t; };
		
		// Initialize EAT, LDT.
		Matrix<double> EAT(n, n), LDT(n, n);
		for (int i = 0; i < n; ++i) EAT[i] = compute_EAT(vrp, i);
		for (int j = 0; j < n; ++j) LDT[j] = compute_LDT(vrp, j);
		// Transpose LDT so LDT[i][j] is latest departure time from i to reach j.
		for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) swap(LDT[i][j], LDT[j][i]);
		
		// Set EAT and LDT in the instance.
		instance["EAT"] = EAT;
		instance["LDT"] = LDT;
		instance["precedence_matrix"] = Matrix<bool>(n, n, false); // P[i][j] = i < j.
		
		// Initialize BEFORE(k) = { i | EAT(k, i) > b_i }.
		vector <vector<Vertex>> BEFORE(n);
		for (Vertex k: V)
			for (Vertex i: exclude(V, {k}))
				if (vrp.D.IncludesArc({i, k}))
					if (epsilon_bigger(EAT[k][i], b(i)) || i == vrp.o) BEFORE[k].push_back(i);
		
		// Initialize AFTER(k) = { j | EAT(j, k) > b_k }.
		vector <vector<Vertex>> AFTER(n);
		for (Vertex k: V)
			for (Vertex j: exclude(V, {k}))
				if (vrp.D.IncludesArc({k, j}))
					if (epsilon_bigger(EAT[j][k], b(k)) || j == vrp.d) AFTER[k].push_back(j);
		
		// Rule 1: (3.12) 	Upper bound adjustment derived from the latest arrival time at node k from its predecessors,
		//					for k \in N - {o, d}.
		for (Vertex k:exclude(V, {o}))
		{
			double max_arrival = -INFTY;
			for (Vertex i: vrp.D.Predecessors(k)) max_arrival = max(max_arrival, latest_arrival(vrp, i, k));
			set_b(k, min(b(k), max(a(k), max_arrival)));
		}
		
		// Rule 2: (3.13)	Lower bound adjustment derived from the earliest departure time from node k to its successors,
		//					for k \in N - {o,d}.
		for (Vertex k:exclude(V, {d}))
		{
			double min_dep = INFTY;
			for (Vertex j: vrp.D.Successors(k)) min_dep = min(min_dep, earliest_departure(vrp, k, j));
			set_a(k, max(a(k), min(b(k), min_dep)));
		}
		
		// Rule 0: (3.11a) 	Lower bound adjustment derived from the earliest arrival time at node k from its predecessors,
		// 					for k \in N - {o,d}.
		for (Vertex k: exclude(V, {o}))
		{
			double max_EAT = -INFTY;
			for (Vertex i: BEFORE[k]) max_EAT = max(max_EAT, EAT[i][k]);
			set_a(k, max(a(k), max_EAT));
		}
		
		// Rule 3: (3.14a) 	Upper bound adjustment derived from the latest departure time from node k to its successors,
		// 					for k \in N - {o,d}.
		for (Vertex k: exclude(V, {d}))
		{
			double min_LDT = INFTY;
			for (Vertex j: AFTER[k]) min_LDT = min(min_LDT, LDT[k][j]);
			set_b(k, min(b(k), max(a(k), min_LDT)));
		}
		
		// Remove infeasible tw arcs.
		for (Arc ij: vrp.D.Arcs())
		{
			int i = ij.tail, j = ij.head;
			if (vrp.ArrivalTime({i, j}, a(i)) == INFTY)
			{
				instance["precedence_matrix"][j][i] = true;
				remove_arc(instance, i, j);
			}
		}
		for (Vertex i: V)
			for (Vertex j: V)
				if (epsilon_bigger(EAT[i][j], b(j)))
					instance["precedence_matrix"][j][i] = true;
		
		// Rule 4: If exists k such that EAT(k, i) > LDT(i, j) and EAT(i, j) > LDT(j, k), then remove arc (i, j).
		for (Vertex i: V)
		{
			for (Vertex j: V)
			{
				for (Vertex k: V)
				{
					if (k == i || k == j || i == j) continue;
					if (epsilon_bigger(EAT[k][i], LDT[i][j]) && epsilon_bigger(EAT[i][j], LDT[j][k]))
					{
						instance["predecence_matrix"][i][k] = instance["predecence_matrix"][k][j] = true;
						remove_arc(instance, i, j);
						break;
					}
				}
			}
		}
		
		for (Vertex i: exclude(V, {o, d}))
		{
			instance["precedence_matrix"][o][i] = true;
			instance["precedence_matrix"][i][d] = true;
		}
		instance["precedence_matrix"][o][d] = true;
	}
}
} // namespace tdtsptw