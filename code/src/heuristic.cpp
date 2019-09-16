//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "heuristic.h"
#include "pwl_domination_function.h"

using namespace std;
using namespace goc;

namespace tdtsptw
{
Route initial_heuristic(const VRPInstance& vrp, vector<Vertex>& P, VertexSet S, double t)
{
	int n = vrp.D.VertexCount();
	Vertex u = P.back();
	if (P.size() == n)
	{
		if (u == vrp.d) return vrp.BestDurationRoute(P);
		else return Route({}, 0.0, INFTY);
	}
	for (Vertex w: vrp.D.Vertices()) if (!S.test(w) && epsilon_smaller(vrp.LDT[u][w], t)) return Route({}, 0.0, INFTY);
	
	for (Vertex v: vrp.D.Successors(P.back()))
	{
		if (S.test(v)) continue;
		P.push_back(v);
		Route r = initial_heuristic(vrp, P, unite(S, {v}), vrp.ArrivalTime({u,v}, t));
		if (r.duration != INFTY) return r;
		P.pop_back();
	}
	return Route({}, 0.0, INFTY);
}
} // namespace tdtsptw