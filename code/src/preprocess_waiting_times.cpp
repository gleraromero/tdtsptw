//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "preprocess_waiting_times.h"

using namespace std;
using namespace goc;
using namespace nlohmann;

namespace tdtsptw
{
void preprocess_waiting_times(json& instance)
{
	Digraph D = instance["digraph"];
	Matrix<PWLFunction> tau = instance["travel_times"];
	auto a = [&] (Vertex i) -> double { return instance["time_windows"][i][0]; };
	auto b = [&] (Vertex i) -> double { return instance["time_windows"][i][1]; };

	/// (i) \tau'_ij(t) = max(a_j, t+\tau_ij(t)) - t for each ij \in A.
	for (Vertex i: D.Vertices())
	{
		for (Vertex j: D.Successors(i))
		{
			if (!tau[i][j].Domain().Includes(a(i))) continue;
			PWLFunction id = PWLFunction::IdentityFunction(tau[i][j].Domain());
			PWLFunction arr = Max(a(j), id+tau[i][j]);
			arr = Max(arr, arr.Value(a(i)));
			arr = arr.RestrictImage({a(j), b(j)});
			arr = arr.RestrictDomain({a(i), b(i)});
			tau[i][j] = arr - id;
		}
	}
	instance["travel_times"] = tau;
}
} // namespace tdtsptw