//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "preprocess_remove_time_windows.h"
#include <vector>
#include "goc/goc.h"

using namespace std;
using namespace goc;
using namespace nlohmann;

namespace tdtsptw
{
void preprocess_remove_time_windows(json& instance)
{
	int n = instance["digraph"]["vertex_count"];
	instance["time_windows"] = vector<Interval>(n);
	for (int i = 0; i < n; ++i) instance["time_windows"][i] = instance["horizon"];
}
} // namespace tdtsptw