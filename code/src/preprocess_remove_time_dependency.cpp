//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "preprocess_remove_time_dependency.h"

using namespace std;
using namespace nlohmann;

namespace tdtsptw
{
void preprocess_remove_time_dependency(json& instance)
{
	for (auto& cluster_row: instance["cluster_speeds"])
		for (auto& speed: cluster_row)
			speed = 1.0;
}
} // namespace tdtsptw