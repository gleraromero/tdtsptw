//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef TDTSPTW_PREPROCESS_TRAVEL_TIMES_H
#define TDTSPTW_PREPROCESS_TRAVEL_TIMES_H

#include <goc/goc.h>

namespace tdtsptw
{
// Takes a JSON instance of the vehicle routing problems with the following attributes:
//	- digraph
//	- distances
//	- clusters
//	- cluster_speeds
//	- speed_zones
//	- service_times (optional)
//	- time_windows (optional).
// Adds the travel_times attribute to a JSON instance with a matrix of piecewise linear functions.
void preprocess_travel_times(nlohmann::json& instance);
} // namespace tdtsptw

#endif //TDTSPTW_PREPROCESS_TRAVEL_TIMES_H
