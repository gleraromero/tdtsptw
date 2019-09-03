//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef TDTSPTW_PREPROCESS_WAITING_TIMES_H
#define TDTSPTW_PREPROCESS_WAITING_TIMES_H

#include <goc/goc.h>

namespace tdtsptw
{
// Takes a JSON instance of the vehicle routing problems with the following attributes:
//	- digraph
//	- travel_times
//	- time_windows (optional)
//	- service_times (optional)
// We transform the instances into ones without waiting times. In order to do this we apply
// the preprocessing technique from Lera-Romero & Miranda-Bront 2019 (2.1).
// 	(i) 	\tau'_ij(t) = max(a_j, t+\tau_ij(t)) - t for each ij \in A.
	void preprocess_waiting_times(nlohmann::json& instance);
} // namespace tdtsptw

#endif //TDTSPTW_PREPROCESS_WAITING_TIMES_H
