//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef TDTSPTW_PREPROCESS_TIME_WINDOWS_H
#define TDTSPTW_PREPROCESS_TIME_WINDOWS_H

#include <goc/goc.h>

namespace tdtsptw
{
// Takes a JSON instance of the vehicle routing problems with the following attributes:
//	- digraph
//	- travel_times
//	- time_windows
//	- start_depot
//	- end_depot
// Assumes preprocess_travel_times was called (i.e. instance has no service nor waiting times).
// Shrinks time windows [a_i, b_i] according to the techinques introduced in:
//	Desrosiers, J., Dumas, Y., Solomon, M. M., & Soumis, F. (1995).
// and removes infeasible arcs.
bool preprocess_time_windows(nlohmann::json& instance);
} // namespace tdtsptw

#endif //TDTSPTW_PREPROCESS_TIME_WINDOWS_H
