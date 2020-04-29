//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef TDTSPTW_PREPROCESS_REMOVE_TIME_DEPENDENCY_H
#define TDTSPTW_PREPROCESS_REMOVE_TIME_DEPENDENCY_H

#include <goc/goc.h>

namespace tdtsptw
{
// Takes a JSON instance of the vehicle routing problems with the following attributes:
//	- cluster_speeds
// Sets all speeds to 1.0.
void preprocess_remove_time_dependency(nlohmann::json& instance);
} // namespace tdtsptw

#endif //TDTSPTW_PREPROCESS_REMOVE_TIME_DEPENDENCY_H
