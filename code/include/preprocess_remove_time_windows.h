//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef TDTSPTW_PREPROCESS_REMOVE_TIME_WINDOWS_H
#define TDTSPTW_PREPROCESS_REMOVE_TIME_WINDOWS_H

#include <goc/goc.h>

namespace tdtsptw
{
// Takes a JSON instance of the vehicle routing problems with the following attributes:
//	- horizon
//	- time_windows
// Sets the time windows equal to the horizon.
void preprocess_remove_time_windows(nlohmann::json& instance);
} // namespace tdtsptw

#endif //TDTSPTW_PREPROCESS_REMOVE_TIME_WINDOWS_H
