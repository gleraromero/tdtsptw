//
// Created by Gonzalo Lera Romero on 12/05/2020.
//

#ifndef TDTSPTW_PREPROCESS_TIME_PRECEDENCE_H
#define TDTSPTW_PREPROCESS_TIME_PRECEDENCE_H

#include <goc/goc.h>

namespace tdtsptw
{
// Takes a JSON instance of the vehicle routing problems with the following attributes:
//	- digraph
//	- start_depot
//	- end_depot
//	- time_windows
//	- LDT
//	- EAT
// Sets the time windows of each vertex given that they are reached at the k-th position of a route.
// TPW[k][v] should be the interval where a route of length k can arrive at vertex v.
void preprocess_time_precedence(nlohmann::json& instance);
} // namespace tdtsptw

#endif //TDTSPTW_PREPROCESS_TIME_PRECEDENCE_H
