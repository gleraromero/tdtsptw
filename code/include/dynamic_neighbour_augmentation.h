//
// Created by Gonzalo Lera Romero on 07/05/2020.
//

#ifndef TDTSPTW_DYNAMIC_NEIGHBOUR_AUGMENTATION_H
#define TDTSPTW_DYNAMIC_NEIGHBOUR_AUGMENTATION_H

#include <vector>
#include "core.h"
#include "ngl_info.h"
#include "vrp_instance.h"
#include "goc/goc.h"
#include "relaxation_solver.h"

namespace tdtsptw
{
// Augments the neighbours from _ngl_info_f_ and _ngl_info_b_ up until size _delta_ by breaking cycles.
// If an elementary route is found then the corresponding UB is set in the output parameter.
// The lb is updated at each iteration since the result gives a lower bound.
// Observation: v \in N[v] for all v \in V, therefore, delta includes such vertex as well.
void dynamic_neighbour_augmentation(const RelaxationSolver& relaxation, const VRPInstance& vrp_f,
									const VRPInstance& vrp_b, NGLInfo& ngl_info_f, NGLInfo& ngl_info_b, int delta,
									const std::vector<double>& penalties, const goc::Duration& time_limit,
									goc::Route* UB, double* lb, nlohmann::json* log);
} // namespace tdtsptw

#endif //TDTSPTW_DYNAMIC_NEIGHBOUR_AUGMENTATION_H
