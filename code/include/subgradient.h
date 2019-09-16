//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef TDTSPTW_SUBGRADIENT_H
#define TDTSPTW_SUBGRADIENT_H

#include <vector>
#include <tuple>

#include "goc/goc.h"

#include "vrp_instance.h"
#include "lbl_ng.h"

namespace tdtsptw
{
// Runs a subradient descent algorithm.
// @param vrp: VRP Instance
// @param NG: structure with NG information.
// @param max_iter: number of iterations.
// @param [in/out] LB: lower bound on the duration of the optimal route.
// @param [in/out] UB: upper bound on the duration of the optimal route.
// @param [out] log: output log to save the execution information.
// @returns the set of optimal routes for each iteration.
std::vector<goc::Route> subgradient(const VRPInstance& vrp, const NGStructure& NG, bool use_td_relaxation, int max_iter, goc::Route& UB, double& LB, goc::CGExecutionLog* log);
} // namespace tdtsptw

#endif //TDTSPTW_SUBGRADIENT_H
