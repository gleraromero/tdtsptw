//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef TDTSPTW_INITIAL_LB_H
#define TDTSPTW_INITIAL_LB_H

#include <vector>
#include <tuple>

#include "goc/goc.h"

#include "vrp_instance.h"
#include "labeling.h"
#include "spf.h"

namespace tdtsptw
{
// Runs an NGL-tour relaxation algorithm to obtain a LB.
// @param vrp: VRP Instance
// @param NG: structure with NG information.
// @param [in/out] LB: lower bound on the duration of the optimal route.
// @param [in/out] UB: upper bound on the duration of the optimal route.
// @param [out] log: output log to save the execution information.
// @returns the optimal route of the NGL-tour relaxation.
std::vector<goc::Route> initial_lb(const VRPInstance& vrp, const NGStructure& NG, const std::string& relaxation, goc::Route& UB, double& LB, goc::CGExecutionLog* log);
} // namespace tdtsptw

#endif //TDTSPTW_INITIAL_LB_H
