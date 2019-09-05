//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef TDTSPTW_LABELING_H
#define TDTSPTW_LABELING_H

#include <vector>

#include "goc/goc.h"
#include "vrp_instance.h"
#include "label.h"

namespace tdtsptw
{
goc::Route initial_heuristic(const VRPInstance& vrp, std::vector<goc::Vertex>& P, VertexSet S, double t);

goc::Route run_labeling(const VRPInstance& vrp, const goc::Duration& time_limit, goc::MLBExecutionLog* log, bool t0_is_zero);

std::vector<goc::Route> run_ng_labeling(const VRPInstance& vrp, const goc::Duration& time_limit, goc::MLBExecutionLog* log, bool t0_is_zero, const std::vector<double>& penalties);
} // namespace tdtsptw

#endif //TDTSPTW_LABELING_H
