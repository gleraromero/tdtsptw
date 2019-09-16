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

namespace tdtsptw
{
goc::Route initial_heuristic(const VRPInstance& vrp, std::vector<goc::Vertex>& P, VertexSet S, double t);
} // namespace tdtsptw

#endif //TDTSPTW_LABELING_H
