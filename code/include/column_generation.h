//
// Created by Gonzalo Lera Romero on 12/05/2020.
//

#ifndef TDTSPTW_COLUMN_GENERATION_H
#define TDTSPTW_COLUMN_GENERATION_H

#include <vector>
#include "ngl_info.h"
#include "vrp_instance.h"
#include "goc/goc.h"
#include "relaxation_solver.h"

namespace tdtsptw
{
void column_generation(const RelaxationSolver& relaxation, const VRPInstance& vrp_f, const VRPInstance& vrp_b,
					   NGLInfo& ngl_info_f, NGLInfo& ngl_info_b, const goc::Duration& time_limit,
					   std::vector<double>* penalties, goc::Route* UB, double* lb, nlohmann::json* log);
} // namespace tdtsptw

#endif //TDTSPTW_COLUMN_GENERATION_H
