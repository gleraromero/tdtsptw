//
// Created by Gonzalo Lera Romero on 07/05/2020.
//

#ifndef TDTSPTW_EXACT_SOLVER_H
#define TDTSPTW_EXACT_SOLVER_H

#include "goc/goc.h"
#include "vrp_instance.h"
#include "ngl_info.h"
#include "bounding_tree.h"
#include "label_sequence_td.h"

namespace tdtsptw
{
goc::MLBStatus run_exact(const VRPInstance& vrp, const NGLInfo& ngl_info, const std::vector<double>& penalties,
						 const BoundingTree<LabelSequenceTD>& B, const goc::Duration &time_limit, double lb,
						 goc::Route* UB, nlohmann::json *log);
}

#endif //TDTSPTW_EXACT_SOLVER_H
