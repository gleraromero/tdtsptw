//
// Created by Gonzalo Lera Romero on 27/04/2020.
//

#ifndef TDTSPTW_RELAXATION_SOLVER_H
#define TDTSPTW_RELAXATION_SOLVER_H

#include <vector>

#include "goc/goc.h"

#include "vrp_instance.h"
#include "ngl_info.h"
#include "bounding_tree.h"

namespace tdtsptw
{
class RelaxationSolver
{
public:
	enum Type { NGLTI, NGLTD };
	enum Direction { Forward, Backward, Bidirectional };
	Type type;
	Direction direction;
	bool asymmetric;

	RelaxationSolver() = default;

	RelaxationSolver(Type type, Direction direction, bool asymmetric);

	goc::BLBStatus Run(const VRPInstance& vrp_f, const VRPInstance& vrp_b, const NGLInfo& ngl_info_f,
					   const NGLInfo& ngl_info_b, const std::vector<double>& penalties, BoundingTree* B,
					   const goc::Duration& time_limit, goc::Route* opt, double* opt_cost, nlohmann::json* log) const;
};
} // namespace tdtsptw

#endif //TDTSPTW_RELAXATION_SOLVER_H
