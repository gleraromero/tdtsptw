//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef TDTSPTW_LBL_EXACT_H
#define TDTSPTW_LBL_EXACT_H

#include <vector>
#include <unordered_map>
#include <tuple>

#include "goc/goc.h"

#include "vrp_instance.h"
#include "pricing_problem.h"
#include "dssr.h"

namespace tdtsptw
{
// Runs an NG labeling algorithm to find negative cost routes.
// @param vrp: VRP Instance
// @param NG: structure with NG information.
// @param B: bounding structure for completion bounds.
// @param lambda: penalties for vertices.
// @param UB: upper bound on the duration of the optimal route.
// @param [out] best_route: route with best cost.
// @param [out] log: output log to save the execution information.
// @returns the optimal tour.
goc::Route run_exact(const VRPInstance& vrp, const NGStructure& NG, BoundingStructure& B,
					 const std::vector<double>& lambda, const goc::Route& UB, double lb, goc::MLBExecutionLog* log);
	

class MinFunc
{
public:
	std::vector<goc::LinearFunction> Merge(const std::vector<goc::LinearFunction>& F);
	
	const goc::LinearFunction& operator[](int index);
	
private:
	goc::PWLFunction pieces_;
};

goc::Route run_exact_piecewise(const VRPInstance& vrp, const goc::GraphPath& L, const std::vector<double>& lambda,
						  double LB, double UB, goc::MLBExecutionLog* log);
} // namespace tdtsptw

#endif //TDTSPTW_LBL_EXACT_H
