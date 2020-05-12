//
// Created by Gonzalo Lera Romero on 12/05/2020.
//

#ifndef TDTSPTW_COLUMN_GENERATION_H
#define TDTSPTW_COLUMN_GENERATION_H

#include <vector>
#include "core.h"
#include "ngl_info.h"
#include "vrp_instance.h"
#include "goc/goc.h"
#include "relaxation_solver.h"
#include "label_sequence_ti.h"
#include "label_sequence_td.h"
#include "spf.h"
#include "pricing_problem.h"

namespace tdtsptw
{
template<typename LS>
void column_generation(const VRPInstance& vrp_f, const VRPInstance& vrp_b, NGLInfo& ngl_info_f,
					   NGLInfo& ngl_info_b, const goc::Duration& time_limit,
					   std::vector<double>* penalties, goc::Route* UB, double* lb, nlohmann::json* log)
{
	SPF spf(vrp_f.D.VertexCount());
	spf.AddRoute(*UB);
	goc::CGSolver cg_solver;
	goc::LPSolver lp_solver;
	cg_solver.time_limit = time_limit;
	cg_solver.lp_solver = &lp_solver;
	cg_solver.screen_output = &std::clog;

	// Use as pricing function the relaxation.
	cg_solver.pricing_function = [&](const std::vector<double>& duals, double incumbent_value,
									 goc::Duration time_limit, goc::CGExecutionLog* cg_execution_log) {
		goc::Route opt;
		double opt_cost;
		nlohmann::json log;
		auto pricing_problem = spf.InterpretDuals(duals);
		auto status = run_relaxation<LS>(vrp_f, vrp_b, ngl_info_f, ngl_info_b, pricing_problem.penalties, nullptr, time_limit, &opt, &opt_cost, &log);
		cg_execution_log->iterations->push_back(log);

		if (status == goc::BLBStatus::TimeLimitReached) { std::clog << "> Time limit reached" << std::endl; return false; }

		// Check if a new lower bound is found by the solution of the relaxation.
		double pp_penalties_sum = goc::sum(pricing_problem.penalties);
		if (opt_cost + pp_penalties_sum > *lb)
		{
			*lb = opt_cost + pp_penalties_sum;
			*penalties = pricing_problem.penalties;
			std::clog << "> Found new lower bound: " << *lb << std::endl;
		}

		// Check if the solution of the relaxation is elementary.
		if (goc::epsilon_equal(*lb, UB->duration))
		{
			std::clog << "> Optimum was found since LB reached UB" << std::endl;
			return false;
		}

		// Check if a negative reduced cost route was found, if so, add it to the SPF.
		if (goc::epsilon_smaller(opt_cost, 0.0)) spf.AddRoute(opt);

		// Keep iterating while a route was added to the SPF.
		return goc::epsilon_smaller(opt_cost, 0.0);
	};

	auto result = cg_solver.Solve(spf.formulation, {goc::CGOption::IterationsInformation});
	*log = result; // Save algorithm result in output variable.
	if (goc::epsilon_equal(UB->duration, *lb)) return; // If optimum was found, do nothing else.
	if (result.status == goc::CGStatus::Optimum) *lb = result.incumbent_value; // Update lb to the incumbent of CG if finished.
	if (result.status == goc::CGStatus::TimeLimitReached) std::clog << "> Time limit reached" << std::endl;
}
} // namespace tdtsptw

#endif //TDTSPTW_COLUMN_GENERATION_H
