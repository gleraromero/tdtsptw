//
// Created by Gonzalo Lera Romero on 12/05/2020.
//

#include "column_generation.h"

#include "spf.h"
#include "pricing_problem.h"

using namespace std;
using namespace goc;
using namespace nlohmann;

namespace tdtsptw
{
void column_generation(const RelaxationSolver& relaxation, const VRPInstance& vrp_f, const VRPInstance& vrp_b,
					   NGLInfo& ngl_info_f, NGLInfo& ngl_info_b, const Duration& time_limit,
					   vector<double>* penalties, Route* UB, double* lb, nlohmann::json* log)
{
	SPF spf(vrp_f.D.VertexCount());
	spf.AddRoute(*UB);
	CGSolver cg_solver;
	LPSolver lp_solver;
	cg_solver.time_limit = time_limit;
	cg_solver.lp_solver = &lp_solver;
	cg_solver.screen_output = &clog;

	// Use as pricing function the relaxation.
	cg_solver.pricing_function = [&](const vector<double>& duals, double incumbent_value,
									 Duration time_limit, CGExecutionLog* cg_execution_log) {
		Route opt;
		double opt_cost;
		nlohmann::json log;
		auto pricing_problem = spf.InterpretDuals(duals);
		auto status = relaxation.Run(vrp_f, vrp_b, ngl_info_f, ngl_info_b, pricing_problem.penalties, nullptr,
									 time_limit, &opt, &opt_cost, &log);
		cg_execution_log->iterations->push_back(log);

		if (status == BLBStatus::TimeLimitReached) { clog << "> Time limit reached" << endl; cg_execution_log->status = CGStatus::TimeLimitReached; return false; }

		// Check if a new lower bound is found by the solution of the relaxation.
		double pp_penalties_sum = sum(pricing_problem.penalties);
		if (opt_cost + pp_penalties_sum > *lb)
		{
			*lb = opt_cost + pp_penalties_sum;
			*penalties = pricing_problem.penalties;
			clog << "> Found new lower bound: " << *lb << endl;
		}

		// Check if the solution of the relaxation is elementary.
		if (epsilon_equal(*lb, UB->duration))
		{
			clog << "> Optimum was found since LB reached UB" << endl;
			return false;
		}

		// Check if a negative reduced cost route was found, if so, add it to the SPF.
		if (epsilon_smaller(opt_cost, 0.0)) spf.AddRoute(opt);

		// Keep iterating while a route was added to the SPF.
		return epsilon_smaller(opt_cost, 0.0);
	};

	auto result = cg_solver.Solve(spf.formulation, {CGOption::IterationsInformation});
	*log = result; // Save algorithm result in output variable.
	if (epsilon_equal(UB->duration, *lb)) return; // If optimum was found, do nothing else.
	if (result.status == CGStatus::Optimum) *lb = result.incumbent_value; // Update lb to the incumbent of CG if finished.
	if (result.status == CGStatus::TimeLimitReached) clog << "> Time limit reached" << endl;
}
} // namespace tdtsptw