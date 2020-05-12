//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include <iostream>
#include <vector>
#include <goc/goc.h>

#include "vrp_instance.h"
#include "preprocess_travel_times.h"
#include "preprocess_time_windows.h"
#include "preprocess_waiting_times.h"
#include "preprocess_remove_time_windows.h"
#include "preprocess_remove_time_dependency.h"
#include "preprocess_time_precedence.h"
#include "ngl_info.h"
#include "pricing_problem.h"
#include "spf.h"
#include "relaxation_solver.h"
#include "exact_solver.h"
#include "label_sequence_ti.h"
#include "label_sequence_td.h"
#include "dynamic_neighbour_augmentation.h"
#include "column_generation.h"

using namespace std;
using namespace goc;
using namespace nlohmann;
using namespace tdtsptw;

namespace
{
Route get_initial_solution(const json& solutions, const string& tag)
{
	Route UB;
	for (auto& solution: solutions)
	{
		set<string> tags = solution["tags"];
		if (includes(tags, tag))
		{
			clog << "> Found: " << solution << endl;
			UB = (Route) solution["routes"][0];
			break;
		}
	}
	return UB;
}

void preprocess_instance(json& instance, const string& objective, bool remove_td, bool remove_time_windows)
{
	// Set departing time from depot equal to 0 if makespan objective.
	if (!has_key(instance, "time_windows") || remove_time_windows) preprocess_remove_time_windows(instance);
	if (objective == "makespan") instance["time_windows"][0] = Interval(0, 0);
	if (remove_td) preprocess_remove_time_dependency(instance);
	if (remove_td) preprocess_remove_time_dependency(instance);
	preprocess_travel_times(instance);

	// Make iterative preprocessing between time windows and waiting times.
	int iter_preprocess = 0;
	preprocess_waiting_times(instance);
	while (preprocess_time_windows(instance))
	{
		clog << "> Preprocessing time windows iteration " << ++iter_preprocess << endl;
		preprocess_waiting_times(instance);
	}
	preprocess_time_precedence(instance);
}
}

int main(int argc, char** argv)
{
	try
	{
		// This line only works in debug mode.
		simulate_runner_input("instances/td-ascheuer", "rbg021.6", "experiments/ascheuer.json", "TI");

		// Read input.
		json experiment, instance, solutions;
		cin >> experiment >> instance >> solutions;
		
		// Parse experiment.
		Duration tl_exact = Duration(value_or_default(experiment, "time_limit_exact", 1200), DurationUnit::Seconds);
		Duration tl_cg = Duration(value_or_default(experiment, "time_limit_cg", 1200), DurationUnit::Seconds);
		Duration tl_dna = Duration(value_or_default(experiment, "time_limit_dna", 1200), DurationUnit::Seconds);
		string objective = value_or_default(experiment, "objective", "duration"); // duration | makespan
		string relaxation = value_or_default(experiment, "relaxation", "NGLTD"); // NGLTD | NGLTI
		bool bounding = value_or_default(experiment, "bounding", true); // true | false
		bool remove_td = value_or_default(experiment, "remove_td", false); // true | false
		bool remove_tw = value_or_default(experiment, "remove_tw", false); // true | false
		string initial_solution_tag = value_or_default(experiment, "initial_solution_tag", "INIT"); // "DFS" | "TILK"
		
		// Show experiment details.
		clog << "Experiment parameters" << endl;
		clog << "> Time limit (sec): (CG) " << tl_cg << " - (DNA) " << tl_dna << " - (EXACT) " << tl_exact << endl;
		clog << "> Objective: " << objective << endl;
		clog << "> Relaxation: " << relaxation << endl;
		clog << "> Bounding: " << bounding << endl;
		clog << "> Remove TD: " << remove_td << endl;
		clog << "> Remove TW: " << remove_tw << endl;
		clog << "> Initial solution tag: " << initial_solution_tag << endl;

		// Preprocess instance JSON.
		clog << "Preprocessing instance..." << endl;
		preprocess_instance(instance, objective, remove_td, remove_tw);
		
		// Parse instance.
		clog << "Parsing instance..." << endl;
		VRPInstance vrp = instance; // forward instance.
		VRPInstance vrp_r = reverse_instance(vrp); // backward instance.

		// Get initial solution.
		clog << "Looking for initial solution with tag " << initial_solution_tag << " in solutions.json..." << endl;
		Route UB = get_initial_solution(solutions, initial_solution_tag);

		// Set output of this process.
		json output;
		BCExecutionLog general_log; // Keep track of general algorithm information.
		general_log.status = BCStatus::Optimum;
		Stopwatch rolex(true), rolex_temp(false); // Keep track of elapsed time.

		// If no initial solution exists, deem the instance as infeasible.
		if (UB.path.empty())
		{
			clog << "No initial solution found. Instance is infeasible." << endl;
			general_log.status = BCStatus::Infeasible;
			general_log.time = rolex.Peek();
			output["General"] = general_log;
			cout << output << endl;
			return 0;
		}
		clog << "> Initial UB: " << UB.duration << endl;

		int n = vrp.D.VertexCount();
		double lb = 0.0;
		vector<double> penalties(n, 0.0); // Initial set of penalties.
		double penalty_sum = 0.0;
		json log;

		// Set relaxation solver for CG and DNA.
		RelaxationSolver cg_relaxation(relaxation == "NGLTI" ? RelaxationSolver::NGLTI : RelaxationSolver::NGLTD, RelaxationSolver::Bidirectional);
		RelaxationSolver dna_relaxation(RelaxationSolver::NGLTD, RelaxationSolver::Bidirectional);

		// Generate NG structures for forward and backward instances.
		NGLInfo ngl_info, ngl_info_r;
		create_default_nginfo(vrp, 3, &ngl_info, &ngl_info_r);

		// Run bounding phase.
		if (bounding)
		{
			// Run initial relaxation with default penalties to get a LB.
			clog << "Running initial relaxation..." << endl;
			Route opt;
			double opt_cost;
			rolex_temp.Resume();
			auto status = cg_relaxation.Run(vrp, vrp_r, ngl_info, ngl_info_r, penalties, nullptr, tl_cg, &opt, &opt_cost, &log);
			rolex_temp.Pause();
			if (status == BLBStatus::TimeLimitReached) clog << "> Time limit reached" << endl;
			lb = max(lb, opt_cost + penalty_sum);
			clog << "> Finished in " << rolex_temp.Peek() << " - LB: " << lb << endl;
			output["initial_relaxation"] = log;
			clog << opt << endl;

			// Solve column generation.
			if (epsilon_different(UB.duration, lb))
			{
				clog << "Running column generation..." << endl;
				rolex_temp.Reset().Resume();
				column_generation(cg_relaxation, vrp, vrp_r, ngl_info, ngl_info_r, tl_cg, &penalties, &UB, &lb, &log);
				rolex_temp.Pause();
				output["column_generation"] = log;
				clog << "> Finished in " << rolex_temp.Peek() << " - LB: " << lb << endl;
				clog << "> Penalties: " << penalties << endl;
			}

			// Solve dynamic neighbour augmentation.
			if (epsilon_different(UB.duration, lb))
			{
				clog << "Running dynamic neighbour augmentation..." << endl;
				rolex_temp.Reset().Resume();
				dynamic_neighbour_augmentation(dna_relaxation, vrp, vrp_r, ngl_info, ngl_info_r, 10, penalties, tl_dna, &UB, &lb, &log);
				rolex_temp.Pause();
				output["dna"] = log;
				clog << "> Finished in " << rolex_temp.Peek() << " - LB: " << lb << endl;
			}
		}

		// Solve exact algorithm.
		if (epsilon_different(UB.duration, lb))
		{
			Stopwatch rolex_exact(true);
			BoundingTree B(&vrp, &ngl_info, penalties);
			if (bounding)
			{
				clog << "Building bounding tree..." << endl;
				rolex_temp.Reset().Resume();
				Route opt;
				double opt_cost;
				dna_relaxation.direction = RelaxationSolver::Backward;
				dna_relaxation.Run(vrp, vrp_r, ngl_info, ngl_info_r, penalties, &B, tl_exact, &opt, &opt_cost, &log);
				rolex_temp.Pause();
				output["bounding_tree_build"] = log;
				clog << "> Finished in " << rolex_temp.Peek() << endl;
			}
			if (!bounding) B.Disable();

			clog << "Running exact algorithm..." << endl;
			rolex_temp.Reset().Resume();
			run_exact(vrp, ngl_info, penalties, B, tl_exact - rolex_exact.Peek(), &lb, &UB, &log);
			rolex_temp.Pause();
			output["exact"] = log;
			clog << "> Finished in " << rolex_temp.Peek() << " - UB: " << UB.duration << endl;
		}

		if (epsilon_equal(lb, INFTY)) clog << "The problem is infeasible." << endl;
		else if (epsilon_equal(lb, UB.duration)) clog << "Optimum found: " << UB << endl;
		else clog << "Gap was not closed. Best LB: " << lb << " - Best UB: " << UB.duration << endl;

		general_log.best_bound = lb;
		general_log.best_int_value = UB.duration;
		general_log.time = rolex.Peek();
		output["general"] = general_log;
		output["best_solution"] = VRPSolution(UB.duration, {UB});
		cout << output << endl;

//		// If there is an initial upper bound, start algorithm
//		if (UB.duration < INFTY)
//		{
//			clog << "Initial UB: " << UB.duration << ", Initial LB: " << LB << endl;
//
//			// Build NG structure that contains all information about neighbourhoods and ngL-tour path.
//			clog << "Building NG structure..." << endl;
//			auto rvrp = reverse_instance(vrp);
//			NGStructure NG(vrp, 3);
//			NGStructure rNG(rvrp, NG.N, reverse(NG.L), NG.delta); // NG for reversed instance.
//
//			vector<double> penalties(vrp.D.VertexCount(), 0.0); // Keep best set of penalties.
//
//			// Run initial iteration of ngL with penalties 0 for initial LB.
//			Stopwatch rolex(true);
//			clog << "Running initial NG with penalties 0" << endl;
//			CGExecutionLog initial_lb_log;
//			initial_lb(vrp, NG, relaxation, UB, LB, &initial_lb_log);
//			output["Initial NG"] = initial_lb_log;
//
//			LPExecutionLog lb_log;
//
//			// Solve CG to obtain best penalties.
//			if (epsilon_smaller(LB, UB.duration))
//			{
//				if (colgen)
//				{
//					clog << "Running CG algorithm..." << endl;
//					// Initialize SPF.
//					SPF spf(vrp.D.VertexCount());
//					spf.AddRoute(UB);
//
//					// Configure CG algorithm.
//					CGSolver cg_solver;
//					LPSolver lp_solver;
//					cg_solver.time_limit = tl_cg;
//					cg_solver.lp_solver = &lp_solver;
//					cg_solver.screen_output = &clog;
//
//					bool early_stop = false; // If any other termination condition is met, early_stop is set to true.
//					cg_solver.pricing_function = [&](const vector<double>& duals, double incumbent_value,
//													 Duration time_limit, CGExecutionLog* cg_execution_log) {
//						if (!colgen) return false;
//						auto pp = spf.InterpretDuals(duals);
//						Route best;
//						double best_cost;
//						if (relaxation == "NGLTI")
//						{
//							MLBExecutionLog iteration_log(true);
//							run_nglti(vrp, NG, pp.penalties, UB.duration, best, best_cost, &iteration_log);
//							cg_execution_log->iterations->push_back(iteration_log);
//						}
//						else if (relaxation == "NGLTD")
//						{
//							MLBExecutionLog iteration_log(true);
//							best = run_ngltd(vrp, NG, pp.penalties, &iteration_log, nullptr, LB, time_limit);
//							best_cost = best.duration - sum<Vertex>(best.path, [&](Vertex v) { return pp.penalties[v]; });
//							cg_execution_log->iterations->push_back(iteration_log);
//						}
//
//						// Compute new LB.
//						if (best_cost + sum(pp.penalties) >= LB)
//						{
//							LB = best_cost + sum(pp.penalties);
//							penalties = pp.penalties;
//						}
//
//						// If LB=UB, we know that UB is optimum, we can stop.
//						if (epsilon_equal(LB, UB.duration))
//						{
//							early_stop = true;
//							clog << "Stop CG (LB = UB)" << endl;
//							return false;
//						}
//
//						// Add the best solution to the RMP if it is not a feasible solution (otherwise it was already added).
//						if (epsilon_smaller(best_cost, 0.0)) spf.AddRoute(best);
//
//						return epsilon_smaller(best_cost, 0.0);
//					};
//
//					// Run CG.
//					auto cg_log = cg_solver.Solve(spf.formulation, {CGOption::IterationsInformation});
//
//					// If gap was not closed, get best LB.
//					if (!early_stop && cg_log.status == CGStatus::Optimum) LB = cg_log.incumbent_value;
//					cg_log.incumbent_value = LB;
//					output["CG"] = cg_log;
//					clog << "Finished CG in " << cg_log.time << "secs with LB: " << LB << endl;
//				}
//
//				bool found_opt = epsilon_equal(LB, UB.duration);
//				if (found_opt)
//				{
//					clog << "Optimality was closed in CG" << endl;
//					lb_log.status = LPStatus::Optimum;
//				}
//
//				clog << "Penalties: " << penalties << endl;
//
//				if (dssr && !found_opt)
//				{
//					clog << "Running DNA to improve bounds." << endl;
//					CGExecutionLog dna_log;
//					auto R = run_dna(vrp, rvrp, NG, rNG, penalties, &dna_log, LB, tl_dna, bidirectional_dna);
//					if (R.duration != INFTY)
//					{
//						UB = R;
//						clog << "\tFound solution " << UB.path << " " << UB.duration << endl;
//						dna_log.status = CGStatus::Optimum;
//						lb_log.status = LPStatus::Optimum;
//						found_opt = true;
//					}
//					clog << "Finished DNA in " << dna_log.time << "s, with LB: " << LB << endl;
//					output["DNA"] = dna_log;
//				}
//
//				// If optimum was not found, run exact algorithm.
//				if (!found_opt)
//				{
//					MLBExecutionLog log_ngl(true);
//					Bounding B(rvrp, rNG, penalties);
//					if (relaxation != "None")
//					{
//						Route R_NG = run_ngltd(rvrp, rNG, penalties, &log_ngl, &B, LB, tl_exact);
//						output["NGL-TD"] = log_ngl;
//						clog << "NGL-TD:   " << LB << "\t" << log_ngl.time << "\t" << log_ngl.processed_count << "\t"
//							 << log_ngl.enumerated_count << endl;
//					}
//
//					MLBExecutionLog log(true);
//					auto r = run_exact_piecewise(vrp, NG.L, penalties, LB, UB.duration, &log,
//												 relaxation == "None" ? nullptr : &B, tl_exact);
//					clog << "Exact: " << r.duration << "\t" << log.time << "\t" << log.processed_count << "\t"
//						 << log.enumerated_count << endl;
//					output["Exact"] = log;
//					if (r.duration < UB.duration) UB = vrp.BestDurationRoute(r.path);
//					lb_log.status = log.status == MLBStatus::TimeLimitReached ? LPStatus::TimeLimitReached : LPStatus::Optimum;
//				}
//			}
//			else
//			{
//				lb_log.status = LPStatus::Optimum;
//			}
//			lb_log.incumbent_value = LB;
//			lb_log.time = rolex.Peek();
//			output["General"] = lb_log;
//
//			// Get best route.
//			if (UB.duration != INFTY)
//			{
//				Route best = UB;
//				clog << "Best solution: " << endl;
//				clog << "\tpath: " << best.path << endl;
//				clog << "\tt0: " << best.t0 << endl;
//				clog << "\tduration: " << best.duration << endl;
//				output["Best solution"] = VRPSolution(best.duration, {best});
//			}
//			output["status"] = *lb_log.status;
//		}
//
//		// Send JSON output to cout.
//		cout << output << endl;
	}
	catch (std::bad_alloc& e)
	{
		return 3;
	}
	return 0;
}