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
#include "relaxation_solver.h"
#include "exact_solver.h"
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
		VRPInstance vrp_f = instance; // forward instance.
		VRPInstance vrp_b = reverse_instance(vrp_f); // backward instance.

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

		int n = vrp_f.D.VertexCount();
		double lb = 0.0;
		vector<double> penalties(n, 0.0); // Initial set of penalties.
		double penalty_sum = 0.0;
		json log;

		// Set relaxation solver for CG and DNA.
		RelaxationSolver cg_relaxation(relaxation == "NGLTI" ? RelaxationSolver::NGLTI : RelaxationSolver::NGLTD, RelaxationSolver::Bidirectional);
		RelaxationSolver dna_relaxation(RelaxationSolver::NGLTD, RelaxationSolver::Bidirectional);

		// Generate NG structures for forward and backward instances.
		NGLInfo ngl_info_f, ngl_info_b;
		create_default_nginfo(vrp_f, 3, &ngl_info_f, &ngl_info_b);

		// Run bounding phase.
		if (bounding)
		{
			// Run initial relaxation with default penalties to get a LB.
			clog << "Running initial relaxation..." << endl;
			Route opt;
			double opt_cost;
			rolex_temp.Resume();
			auto status = cg_relaxation.Run(vrp_f, vrp_b, ngl_info_f, ngl_info_b, penalties, nullptr, tl_cg, &opt, &opt_cost, &log);
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
				column_generation(cg_relaxation, vrp_f, vrp_b, ngl_info_f, ngl_info_b, tl_cg, &penalties, &UB, &lb, &log);
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
				dynamic_neighbour_augmentation(dna_relaxation, vrp_f, vrp_b, ngl_info_f, ngl_info_b, 15, penalties, tl_dna, &UB, &lb, &log);
				rolex_temp.Pause();
				output["dna"] = log;
				clog << "> Finished in " << rolex_temp.Peek() << " - LB: " << lb << endl;
			}
		}

		// Solve exact algorithm.
		if (epsilon_different(UB.duration, lb))
		{
			Stopwatch rolex_exact(true);
			BoundingTree B(&vrp_f, &ngl_info_f, penalties);
			if (bounding)
			{
				clog << "Building bounding tree..." << endl;
				rolex_temp.Reset().Resume();
				Route opt;
				double opt_cost;
				dna_relaxation.direction = RelaxationSolver::Backward;
				dna_relaxation.Run(vrp_f, vrp_b, ngl_info_f, ngl_info_b, penalties, &B, tl_exact, &opt, &opt_cost, &log);
				rolex_temp.Pause();
				output["bounding_tree_build"] = log;
				clog << "> Finished in " << rolex_temp.Peek() << endl;
			}
			if (!bounding) B.Disable();

			clog << "Running exact algorithm..." << endl;
			rolex_temp.Reset().Resume();
			run_exact(vrp_f, ngl_info_f, penalties, B, tl_exact - rolex_exact.Peek(), &lb, &UB, &log);
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
	}
	catch (std::bad_alloc& e)
	{
		return 3; // Memory limit exceeded code for runner.
	}
	return 0;
}