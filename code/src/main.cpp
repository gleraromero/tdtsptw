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
		simulate_runner_input("instances/adamo_et_al_2019", "15_80_A_A5", "experiments/adamo.json", "TI");

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
		RelaxationSolver cg_relaxation(relaxation == "NGLTI" ? RelaxationSolver::NGLTI : RelaxationSolver::NGLTD, RelaxationSolver::Forward);
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
//				column_generation(cg_relaxation, vrp_f, vrp_b, ngl_info_f, ngl_info_b, tl_cg, &penalties, &UB, &lb, &log);
				penalties = {0, 17052.9, 0, 0, 0, 0, 0, 7.30061e-09, 7.9142e-08, 0, 0, 6.65558e-08, 2.80417e-08, 2.65634e-08, 2.80412e-08, 9.54469e-10, 0, -1.8725e-08, 4.64024e-08, 6.1664e-08, 2.94141e-08, 4.63357e-08, 6.65135e-08, 0, 5.47962e-08, 4.2529e-08, 3.5181e-08, 7.60384e-08, 0, 5.34757e-08, 8.98529e-08, 2.68107e-08, 8.75151e-08, 4.3357e-08, 3.38174e-08, 7.75771e-08, 1.00652e-07, 0, 8.66318e-08, 7.70041e-08, 4.6384e-08, 5.29522e-08, 8.98601e-08, 9.98032e-08, 4.78732e-08, 2.81475e-08, 8.33951e-08, 3.80988e-08, 3.26818e-08, 3.54091e-08, 6.46085e-08, 5.25569e-08, 0, -2.67596e-09, -1.61579e-08, 2.01494e-09, 6.74722e-08, 1.32652e-09, 4.47528e-08, 0, 4.34219e-08, 4.6421e-08, 6.81845e-08, 5.67269e-08, 2.6306e-08, 4.72839e-08, 1.71332e-08, 3.65428e-08, 2.28202e-08, 8.94036e-08, 8.76387e-08, 1.34991e-08, 0, 6.99816e-08, 5.50127e-08, 2.91022e-08, -3.77274e-09, 9.89243e-09, -1.14773e-08, -2.02828e-08, 2.13529e-08, 0, 5.04662e-08, -2.86895e-08, -7.23825e-08, 1.57487e-09, 0, 0, 0, -3.1207e-08, 9.64077e-10, -2.99348e-08, 4.76065e-08, -1.38276e-08, -5.68819e-08, 0, -1.09439e-07, -3.39099e-08, -9.77524e-08, 2.88166e-08, -1.06677e-08, 0, -1.7618e-08, -1.99955e-09, 7.90731e-08, 8.15821e-09, -4.24225e-08, -7.18517e-09, -1.13148e-07, 3.14721e-08, -3.03218e-08, -3.16865e-08, 3.21637e-08, -9.1292e-10, 0, -3.01908e-08, -1.10079e-08, 0, 7.97336e-08, 5.52414e-08, 3.69299e-08, 1.18646e-09, 1.19221e-07, 1.07491e-07, 3.47092e-08, 0, 6.94508e-08, 8.35403e-08, 7.63492e-08, -1.72996e-10, 8.41931e-08, 1.00925e-07, 4.85605e-08, 2.67505e-08, 8.18936e-08, 4.71766e-08, 5.14041e-08, 5.99016e-08, 0, 9.94796e-08, 8.45221e-08, 4.82755e-08, 8.22423e-08, 7.65865e-08, 1.285e-07, 1.43262e-07, 1.35685e-07, 8.15154e-08, 8.45311e-08, 7.60527e-08, 1.6158e-07, 6.06689e-08, 1.29226e-07, 0, 7.60454e-08, 4.31036e-07, 2.55556, 6.72222, 0.777778, 0.166668, 11.0556, 2.55556, 0, 11.0556, 13.8333, 20.7222, 17.3889, 19.9444, 12.8333, 5.94445, 0, 6.16667, 28.8333, 0};
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