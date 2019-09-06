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
#include "labeling.h"
#include "pricing_problem.h"
#include "spf.h"

using namespace std;
using namespace goc;
using namespace nlohmann;
using namespace tdtsptw;

int main(int argc, char** argv)
{
	json output; // STDOUT output will go into this JSON.
	
	simulate_runner_input("instances/afg_1995", "rbg021.6", "experiments/afg.json", "duration");
	
	json experiment, instance, solutions;
	cin >> experiment >> instance >> solutions;
	
	// Parse experiment.
	Duration time_limit = Duration(value_or_default(experiment, "time_limit", 7200), DurationUnit::Seconds);
	string objective = value_or_default(experiment, "objective", "makespan");
	
	// Show experiment details.
	clog << "Time limit: " << time_limit << "sec." << endl;
	clog << "Objective: " << objective << endl;
	
	// Preprocess instance JSON.
	preprocess_travel_times(instance);
	preprocess_waiting_times(instance);
	preprocess_time_windows(instance);
	
	// Parse instance.
	VRPInstance vrp = instance;
	SPF spf(vrp.D);
	vector<Vertex> P = {vrp.o};
	Route UB = initial_heuristic(vrp, P, create_bitset<MAX_N>({vrp.o}), vrp.tw[vrp.o].left);
	clog << "UB: " << UB << endl;
	spf.AddRoute(UB);
	CGSolver cg_solver;
	LPSolver lp_solver;
	cg_solver.time_limit = 100.0_sec;
	cg_solver.lp_solver = &lp_solver;
	cg_solver.screen_output = &clog;
	cg_solver.pricing_function = [&] (const vector<double>& duals, double incumbent_value, Duration time_limit, CGExecutionLog* cg_execution_log)
	{
		MLBExecutionLog iteration_log(true);
		auto R = run_ng_labeling(vrp, time_limit, &iteration_log, objective == "makespan", spf.InterpretDuals(duals).penalties);
		for (auto& r: R) spf.AddRoute(r);
		cg_execution_log->iterations->push_back(iteration_log);
		return !R.empty();
	};
	CGExecutionLog cg_log;
	cg_log = cg_solver.Solve(spf.formulation, {CGOption::IterationsInformation});
	clog << "Finished CG in " << cg_log.time << "secs." << endl;
	
	// Solve Exact.
	MLBExecutionLog exact_log;
	Route best = run_labeling(vrp, time_limit, &exact_log, objective == "makespan");
	
	// Show results.
	clog << "Time: " << exact_log.time << endl;
	clog << "Status: " << exact_log.status << endl;
	if (!best.path.empty())
	{
		clog << "Best route: " << best.path << ", departing at: " << best.t0 << ", duration: " << best.duration << endl;
		output["Best solution"] = VRPSolution(best.duration, {best});
	}
	output["Exact"] = exact_log;
	output["NG"] = cg_log;

	// Send JSON output to cout.
	cout << output << endl;
	return 0;
}