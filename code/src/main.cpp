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
#include "lbl_ng.h"
#include "lbl_exact.h"
#include "pricing_problem.h"
#include "subgradient.h"
#include "spf.h"
#include "heuristic.h"
#include "dssr.h"

using namespace std;
using namespace goc;
using namespace nlohmann;
using namespace tdtsptw;

int main(int argc, char** argv)
{
	json output; // STDOUT output will go into this JSON.
	
	simulate_runner_input("instances/lms_2019", "rbg132.2", "experiments/lms.json", "NGL");
	
	json experiment, instance, solutions;
	cin >> experiment >> instance >> solutions;
	
	// Parse experiment.
	Duration time_limit = Duration(value_or_default(experiment, "time_limit", 7200), DurationUnit::Seconds);
	string objective = value_or_default(experiment, "objective", "duration");
	string relaxation = value_or_default(experiment, "relaxation", "NGL-TD");
	bool colgen = value_or_default(experiment, "colgen", true);
	bool dssr = value_or_default(experiment, "dssr", true);
	
	// Show experiment details.
	clog << "Time limit: " << time_limit << "sec." << endl;
	clog << "Objective: " << objective << endl;
	clog << "Relaxation: " << relaxation << endl;
	clog << "Colgen: " << colgen << endl;
	clog << "DSSR: " << colgen << endl;
	
	// Set departing time from depot equal to 0 if makespan objective.
	if (objective == "makespan") instance["time_windows"][0] = Interval(0, 0);
	
	// Preprocess instance JSON.
	preprocess_travel_times(instance);
	preprocess_waiting_times(instance);
	preprocess_time_windows(instance);
	preprocess_waiting_times(instance);
	
	// Parse instance.
	clog << "Parsing instance..." << endl;
	VRPInstance vrp = instance;
	
	// Get UB.
	double LB = 0.0;
	vector<Vertex> P = {vrp.o};
	Route UB = initial_heuristic(vrp, P, create_bitset<MAX_N>({vrp.o}), vrp.tw[vrp.o].left);
	
//	Route UB = vrp.BestDurationRoute({0,5,16,17,19,10,11,8,2,3,4,12,1,15,13,6,14,9,18,7,20});
	if (UB.duration == INFTY)
	{
		output["status"] = "Infeasible";
		clog << "Infeasible" << endl;
	}
	
	if (UB.duration < INFTY)
	{
		clog << "Initial UB: " << UB.duration << ", Initial LB: " << LB << endl;
		
		// Build NG structure.
		clog << "Building NG structure..." << endl;
		NGStructure NG(vrp, 3);
		
		{
			vector<double> penalties(vrp.D.VertexCount(), 0.0);
			MLBExecutionLog log(true);
			auto r = run_exact_piecewise(vrp, NG.L, penalties, 0.0, UB.duration, &log);
			json output;
			output["Exact"] = log;
			cout << output << endl;
			clog << r << endl;
			exit(0);
		}
		
		// Run subgradient.
		vector<Route> sg_routes;
//		CGExecutionLog subgradient_log;
//		sg_routes = subgradient(vrp, NG, relaxation == "NGL-TD", 10, UB, LB, &subgradient_log);
//		output["Subgradient"] = subgradient_log;
//
		// Solve CG to obtain best penalties.
		if (epsilon_smaller(LB, UB.duration))
		{
			vector<double> penalties(vrp.D.VertexCount(), 0.0); // Keep best set of penalties.
			
			if (colgen)
			{
				clog << "Running CG algorithm..." << endl;
				// Initialize SPF.
				SPF spf(vrp.D.VertexCount());
				spf.AddRoute(UB);
				for (auto& r: sg_routes) spf.AddRoute(r);
				
				// Configure CG algorithm.
				CGSolver cg_solver;
				LPSolver lp_solver;
				cg_solver.time_limit = time_limit;
				cg_solver.lp_solver = &lp_solver;
				cg_solver.screen_output = &clog;
				
				bool early_stop = false; // If any other termination condition is met, early_stop is set to true.
				cg_solver.pricing_function = [&](const vector<double>& duals, double incumbent_value,
												 Duration time_limit, CGExecutionLog* cg_execution_log) {
					if (!colgen) return false;
					auto pp = spf.InterpretDuals(duals);
					MLBExecutionLog iteration_log(true);
					Route best;
					double best_cost;
					vector<Route> R;
					if (relaxation == "NGL-TD")
						R = run_ng(vrp, NG, pp.penalties, UB.duration, &best, &best_cost, &iteration_log);
					else if (relaxation == "NGL")
						R = run_ng_td(vrp, NG, pp.penalties, UB.duration, &best, &best_cost, &iteration_log);
					
					// Compute new LB.
					LB = max(LB, best_cost + sum(pp.penalties));
					
					// If LB=UB, we know that UB is optimum, we can stop.
					if (epsilon_equal(LB, UB.duration))
					{
						early_stop = true;
						penalties = pp.penalties;
						clog << "Stop CG (LB = UB)" << endl;
						return false;
					}
					
					// Add the best solution to the RMP if it is not a feasible solution (otherwise it was already added).
					if (!R.empty()) spf.AddRoute(best);
					
					// Log iteration.
					cg_execution_log->iterations->push_back(iteration_log);
					
					// Keep best penalties if it is the last iteration.
					if (R.empty()) penalties = pp.penalties;
					return !R.empty();
				};
				
				// Run CG.
				auto cg_log = cg_solver.Solve(spf.formulation, {CGOption::IterationsInformation});
				output["CG"] = cg_log;
				
				// If gap was not closed, get best LB.
				if (!early_stop && cg_log.status == CGStatus::Optimum) LB = cg_log.incumbent_value;
				clog << "Finished CG in " << cg_log.time << "secs with LB: " << LB << endl;
			}
			bool found_opt = epsilon_equal(LB, UB.duration);
			if (found_opt) clog << "Optimality was closed in CG" << endl;
			
			// If optimum was not found, run exact algorithm.
			if (!found_opt)
			{
				CGExecutionLog dssr_log;
				MLBExecutionLog exact_log(true);
				run_dssr(vrp, NG, relaxation == "None" ? -1 : (dssr ? 5 : 0), penalties, UB, LB, &dssr_log, &exact_log);
				output["DSSR"] = dssr_log;
				output["Exact"] = exact_log;
			}
		}
		LPExecutionLog lb_log;
		lb_log.incumbent_value = LB;
		output["LB"] = lb_log;
		
		// Get best route.
		if (UB.duration != INFTY)
		{
			Route best = UB;
			clog << "Best solution: " << endl;
			clog << "\tpath: " << best.path << endl;
			clog << "\tt0: " << best.t0 << endl;
			clog << "\tduration: " << best.duration << endl;
			output["Best solution"] = VRPSolution(best.duration, {best});
			output["status"] = "Optimum";
		}
		output["status"] = "TimeLimitReached";
	}
	
	// Send JSON output to cout.
	cout << output << endl;
	return 0;
}