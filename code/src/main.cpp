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

using namespace std;
using namespace goc;
using namespace nlohmann;
using namespace tdtsptw;

VRPInstance reverse_instance(const VRPInstance& vrp)
{
	int n = vrp.D.VertexCount();
	
	VRPInstance rev;
	rev.D = vrp.D.Reverse();
	rev.o = vrp.d, rev.d = vrp.o;
	rev.T = vrp.T;
	for (Vertex v: vrp.D.Vertices()) rev.tw.push_back({rev.T-vrp.tw[v].right, rev.T-vrp.tw[v].left});
	for (Vertex v: vrp.D.Vertices()) rev.a.push_back(rev.tw[v].left);
	for (Vertex v: vrp.D.Vertices()) rev.b.push_back(rev.tw[v].right);
	rev.prec = Matrix<bool>(n,n, false);
	rev.prec_count = vector<int>(n, 0);
	rev.suc_count = vector<int>(n, 0);
	for (Vertex v: vrp.D.Vertices())
	{
		for (Vertex w: vrp.D.Vertices())
		{
			if (vrp.prec[w][v])
			{
				rev.prec[v][w] = true;
				rev.prec_count[w]++;
				rev.suc_count[v]++;
			}
		}
	}
	rev.LDT = rev.EAT = Matrix<double>(n, n);
	for (Vertex v: vrp.D.Vertices())
	{
		for (Vertex w: vrp.D.Vertices())
		{
			rev.EAT[v][w] = vrp.T - vrp.LDT[w][v];
			rev.LDT[v][w] = vrp.T - vrp.EAT[w][v];
		}
	}
	rev.arr = rev.tau = rev.dep = rev.pretau = Matrix<PWLFunction>(n, n);
	for (Vertex u: vrp.D.Vertices())
	{
		for (Vertex v: vrp.D.Successors(u))
		{
			// Compute reverse travel functions.
			rev.arr[v][u] = vrp.T - vrp.dep[u][v].Compose(vrp.T - PWLFunction::IdentityFunction({0.0, vrp.T}));
			rev.arr[v][u] = Min(PWLFunction::ConstantFunction(min(img(rev.arr[v][u])), {min(rev.tw[v]), min(dom(rev.arr[v][u]))}), rev.arr[v][u]);
			rev.tau[v][u] = rev.arr[v][u] - PWLFunction::IdentityFunction({0.0, vrp.T});
			rev.dep[v][u] = rev.arr[v][u].Inverse();
			rev.pretau[v][u] = PWLFunction::IdentityFunction(dom(rev.dep[v][u])) - rev.dep[v][u];
		}
	}
	// Add travel functions for (i, i) (for boundary reasons).
	for (Vertex u: rev.D.Vertices())
	{
		rev.tau[u][u] = rev.pretau[u][u] = PWLFunction::ConstantFunction(0.0, rev.tw[u]);
		rev.dep[u][u] = rev.arr[u][u] = PWLFunction::IdentityFunction(rev.tw[u]);
	}
	return rev;
}

int main(int argc, char** argv)
{
	json output; // STDOUT output will go into this JSON.
	
	simulate_runner_input("instances/guerriero_et_al_2014", "20_70_B_C6", "experiments/easy.json", "makespan");
	
	json experiment, instance, solutions;
	cin >> experiment >> instance >> solutions;
	
	// Parse experiment.
	Duration time_limit = Duration(value_or_default(experiment, "time_limit", 7200), DurationUnit::Seconds);
	string objective = value_or_default(experiment, "objective", "duration");
	
	// Show experiment details.
	clog << "Time limit: " << time_limit << "sec." << endl;
	clog << "Objective: " << objective << endl;
	
	// Set departing time from depot equal to 0 if makespan objective.
	if (objective == "makespan") instance["time_windows"][0] = Interval(0, 0);
	
	// Preprocess instance JSON.
	preprocess_travel_times(instance);
	preprocess_waiting_times(instance);
	preprocess_time_windows(instance);
	
	// Parse instance.
	clog << "Parsing instance..." << endl;
	VRPInstance vrp = instance;
	VRPInstance back_vrp = reverse_instance(vrp);
	
	// Get UB.
	double LB = 0.0;
	vector<Vertex> P = {vrp.o};
	Route UB = initial_heuristic(vrp, P, create_bitset<MAX_N>({vrp.o}), vrp.tw[vrp.o].left);
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
		
		// Run subgradient.
		CGExecutionLog subgradient_log;
		auto sg_routes = subgradient(vrp, NG, 10, UB, LB, &subgradient_log);
		output["Subgradient"] = subgradient_log;
		
		// Solve CG to obtain best penalties.
		if (epsilon_smaller(LB, UB.duration))
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
			
			vector<double> penalties(vrp.D.VertexCount(), 0.0); // Keep best set of penalties.
			bool early_stop = false; // If any other termination condition is met, early_stop is set to true.
			cg_solver.pricing_function = [&](const vector<double>& duals, double incumbent_value, Duration time_limit,
											 CGExecutionLog* cg_execution_log) {
				auto pp = spf.InterpretDuals(duals);
				MLBExecutionLog iteration_log(true);
				Route best;
				double best_cost;
				auto R = run_ng(vrp, NG, pp.penalties, UB.duration, &best, &best_cost, &iteration_log, nullptr);
				
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
			bool found_opt = epsilon_equal(LB, UB.duration);
			if (found_opt) clog << "Optimality was closed in CG" << endl;
			
			// If optimum was not found, run exact algorithm.
			if (!found_opt)
			{
				clog << "Obtaining bound labels..." << endl;
				Route best;
				double best_cost;
				MLBExecutionLog bound_log(true);
				BoundingStructure B(&vrp, &NG, penalties);
				run_ng(vrp, NG, penalties, UB.duration, &best, &best_cost, &bound_log, &B);
				
				clog << "Running exact algorithm..." << endl;
				MLBExecutionLog exact_log(true);
				UB = Route(reverse(UB.path), vrp.T - UB.t0, UB.duration);
				UB = run_exact(back_vrp, NG, &B, penalties, UB, LB, &exact_log);
				UB = Route(reverse(UB.path), vrp.T - UB.t0, UB.duration);
				output["Exact"] = exact_log;
			}
		}
		
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