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

namespace
{
VRPInstance reverse_instance(const VRPInstance& vrp)
{
	int n = vrp.D.VertexCount();
	
	VRPInstance rev;
	rev.D = vrp.D.Reverse();
	rev.o = vrp.d, rev.d = vrp.o;
	rev.T = vrp.T;
	for (Vertex v: vrp.D.Vertices()) rev.tw.push_back({-vrp.tw[v].right, -vrp.tw[v].left});
	for (Vertex v: vrp.D.Vertices()) rev.a.push_back(rev.tw[v].left);
	for (Vertex v: vrp.D.Vertices()) rev.b.push_back(rev.tw[v].right);
	rev.prec = Matrix<bool>(n, n, false);
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
			rev.EAT[v][w] = -vrp.LDT[w][v];
			rev.LDT[v][w] = -vrp.EAT[w][v];
		}
	}
	rev.arr = rev.tau = rev.dep = rev.pretau = Matrix<PWLFunction>(n, n);
	for (Vertex u: vrp.D.Vertices())
	{
		for (Vertex v: vrp.D.Successors(u))
		{
			// Compute reverse travel functions.
			rev.tau[v][u] = vrp.pretau[u][v].Compose(PWLFunction::IdentityFunction({-vrp.T, 0.0}) * -1);
			auto init_piece = LinearFunction({-vrp.T, rev.tau[v][u](min(dom(rev.tau[v][u]))) + min(dom(rev.tau[v][u])) + vrp.T}, {min(dom(rev.tau[v][u])), rev.tau[v][u](min(dom(rev.tau[v][u])))});
			rev.tau[v][u] = Min(rev.tau[v][u], PWLFunction({init_piece}));
			rev.arr[v][u] = rev.tau[v][u] + PWLFunction::IdentityFunction(rev.tau[v][u].Domain());
			rev.dep[v][u] = rev.arr[v][u].Inverse();
			rev.pretau[v][u] = PWLFunction::IdentityFunction(rev.dep[v][u].Domain()) - rev.dep[v][u];
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
}

int main(int argc, char** argv)
{
	try
	{
		json output; // STDOUT output will go into this JSON.
		
		simulate_runner_input("instances/td-ascheuer", "rbg048a", "experiments/ascheuer.json", "CG-NGLTI-DNA");
		
		json experiment, instance, solutions;
		cin >> experiment >> instance >> solutions;
		
		if (!has_key(instance, "time_windows"))
		{
			int n = instance["digraph"]["vertex_count"];
			instance["time_windows"] = vector<Interval>(n);
			for (int i = 0; i < n; ++i) instance["time_windows"][i] = instance["horizon"];
		}
		
		// Parse experiment.
		Duration time_limit = Duration(value_or_default(experiment, "time_limit", 7200), DurationUnit::Seconds);
		string objective = value_or_default(experiment, "objective", "duration");
		string relaxation = value_or_default(experiment, "relaxation", "NGLTD");
		bool colgen = value_or_default(experiment, "colgen", true);
		bool dssr = value_or_default(experiment, "dssr", false);
		bool time_independent = value_or_default(experiment, "time_independent", false);
		bool bidirectional_dna = value_or_default(experiment, "bidirectional_dna", false);
		bool bidirectional_cg = value_or_default(experiment, "bidirectional_cg", false);

		// Time-independentize instance.
		if (time_independent)
			for (auto& cluster_row: instance["cluster_speeds"])
				for (auto& speed: cluster_row)
					speed = 1.0;
		
		// Show experiment details.
		clog << "Time limit: " << time_limit << "sec." << endl;
		clog << "Objective: " << objective << endl;
		clog << "Relaxation: " << relaxation << endl;
		clog << "Colgen: " << colgen << endl;
		clog << "DSSR: " << dssr << endl;
		clog << "Bidirectional DNA: " << bidirectional_dna << endl;
		clog << "Bidirectional CG: " << bidirectional_cg << endl;
		
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
			auto rvrp = reverse_instance(vrp);
			NGStructure NG(vrp, 3);
			NGStructure rNG(rvrp, NG.N, reverse(NG.L), NG.delta);

			vector<double> penalties(vrp.D.VertexCount(), 0.0); // Keep best set of penalties.

			// Run subgradient.
			Stopwatch rolex(true);
			vector<Route> sg_routes;
//			CGExecutionLog subgradient_log;
//			sg_routes = subgradient(vrp, NG, relaxation == "NGL", 10, UB, LB, penalties, &subgradient_log);
//			output["Subgradient"] = subgradient_log;
			
			// Solve CG to obtain best penalties.
			if (epsilon_smaller(LB, UB.duration))
			{
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
					cg_solver.time_limit = 1000.0_sec;
					cg_solver.lp_solver = &lp_solver;
					cg_solver.screen_output = &clog;
					
					bool early_stop = false; // If any other termination condition is met, early_stop is set to true.
					cg_solver.pricing_function = [&](const vector<double>& duals, double incumbent_value,
													 Duration time_limit, CGExecutionLog* cg_execution_log) {
						if (!colgen) return false;
						auto pp = spf.InterpretDuals(duals);
						Route best;
						double best_cost;
						if (relaxation == "NGLTI2RES")
						{
							MLBExecutionLog iteration_log(true);
							run_ngl2res(vrp, NG, pp.penalties, UB.duration, &best, &best_cost, &iteration_log);
							cg_execution_log->iterations->push_back(iteration_log);
						}
						else if (relaxation == "NGLTI")
						{
							if (!bidirectional_cg)
							{
								MLBExecutionLog iteration_log(true);
								run_nglti(vrp, NG, pp.penalties, UB.duration, best, best_cost, &iteration_log);
								cg_execution_log->iterations->push_back(iteration_log);
							}
							else
							{
								BLBExecutionLog blb_log(true);
								bidirectional_run_nglti(vrp, rvrp, NG, rNG, pp.penalties, UB.duration, best, best_cost, &blb_log);
							}
						}
						else if (relaxation == "NGLTD")
						{
							if (!bidirectional_cg)
							{
								MLBExecutionLog iteration_log(true);
								best = run_ngltd(vrp, NG, pp.penalties, &iteration_log, nullptr, LB, time_limit);
								best_cost = best.duration - sum<Vertex>(best.path, [&](Vertex v) { return pp.penalties[v]; });
								cg_execution_log->iterations->push_back(iteration_log);
							}
							else
							{
								BLBExecutionLog iteration_log(true);
								best = run_ngltd_bidirectional(vrp, rvrp, NG, rNG, pp.penalties, &iteration_log, LB, time_limit);
								best_cost = best.duration - sum<Vertex>(best.path, [&](Vertex v) { return pp.penalties[v]; });
								cg_execution_log->iterations->push_back(iteration_log);
							}
						}
						
						// Compute new LB.
						if (best_cost + sum(pp.penalties) >= LB)
						{
							LB = best_cost + sum(pp.penalties);
							penalties = pp.penalties;
						}
						
						// If LB=UB, we know that UB is optimum, we can stop.
						if (epsilon_equal(LB, UB.duration))
						{
							early_stop = true;
							clog << "Stop CG (LB = UB)" << endl;
							return false;
						}
						
						// Add the best solution to the RMP if it is not a feasible solution (otherwise it was already added).
						if (epsilon_smaller(best_cost, 0.0)) spf.AddRoute(best);

						return epsilon_smaller(best_cost, 0.0);
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

				clog << "Penalties: " << penalties << endl;
				
				if (dssr && !found_opt)
				{
					clog << "Running DNA to improve bounds." << endl;
					CGExecutionLog dna_log;
					auto R = run_dna(vrp, rvrp, NG, rNG, penalties, &dna_log, LB, 200.0_sec, bidirectional_dna);
					if (R.duration != INFTY)
					{
						UB = R;
						clog << "\tFound solution " << UB.path << " " << UB.duration << endl;
						dna_log.status = CGStatus::Optimum;
						found_opt = true;
					}
					output["DNA"] = dna_log;
				}
				
				// If optimum was not found, run exact algorithm.
				if (!found_opt)
				{
					MLBExecutionLog log_ngl(true);
					Bounding B(vrp, NG, penalties);
					if (relaxation != "None")
					{
						Route R_NG = run_ngltd(vrp, NG, penalties, &log_ngl, &B, LB, time_limit - rolex.Peek());
						output["NGL-TD"] = log_ngl;
						clog << "NGL-TD:   " << LB << "\t" << log_ngl.time << "\t" << log_ngl.processed_count << "\t"
							 << log_ngl.enumerated_count << endl;
					}
					
					MLBExecutionLog log(true);
					auto r = run_exact_piecewise(rvrp, reverse(NG.L), penalties, LB, UB.duration, &log,
												 relaxation == "None" ? nullptr : &B, time_limit - rolex.Peek());
					clog << "Exact: " << r.duration << "\t" << log.time << "\t" << log.processed_count << "\t"
						 << log.enumerated_count << endl;
					output["Exact"] = log;
					if (r.duration < UB.duration) UB = vrp.BestDurationRoute(reverse(r.path));
				}
			}
			LPExecutionLog lb_log;
			lb_log.incumbent_value = LB;
			lb_log.time = rolex.Peek();
			lb_log.status = rolex.Peek() > time_limit ? LPStatus::TimeLimitReached : LPStatus::Optimum;
			output["General"] = lb_log;
			
			// Get best route.
			if (UB.duration != INFTY)
			{
				Route best = UB;
				clog << "Best solution: " << endl;
				clog << "\tpath: " << best.path << endl;
				clog << "\tt0: " << best.t0 << endl;
				clog << "\tduration: " << best.duration << endl;
				output["Best solution"] = VRPSolution(best.duration, {best});
			}
			output["status"] = rolex.Peek() < time_limit ? "Optimum" : "TimeLimitReached";
		}
		
		// Send JSON output to cout.
		cout << output << endl;
	}
	catch (std::bad_alloc& e)
	{
		return 3;
	}
	return 0;
}