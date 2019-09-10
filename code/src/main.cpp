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

vector<Route> subgradient(const VRPInstance& vrp, CGExecutionLog* log, bool is_t0_zero)
{
	double LB = 0.0;
	vector<double> lambda(vrp.D.VertexCount(), 0.0);
	vector<double> lambda_star = lambda;
	auto& V = vrp.D.Vertices();
	map<GraphPath, Route> R;
	log->iteration_count = 0;
	log->iterations = vector<json>();
	for (int t = 0; t < 10; ++t)
	{
		clog << "Iteration " << (t+1) << endl;
		clog << "\tlambda: " << lambda << endl;
		double Lambda = sum<Vertex>(vrp.D.Vertices(), [&] (Vertex v) { return lambda[v]; });
		MLBExecutionLog it_log(true);
		Route best;
		run_ng_labeling(vrp, 100.0_sec, &it_log, is_t0_zero, lambda, false, &best);
		log->iteration_count++;
		log->iterations->push_back(it_log);
		double LB_best = best.duration - sum<Vertex>(V, [&] (Vertex v) { return lambda[v]; }) + Lambda;
		clog << "\tBest: " << best << endl;
		clog << "\tLB_Best: " << LB_best << endl;
		clog << "\tLB: " << LB << endl;
		if (LB_best > LB)
		{
			LB = LB_best;
			lambda_star = lambda;
		}
		vector<int> delta(vrp.D.VertexCount(), 0);
		for (Vertex v: best.path) delta[v]++;
		vector<double> lambda_prime = lambda;
		for (Vertex j: vrp.D.Vertices())
			lambda_prime[j] = (lambda[j] - (LB - 1.2 * (LB + sum<Vertex>(V, [&] (Vertex l) { return delta[l] * lambda[l]; })))) * (2.0 - 2.0 * delta[j]) / (sum<Vertex>(V, [&] (Vertex l) { return (2.0 - 2.0 * delta[l]) * (2.0 - 2.0 * delta[l]); }));
		lambda = lambda_prime;
		R[best.path] = best;
	}
	vector<Route> r;
	for (auto& route: R) r.push_back(route.second);
	return r;
}

int main(int argc, char** argv)
{
	json output; // STDOUT output will go into this JSON.
	
	simulate_runner_input("instances/afg_1995", "rbg050a", "experiments/afg.json", "duration");
	
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

	CGExecutionLog subgradient_log;
//	auto SR = subgradient(vrp, &subgradient_log, objective == "makespan");
	SPF spf(vrp.D);
//	for (auto& r: SR) spf.AddRoute(r);
	vector<Vertex> P = {vrp.o};
	Route UB = initial_heuristic(vrp, P, create_bitset<MAX_N>({vrp.o}), vrp.tw[vrp.o].left);
//	Route UB({0, 1, 3, 4, 11, 9, 6, 2, 16, 10, 12, 19, 17, 18, 13, 14, 15, 8, 7, 5, 20}, 238, 4484);
//	Route UB({0, 1, 2, 3, 4, 5, 7, 6, 9, 8, 10, 12, 11, 13, 14, 15, 16, 17, 18, 20, 19, 22, 21, 24, 23, 26, 27, 25, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 42, 43, 50, 41, 49, 46, 45, 48, 51, 44, 52, 47, 54, 53, 55, 56}, 562, 6367);
	clog << "UB: " << UB << endl;
	spf.AddRoute(UB);
	CGSolver cg_solver;
	LPSolver lp_solver;
	cg_solver.time_limit = 100.0_sec;
	cg_solver.lp_solver = &lp_solver;
	cg_solver.screen_output = &clog;
	bool relaxation = true;
	cg_solver.pricing_function = [&] (const vector<double>& duals, double incumbent_value, Duration time_limit, CGExecutionLog* cg_execution_log)
	{
		MLBExecutionLog iteration_log(true);
		auto pp = spf.InterpretDuals(duals);
//		auto R = run_ng_labeling_non_hamiltonian(vrp, time_limit, &iteration_log, objective == "makespan", pp.penalties, pp.sigma, relaxation);
		auto R = run_ng_labeling(vrp, time_limit, &iteration_log, objective == "makespan", spf.InterpretDuals(duals).penalties, relaxation);
		for (auto& r: R) spf.AddRoute(r);
		cg_execution_log->iterations->push_back(iteration_log);
		if (R.empty() && relaxation) { relaxation = false; clog << "Relaxation off" << endl; return true; }
		return !R.empty();
	};
	CGExecutionLog cg_log;
//	cg_log = cg_solver.Solve(spf.formulation, {CGOption::IterationsInformation});
//	clog << "Finished CG in " << cg_log.time << "secs." << endl;
	
	// Solve Exact.
	MLBExecutionLog exact_log(true);
	Route best;
	best = run_labeling(vrp, time_limit, &exact_log, objective == "makespan");
	
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