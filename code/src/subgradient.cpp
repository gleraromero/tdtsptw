//
// Created by Gonzalo Lera Romero on 16/09/2019.
//

#include "subgradient.h"

#include "lbl_exact.h"

using namespace std;
using namespace goc;
using namespace nlohmann;

namespace tdtsptw
{
vector<Route> subgradient(const VRPInstance& vrp, const NGStructure& NG, const string& relaxation, int max_iter, Route& UB, double& LB, vector<double>& penalties, CGExecutionLog* log)
{
	vector<double> lambda(vrp.D.VertexCount(), 0.0);
	map<GraphPath, Route> route_set; // To avoid returning multiple times the same solution, we index them by the route.
	log->iteration_count = 0;
	log->iterations = vector<json>();
	Stopwatch rolex(true);
	for (int t = 0; t < 1; ++t)
	{
		clog << "Iteration " << (t + 1) << endl;
		// Run NG pricing problem.
		double Lambda = sum(lambda);
		Route r;
		double cost_r;
		double LB_r;
		if (relaxation == "NGLTI2RES")
		{
			MLBExecutionLog iteration_log(true);
			run_ngl2res(vrp, NG, lambda, UB.duration, &r, &cost_r, &iteration_log);
			log->iterations->push_back(iteration_log);
		} else if (relaxation == "NGLTI")
		{
			MLBExecutionLog iteration_log(true);
			run_nglti(vrp, NG, lambda, UB.duration, r, cost_r, &iteration_log);
			log->iterations->push_back(iteration_log);
		} else if (relaxation == "NGLTD")
		{
			MLBExecutionLog iteration_log(true);
			r = run_ngltd(vrp, NG, lambda, &iteration_log, nullptr, LB, 1000.0_sec);
			log->iterations->push_back(iteration_log);
			cost_r = r.duration - sum<Vertex>(r.path, [&](Vertex v) { return lambda[v]; });
		}
		LB_r = cost_r + Lambda;
		clog << "New cost: " << cost_r << " " << Lambda << " with LB: " << LB_r << endl;
		if (cost_r + Lambda > LB)
		{
			LB = cost_r + Lambda;
			penalties = lambda;
		}

		// Log iteration information.
		log->iteration_count++;

		if (epsilon_equal(LB, UB.duration))
		{
			clog << "\tOptimality gap closed in Initial NG." << endl;
			return {};
		}

		// Compute new penalties.
		vector<int> delta(vrp.D.VertexCount(), 0);
		for (Vertex v: r.path) delta[v]++;
		if (all_of(r.path.begin(), r.path.end(), [&](Vertex v) { return delta[v] == 1; }))
		{
			clog << "\tFound elementary route in Initial NG." << endl;
			break;
		}

		double norm_square_gk = sum<Vertex>(vrp.D.Vertices(),
											[&](Vertex v) { return (1.0 * delta[v] - 1.0) * (1.0 * delta[v] - 1.0); });
		double step_size = (0.2 * LB_r) / norm_square_gk;
		for (Vertex j: vrp.D.Vertices()) lambda[j] += step_size * (delta[j] - 1.0);
		route_set[r.path] = r;
	}
	clog << "LB initial NG: " << LB << endl;
	log->incumbent_value = LB;
	log->time = rolex.Peek();
	log->status = CGStatus::Optimum;
	vector<Route> r;
//	for (auto &route: route_set) r.push_back(route.second);
	return r;
}
} // namespace tdtsptw