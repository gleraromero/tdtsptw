//
// Created by Gonzalo Lera Romero on 16/09/2019.
//

#include "initial_lb.h"

#include "labeling.h"

using namespace std;
using namespace goc;
using namespace nlohmann;

namespace tdtsptw
{
vector<Route> initial_lb(const VRPInstance& vrp, const NGStructure& NG, const string& relaxation, Route& UB, double& LB, CGExecutionLog* log)
{
	vector<double> lambda(vrp.D.VertexCount(), 0.0);
	log->iteration_count = 0;
	log->iterations = vector<json>();
	Stopwatch rolex(true);
	// Run NG pricing problem.
	double Lambda = sum(lambda);
	Route r;
	double cost_r;
	double LB_r;
	if (relaxation == "NGLTI")
	{
		MLBExecutionLog iteration_log(true);
		run_nglti(vrp, NG, lambda, UB.duration, r, cost_r, &iteration_log);
		log->iterations->push_back(iteration_log);
	}
	else if (relaxation == "NGLTD")
	{
		MLBExecutionLog iteration_log(true);
		r = run_ngltd(vrp, NG, lambda, &iteration_log, nullptr, LB, 1000.0_sec);
		log->iterations->push_back(iteration_log);
		cost_r = r.duration - sum<Vertex>(r.path, [&](Vertex v) { return lambda[v]; });
	}
	if (cost_r + Lambda > LB)
	{
		LB = cost_r + Lambda;
	}

	// Log iteration information.
	log->iteration_count++;

	if (epsilon_equal(LB, UB.duration))
	{
		clog << "\tOptimality gap closed in Initial NG." << endl;
		return {};
	}

	clog << "LB initial NG: " << LB << endl;
	log->incumbent_value = LB;
	log->time = rolex.Peek();
	log->status = CGStatus::Optimum;
	return {};
}
} // namespace tdtsptw