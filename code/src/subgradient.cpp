//
// Created by Gonzalo Lera Romero on 16/09/2019.
//

#include "subgradient.h"

using namespace std;
using namespace goc;
using namespace nlohmann;

namespace tdtsptw
{
namespace
{
// Indicates if the path is fesasible tour.
// @returns: if the path is feasible (with respect to resources and elementarity).
bool is_feasible_solution(const VRPInstance& vrp, const GraphPath& path)
{
	VertexSet V;
	for (auto &v: path)
	{
		if (V.test(v)) return false;
		V.set(v);
	}
	return vrp.ReadyTime(path, 0.0) != INFTY;
}
}

vector<Route> subgradient(const VRPInstance& vrp, const NGStructure& NG, bool use_td_relaxation, int max_iter, Route& UB, double& LB, CGExecutionLog* log)
{
	vector<double> lambda(vrp.D.VertexCount(), 0.0);
	map<GraphPath, Route> route_set; // To avoid returning multiple times the same solution, we index them by the route.
	log->iteration_count = 0;
	log->iterations = vector<json>();
	Stopwatch rolex(true);
	for (int t = 0; t < max_iter; ++t)
	{
		clog << "Iteration " << (t+1) << endl;
		// Run NG pricing problem.
		double Lambda = sum(lambda);
		MLBExecutionLog it_log(true);
		Route best;
		double best_LB;
		vector<Route> Routes;
		if (!use_td_relaxation) Routes = run_ng_td(vrp, NG, lambda, UB.duration, &best, &best_LB, &it_log, nullptr);
		else Routes = run_ng(vrp, NG, lambda, UB.duration, &best, &best_LB, &it_log, nullptr);
		LB = max(LB, best_LB + Lambda);
		
		// Log iteration information.
		log->iteration_count++;
		log->iterations->push_back(it_log);
		
		if (epsilon_equal(LB, UB.duration))
		{
			clog << "\tOptimality gap closed in Subgradient Algorithm." << endl;
			return {};
		}
		
		// Compute new penalties.
		vector<int> delta(vrp.D.VertexCount(), 0);
		for (Vertex v: best.path) delta[v]++;
		if (all_of(best.path.begin(), best.path.end(), [&] (Vertex v) { return delta[v] == 1; })) break;
		vector<double> lambda_prime = lambda;
		for (Vertex j: vrp.D.Vertices())
			lambda_prime[j] = (lambda[j] - (LB - 1.2 * (LB + sum<Vertex>(vrp.D.Vertices(), [&] (Vertex l) { return delta[l] * lambda[l]; })))) * (2.0 - 2.0 * delta[j]) / (sum<Vertex>(vrp.D.Vertices(), [&] (Vertex l) { return (2.0 - 2.0 * delta[l]) * (2.0 - 2.0 * delta[l]); }));
		lambda = lambda_prime;
		for (auto& r: Routes) route_set[r.path] = r;
		route_set[best.path] = best;
	}
	clog << "LB subgradient: " << LB << endl;
	log->incumbent_value = LB;
	log->time = rolex.Peek();
	log->status = CGStatus::Optimum;
	vector<Route> r;
	for (auto& route: route_set) r.push_back(route.second);
	return r;
}
} // namespace tdtsptw