//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "labeling.h"
#include "pwl_domination_function.h"

using namespace std;
using namespace goc;

namespace tdtsptw
{
Route initial_heuristic(const VRPInstance& vrp, vector<Vertex>& P, VertexSet S, double t)
{
	int n = vrp.D.VertexCount();
	Vertex u = P.back();
	if (P.size() == n)
	{
		if (u == vrp.d) return Route(P, 0.0, t);
		else return Route({}, 0.0, INFTY);
	}
	for (Vertex w: vrp.D.Vertices()) if (!S.test(w) && epsilon_smaller(vrp.LDT[u][w], t)) return Route({}, 0.0, INFTY);
	
	for (Vertex v: vrp.D.Successors(P.back()))
	{
		if (S.test(v)) continue;
		P.push_back(v);
		Route r = initial_heuristic(vrp, P, unite(S, {v}), vrp.ArrivalTime({u,v}, t));
		if (r.duration != INFTY) return r;
		P.pop_back();
	}
	return Route({}, 0.0, INFTY);
}

Route run_labeling(const VRPInstance& vrp, const Duration& time_limit, MLBExecutionLog* log, bool t0_is_zero)
{
	int n = vrp.D.VertexCount();
	
	*log = MLBExecutionLog(true);
	Stopwatch rolex(true), rolex2(false), rolex3(false);
	auto comp = [&] (Label* l, Label* m) { return l->cost > m->cost; };
	priority_queue<Label*, vector<Label*>, decltype(comp)> q(comp);
	q.push(new Label{nullptr, vrp.o, 1, create_bitset<MAX_N>({vrp.o}), PWLFunction::ConstantFunction(0.0, {vrp.tw[vrp.o].left, t0_is_zero ? vrp.tw[vrp.o].left : vrp.tw[vrp.o].right}), vrp.tw[vrp.o].left, 0.0, 0.0});

	Matrix<unordered_map<VertexSet, vector<Label*>>> D(n, n+1);
	Route best({}, 0.0, INFTY);
	while (!q.empty())
	{
		if (rolex.Peek() >= time_limit) { log->status = MLBStatus::TimeLimitReached; break; }

		rolex2.Reset().Resume();
		Label* l = q.top();
		q.pop();
		*log->queuing_time += rolex2.Peek();

		// If l is a tour, then check if it is the best solution.
		if (l->k == n)
		{
			if (l->cost < best.duration) best = Route(l->Path(), l->D.PreValue(l->cost)-l->cost, l->cost);
			delete l;
			continue;
		}

		// Bounding step.
		if (epsilon_bigger_equal(l->cost, best.duration)) { delete l; continue; }

		// Domination step.
		rolex2.Reset().Resume();
		PWLDominationFunction l_D = l->D;
		for (Label* m: D[l->v][l->k][l->S])
		{
			if (epsilon_smaller(max(img(l->D)), m->cost)) break;
			if (l_D.DominatePieces(m->D)) break;
		}
		bool is_dominated = l_D.Empty();
		if (!is_dominated) { l->D = (PWLFunction)l_D; l->t = min(dom(l->D)); l->cost = min(img(l->D)); }
		*log->domination_time += rolex2.Peek();
		if (is_dominated) { *log->positive_domination_time += rolex2.Peek(); log->dominated_count++; }
		else { *log->negative_domination_time += rolex2.Peek(); log->processed_count++; }
		if (is_dominated) { delete l; continue; } // If D(l) == \emptyset, then l is dominated.

		// l is not dominated, add it to the domination structure.
		rolex2.Reset().Resume();
		insert_sorted(D[l->v][l->k][l->S], l, [] (Label* l1, Label* l2) { return l1->cost < l2->cost; });
		*log->process_time += rolex2.Peek();
		stretch_to_size(*log->count_by_length, l->k+1, 0);
		(*log->count_by_length)[l->k]++;

		// Extend label l.
		rolex2.Reset().Resume();
		for (Vertex w: vrp.D.Successors(l->v))
		{
			// Check feasibility.
			if (l->S.test(w)) continue;
			if (epsilon_bigger(l->t, vrp.LDT[l->v][w])) continue;
			TimeUnit t_w = vrp.ArrivalTime({l->v, w}, l->t);
			bool any_unreachable = false;
			for (Vertex u: vrp.D.Vertices()) if (!l->S.test(u) && u != w && epsilon_smaller(vrp.LDT[w][u], t_w)) { any_unreachable = true; break; }
			if (any_unreachable) continue;
			PWLFunction D_w = epsilon_smaller(max(dom(l->D)), min(img(vrp.dep[l->v][w])))
						   ? PWLFunction::ConstantFunction(l->D(max(dom(l->D))) + min(vrp.tw[w]) - max(dom(l->D)), {min(vrp.tw[w]), min(vrp.tw[w])})
						   : (l->D + vrp.tau[l->v][w]).Compose(vrp.dep[l->v][w]);
			if (D_w.Empty()) continue;
			rolex2.Pause();
			rolex3.Reset().Resume();
			q.push(new Label(l, w, l->k+1, unite(l->S, {w}), D_w, min(dom(D_w)), min(img(D_w)), 0.0));
			*log->queuing_time += rolex3.Peek();
			rolex2.Resume();
			log->extended_count++;
		}
		*log->extension_time += rolex2.Peek();
	}
	log->time = rolex.Peek();
	if (log->status == MLBStatus::DidNotStart) log->status = MLBStatus::Finished;

	return best;
}

vector<VertexSet> NGSet(const VRPInstance& vrp, int delta)
{
	vector<VertexSet> N(vrp.D.VertexCount());
	for (Vertex i: vrp.D.Vertices())
	{
		vector<pair<double, Vertex>> N_by_dist;
		for (Vertex j: vrp.D.Vertices())
		{
			if (i == j || j == vrp.o || j == vrp.d) continue;
			if (epsilon_bigger(vrp.EAT[j][i], vrp.LDT[i][j])) continue;
			N_by_dist.push_back({vrp.TravelTime({j,i}, vrp.tw[j].left) , j});
		}
		sort(N_by_dist.begin(), N_by_dist.end());
		for (int k = 0; k < min(delta, (int)N_by_dist.size()); ++k) N[i].set(N_by_dist[k].second);
	}
	return N;
}

vector<Route> run_ng_labeling(const VRPInstance& vrp, const Duration& time_limit, MLBExecutionLog* log, bool t0_is_zero, const vector<double>& penalties)
{
	int n = vrp.D.VertexCount();
	
	// Compute neighbours.
	vector<VertexSet> N = NGSet(vrp, 3);
	
	*log = MLBExecutionLog(true);
	Stopwatch rolex(true), rolex2(false), rolex3(false);
	auto comp = [&] (Label* l, Label* m) { return make_tuple(l->k, l->cost) > make_tuple(m->k, m->cost); };
	priority_queue<Label*, vector<Label*>, decltype(comp)> q(comp);
	q.push(new Label{nullptr, vrp.o, 1, create_bitset<MAX_N>({vrp.o}), PWLFunction::ConstantFunction(0.0, {vrp.tw[vrp.o].left, t0_is_zero ? vrp.tw[vrp.o].left : vrp.tw[vrp.o].right}), vrp.tw[vrp.o].left, 0.0, penalties[vrp.o]});
	
	Matrix<vector<Label*>> D(n, n+1);
	Route best({}, 0.0, INFTY);
	double best_cost = INFTY;
	vector<Route> solutions;
	int count = 0, count2 = 0;
	while (!q.empty())
	{
		if (rolex.Peek() >= time_limit) { log->status = MLBStatus::TimeLimitReached; break; }
		
		rolex2.Reset().Resume();
		Label* l = q.top();
		q.pop();
		*log->queuing_time += rolex2.Peek();
		
		// If l is a tour, then check if it is the best solution.
		if (l->k == n && l->v == vrp.d)
		{
			if (l->cost < best_cost)
			{
				best = Route(l->Path(), l->D.PreValue(min(img(l->D)))-min(img(l->D)), min(img(l->D)));
				best_cost = l->cost;
			}
			if (epsilon_smaller(l->cost, 0.0))
			{
				solutions.push_back(Route(l->Path(), l->D.PreValue(min(img(l->D)))-min(img(l->D)), min(img(l->D))));
				if (solutions.size() == 3000) break;
			}
			delete l;
			continue;
		}
		
		// Domination step.
		rolex2.Reset().Resume();
		PWLDominationFunction l_D = l->D;
		
		bool is_dominated = false;
		for (Label* m: D[l->v][l->k])
		{
			count++;
			if (epsilon_smaller(max(img(l->D))-l->P, m->cost)) break;
			if (!is_subset(m->S, l->S)) continue;
			count2++;
//			if (is_dominated = (epsilon_smaller_equal(m->t, l->t) && epsilon_bigger_equal(m->P, l->P))) break;
			if (is_dominated = l_D.DominatePieces(m->D, l->P-m->P)) break;
		}
//		if (count % 10000 == 0) clog << count << " vs " << count2 << endl;
		if (!is_dominated) { l->D = (PWLFunction)l_D; l->t = min(dom(l->D)); l->cost = min(img(l->D))-l->P; }
		*log->domination_time += rolex2.Peek();
		if (is_dominated) { *log->positive_domination_time += rolex2.Peek(); log->dominated_count++; }
		else { *log->negative_domination_time += rolex2.Peek(); log->processed_count++; }
		if (is_dominated) { delete l; continue; } // If D(l) == \emptyset, then l is dominated.
		
		// l is not dominated, add it to the domination structure.
		rolex2.Reset().Resume();
		insert_sorted(D[l->v][l->k], l, [] (Label* l1, Label* l2) { return l1->cost < l2->cost; });
		*log->process_time += rolex2.Peek();
		stretch_to_size(*log->count_by_length, l->k+1, 0);
		(*log->count_by_length)[l->k]++;
		
		// Extend label l.
		rolex2.Reset().Resume();
		for (Vertex w: vrp.D.Successors(l->v))
		{
			// Check feasibility.
			if (l->S.test(w)) continue;
			if (epsilon_bigger(l->t, vrp.LDT[l->v][w])) continue;
			if (vrp.prec_count[w] > l->k) continue;
			if (l->k == n-1 && w != vrp.d) continue;
			if (l->k < n-1 && w == vrp.d) continue;
			
			TimeUnit t_w = vrp.ArrivalTime({l->v, w}, l->t);
			int unreachable_count = 0;
			for (Vertex u: vrp.D.Vertices()) if (!l->S.test(u) && u != w && epsilon_smaller(vrp.LDT[w][u], t_w)) unreachable_count++;
			if (unreachable_count > l->k) continue;
			PWLFunction D_w = epsilon_smaller(max(dom(l->D)), min(img(vrp.dep[l->v][w])))
							  ? PWLFunction::ConstantFunction(l->D(max(dom(l->D))) + min(vrp.tw[w]) - max(dom(l->D)), {min(vrp.tw[w]), min(vrp.tw[w])})
							  : (l->D + vrp.tau[l->v][w]).Compose(vrp.dep[l->v][w]);
			if (D_w.Empty()) continue;
			rolex2.Pause();
			rolex3.Reset().Resume();
			q.push(new Label(l, w, l->k+1, unite(intersection(l->S, N[w]), {w}), D_w, min(dom(D_w)), min(img(D_w))-l->P-penalties[w], l->P + penalties[w]));
			*log->queuing_time += rolex3.Peek();
			rolex2.Resume();
			log->extended_count++;
		}
		*log->extension_time += rolex2.Peek();
	}
	log->time = rolex.Peek();
	if (log->status == MLBStatus::DidNotStart) log->status = MLBStatus::Finished;
	clog << "Cost: " << best_cost << endl;
	clog << "Best bound: " << best_cost + sum<Vertex>(vrp.D.Vertices(), [&] (Vertex v) { return penalties[v]; }) << endl;
	return solutions;
}
} // namespace tdtsptw