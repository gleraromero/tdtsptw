//
// Created by Gonzalo Lera Romero on 07/05/2020.
//

#include "exact_solver.h"

using namespace std;
using namespace goc;
using namespace nlohmann;

namespace tdtsptw
{
// Reconstructs a path from:
// 	vrp: instance to work with.
//	penalties: penalties of the vertices.
//	L: map (v, S) -> LabelSequence of non-dominated labels.
goc::GraphPath reconstruct_path(const VRPInstance& vrp,
								const vector<double>& penalties,
								vector<std::unordered_map<VertexSet, LabelSequenceTD>>& L)
{
	// Variables that contain information about the path to be built.
	Vertex v = vrp.d;
	VertexSet S;
	for (Vertex w: vrp.D.Vertices()) S.set(w); // S = V.
	double cost, time;
	L[v][S].LowestCostArrival(&time, &cost);
	double route_duration = cost + sum(penalties);

	// Get path.
	int n = vrp.D.VertexCount();
	goc::GraphPath path = {v};
	for (int k = n-1; k > 0; --k)
	{
		bool found_next = false;
		auto Su = S;
		Su.reset(v);
		for (goc::Vertex u: vrp.D.Vertices())
		{
			if (u == v) continue;
			if (!S.test(u)) continue;
			if (!includes_key(L[u], Su)) continue;
			auto& l = L[u][Su];
			double time_u = vrp.DepartureTime({u, v}, time);
			double cost_u = l.CostAt(time_u);
			if (time_u == INFTY || cost_u == INFTY) continue;
			if (epsilon_smaller(cost_u, cost - (time - time_u) + penalties[v]))
				fail("Smaller cost when reconstructing exact path.");
			if (epsilon_bigger(cost_u, cost - (time - time_u) + penalties[v])) continue;
			// Add u to path.
			path.push_back(u);
			S = Su;
			v = u;
			cost = cost_u;
			time = time_u;
			found_next = true;
			break;
		}
		if (!found_next) fail("Unable to reconstruct exact path.");
	}

	// Check that the expected duration of the path matches what we reconstructed.
	if (epsilon_different(vrp.BestDurationRoute(reverse(path)).duration, route_duration)) fail("Path duration is not as expected.");
	return reverse(path);
}

void add_to_queue(vector<vector<vector<unordered_set<VertexSet>>>>& q, int k, Vertex v, VertexSet S, LabelSequenceTD& l,
				  int p, int base)
{
	double next_p = l.NextBound(p+base);
	if (next_p == INFTY) return;
	q[next_p-base][k][v].insert(S);
}

MLBStatus run_exact(const VRPInstance& vrp, const NGLInfo& ngl_info, const vector<double>& penalties,
					const BoundingTree& B, const Duration &time_limit, double* lb, Route* UB, json *log)
{
	int base = floor(*lb+EPS);
	int gap = floor(UB->duration+EPS) - base;
	MLBStatus status = MLBStatus::Finished;
	MLBExecutionLog mlb_log(true);

	goc::Stopwatch rolex(true), rolex_temp(false);

	// Map L: Core (k, r, v, S) -> LabelSequence
	VertexSet all_vertices;
	for (Vertex v: vrp.D.Vertices()) all_vertices.set(v);
	int n = vrp.D.VertexCount();
	vector<unordered_map<VertexSet, LabelSequenceTD>> L(n);
	vector<vector<vector<unordered_set<VertexSet>>>> q(gap+1, vector<vector<unordered_set<VertexSet>>>(n+1, vector<unordered_set<VertexSet>>(n)));
	L[vrp.o][goc::create_bitset<MAX_N>({vrp.o})] = LabelSequenceTD({LabelSequenceTD::Initial(vrp.tw[vrp.o], -penalties[vrp.o], *lb)});
	q[0][1][vrp.o].insert(goc::create_bitset<MAX_N>({vrp.o}));

	TableStream output(&clog, 1.0);
	output.AddColumn("LB", 10);
	output.AddColumn("k", 10);
	output.WriteHeader();
	for (int p = 0; p < gap+1; ++p)
	{
		*lb = max(*lb, (double)base+p);
		for (int k = 1; k < n; ++k)
		{
			bool did_some_processing = false;
			for (Vertex v: vrp.D.Vertices())
			{
				for (auto& S: q[p][k][v])
				{
					did_some_processing = true;
					if (rolex.Peek() > time_limit) { mlb_log.status = MLBStatus::TimeLimitReached; break; } // Check time limit.
					mlb_log.processed_count++;

					auto& l = L[v][S]; // Get label sequence with Core (v, S).
					rolex_temp.Reset().Resume();
					B.Bound(S, v, l, base+p, UB->duration); // Apply bounds to the labels in l.
					*mlb_log.bounding_time += rolex_temp.Pause();

					add_to_queue(q, k, v, S, l, p + 1, base);
					auto l_p = l.WithCompletionBound(base + p); // Get labels with bound in [base+p, base+p+1).
					if (l_p.Empty()) continue;
					mlb_log.extended_count += l_p.Count();

					// Extend labels in l_p.
					for (Vertex w: vrp.D.Vertices())
					{
						// Check that extension to w is feasible.
						if (!vrp.D.IncludesArc({v, w})) continue;
						if (S.test(w)) continue;

						// Get latest departere time from w to reach the remaining vertices.
						double LDT_w = INFTY;
						for (Vertex u: vrp.D.Vertices()) if (u != w && !S.test(u)) LDT_w = min(LDT_w, vrp.LDT[w][u]);
						double LDTw_at_v = vrp.DepartureTime({v,w}, LDT_w);
						if (LDTw_at_v == INFTY) continue; // If it is infeasible to reach w before LDT_w, then stop.

						rolex_temp.Reset().Resume();
						auto l_w = l_p.Extend(vrp, ngl_info, S.count(), v, w, penalties[w], LDTw_at_v); // r is not used in extension.
						*mlb_log.extension_time += rolex_temp.Pause();
						if (l_w.Empty()) continue; // If no extension is feasible, continue.

						// Add extension to queue and L.
						auto S_w = S;
						S_w.set(w);
						rolex_temp.Reset().Resume();
						bool exist_new_non_dominated = L[w][S_w].DominateBy(l_w, true);
						*mlb_log.domination_time += rolex_temp.Pause();
						if (exist_new_non_dominated) q[p][k+1][w].insert(S_w);
					}
				}
			}
			if (did_some_processing) output.WriteRow({STR(base+p), STR(k)});
		}

		// Check if a full route was enumerated.
		if (includes_key(L[vrp.d], all_vertices))
		{
			clog << "Found exact solution." << endl;
			break;
		}
	}

	auto& solution = L[vrp.d][all_vertices];

	if (solution.Empty())
	{
		clog << "The problem is infeasible" << endl;
		*lb = INFTY;
		UB->duration = INFTY;
	}
	else
	{
		GraphPath opt_path = reconstruct_path(vrp, penalties, L);
		*UB = vrp.BestDurationRoute(opt_path);
		*lb = UB->duration;
	}
	mlb_log.time = rolex.Peek();
	*log = mlb_log;

	return status;
}
}