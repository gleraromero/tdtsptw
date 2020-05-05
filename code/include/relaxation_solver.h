//
// Created by Gonzalo Lera Romero on 27/04/2020.
//

#ifndef TDTSPTW_RELAXATION_SOLVER_H
#define TDTSPTW_RELAXATION_SOLVER_H

#include <vector>
#include <unordered_map>
#include <tuple>

#include "goc/goc.h"

#include "vrp_instance.h"
#include "ngl_info.h"
#include "core.h"
#include "bounding_tree.h"

namespace tdtsptw
{
// Reconstructs a path from a .
// 	vrp: instance to work with.
//	ngl_info: information about the NG sets and L path.
//	L: map (k, r, v, S) -> LabelSequence of non-dominated labels.
//	penalties: penalties of the vertices.
// 	c: Core of the label to reconstruct.
//	cost: Cost of the label to reconstruct.
//	time: Arrival time of the label to reconstruct.
//
template<class LS>
goc::GraphPath reconstruct_path(const VRPInstance& vrp, const NGLInfo& ngl_info,
								const std::vector<std::vector<std::vector<std::unordered_map<VertexSet, LS>>>>& L,
								const std::vector<double>& penalties, const Core& c, double cost, double time)
{
	// Variables that contain information about the path to be built.
	goc::Vertex v = c.v;
	int r = ngl_info.L[c.r] == v ? c.r - 1 : c.r;
	VertexSet S = c.S;

	// Get next k-1 vertices.
	goc::GraphPath path = {v};
	for (int k = c.k-1; k > 0; --k)
	{
		bool found_next = false;
		for (goc::Vertex u: vrp.D.Vertices())
		{
			if (u == v) continue;
			if (!goc::is_subset(S, ngl_info.ExtendNG(ngl_info.N[u], v))) continue; // Check that after visiting (u, v) we can have S as the NGset.
			if (!ngl_info.V[r].test(u) && ngl_info.L[r] != u) continue;

			// Check if some sequence with core (k, r, v, _) has a label with  cost=_cost_.
			for (auto& s_l: L[k][r][u])
			{
				auto& s = s_l.first;
				auto& l = s_l.second;
				if (s.test(v)) continue; // If v \in s, then adding (u, v) forms an ng-cycle.
				// Verify that (S_u \cap N[v]) \cup {v} = S_v (i.e. extension of L through v gives the last step).
				if (ngl_info.ExtendNG(s, v) != S) continue;
				double t = vrp.DepartureTime({u, v}, time);
				if (t == goc::INFTY) continue;
				double cost_t = cost - (time - t) + penalties[v];
				double actual_cost = l.CostAt(t);
				if (goc::epsilon_smaller(actual_cost, cost_t)) goc::fail("Smaller cost than expected");
				if (goc::epsilon_equal(actual_cost, cost_t))
				{
					// Move to next label.
					found_next = true;
					path.push_back(u);
					S = s;
					if (ngl_info.L[r] == u) r--;
					v = u;
					time = t;
					cost = cost_t;
					break;
				}
			}

			if (found_next) break; // If next vertex was found, continue.
		}
		if (!found_next)
		{
			goc::fail("Could not reconstruct path");
		}
	}

	return goc::reverse(path);
}

// LS is the class of the LabelSequence used for the relaxation (it can be LabelSequenceTD or LabelSequenceTI).
// Runs a relaxed TSP algorithm with the relaxation specified in LS.
// 	vrp_f: forward instance
// 	vrp_b: backward instance
//	ngl_info_f: NGL structure for forward instance
//	ngl_info_b: NGL structure for backward instance
//	penalties: vertex penalties
//	B: bounding tree to keep non-dominated labels for bounding (nullptr if no bounds need to be kept)
//	time_limit: maximum execution time
//	opt: [output] best route
//	opt_cost: [output] best route cost
//	log: [output] log of the execution.
template<class LS>
goc::BLBStatus run_relaxation(const VRPInstance& vrp_f, const VRPInstance& vrp_b, const NGLInfo& ngl_info_f,
					const NGLInfo& ngl_info_b, const std::vector<double>& penalties, BoundingTree<LS>* B, const goc::Duration& time_limit,
					goc::Route* opt, double* opt_cost, nlohmann::json* log)
{
	// Create execution logs.
	goc::BLBExecutionLog blb_log(true);
	blb_log.status = goc::BLBStatus::Finished;
	blb_log.forward_log = goc::MLBExecutionLog(true);
	blb_log.backward_log = goc::MLBExecutionLog(true);
	goc::MLBExecutionLog* mlb_log[2] = { &(blb_log.forward_log.Value()), &(blb_log.backward_log.Value()) };
	goc::stretch_to_size(*mlb_log[0]->count_by_length, vrp_f.D.VertexCount(), 0);
	goc::stretch_to_size(*mlb_log[1]->count_by_length, vrp_f.D.VertexCount(), 0);

	goc::Stopwatch rolex(true), rolex_temp(false);

	// *** First step: Run forward and backward algorithms until they meet (k_forward + k_backward = n+1).

	// Map L: Core (k, r, v, S) -> LabelSequence
	int n = vrp_f.D.VertexCount();
	std::vector<std::vector<std::vector<std::unordered_map<VertexSet, LS>>>> L[2] = {
			std::vector<std::vector<std::vector<std::unordered_map<VertexSet, LS>>>>(n+1, std::vector<std::vector<std::unordered_map<VertexSet, LS>>>(ngl_info_f.L.size(), std::vector<std::unordered_map<VertexSet, LS>>(n))),
			std::vector<std::vector<std::vector<std::unordered_map<VertexSet, LS>>>>(n+1, std::vector<std::vector<std::unordered_map<VertexSet, LS>>>(ngl_info_f.L.size(), std::vector<std::unordered_map<VertexSet, LS>>(n)))
	}; // Keep one for each direction (0=forward, 1=backward).
	int k_dir[2] = {1, 1}; // k_dir[d] indicates the maximum length of labels enumerated for direction d.
	int c_dir[2] = {1, 1}; // c_dir[d] indicates the number of sequences to extend next in direction d.

	// Add initial routes for forward and backward.
	L[0][1][0][vrp_f.o][goc::create_bitset<MAX_N>({vrp_f.o})] = LS({LS::Initial(vrp_f.tw[vrp_f.o], -penalties[vrp_f.o])});
	L[1][1][0][vrp_b.o][goc::create_bitset<MAX_N>({vrp_b.o})] = LS({LS::Initial(vrp_b.tw[vrp_b.o], -penalties[vrp_b.o])});

	// Keep extending while merging is not possible (merge should yield complete routes with n vertices).
	while (k_dir[0] + k_dir[1] < n+1)
	{
		if (rolex.Peek() > time_limit) { blb_log.status = goc::BLBStatus::TimeLimitReached; break; } // Check time limit.
		int d = c_dir[0] < c_dir[1] ? 0 : 1; // Get the direction which needs to extend the least labels.

		// Extend direction _best_dir_ to one more vertex.
		int k = k_dir[d];
		auto& vrp = d == 0 ? vrp_f : vrp_b;
		auto& ngl_info = d == 0 ? ngl_info_f : ngl_info_b;

		for (int r = 0; r < ngl_info.L.size(); ++r)
		{
			if (rolex.Peek() > time_limit) { blb_log.status = goc::BLBStatus::TimeLimitReached; break; } // Check time limit.
			for (goc::Vertex v: vrp.D.Vertices())
			{
				if (rolex.Peek() > time_limit) { blb_log.status = goc::BLBStatus::TimeLimitReached; break; } // Check time limit.
				for (auto& s_l: L[d][k][r][v])
				{
					if (rolex.Peek() > time_limit) { blb_log.status = goc::BLBStatus::TimeLimitReached; break; } // Check time limit.

					auto& s1 = s_l.first; // ng memory.
					auto& l1 = s_l.second; // label sequence.

					// Dominate L1 by those L2 whose core Core2 <= Core1.
					rolex_temp.Reset().Resume();
					for (auto& s_l2: L[d][k][r][v])
					{
						auto& s2 = s_l2.first;
						auto& l2 = s_l2.second;

						// Check that Core2 <= Core1.
						if (s1 == s2 || !goc::is_subset(s2, s1)) continue;
						l1.DominateBy(l2);
						if (l1.Empty()) break;
					}
					*mlb_log[d]->domination_time += rolex_temp.Pause();
					if (l1.Empty()) continue; // If all pieces were dominated, then skip.

					Core c(k, r, v, s1);
					rolex_temp.Reset().Resume();
					if (B) B->Add(c, l1); // Add label sequence to the bounding tree if existing.
					*mlb_log[d]->bounding_time += rolex_temp.Pause();

					// Extend l1.
					mlb_log[d]->processed_count += l1.Count();
					mlb_log[d]->extended_count++;
					rolex_temp.Reset().Resume();
					for (goc::Vertex w: vrp.D.Vertices())
					{
						// Check that extension to w is feasible.
						if (s1.test(w)) continue;
						if (!ngl_info.V[r].test(w)) continue;
						if (vrp.prec_count[w] > k) continue;
						if (vrp.suc_count[w] > n - k - 1) continue;
						if (w != ngl_info.L[r + 1] && vrp.suc_count[ngl_info.L[r + 1]] > n - k - 2) continue;

						// Extend l1 to w.
						LS l_w = l1.Extend(vrp, ngl_info, c, w, penalties[w]);
						if (l_w.Empty()) continue; // If no extension is feasible, continue.

						// Add extension to queue L, and perform a fusion with the existing sequence with the same
						// core if it exists.
						Core c_w(k+1, r + (ngl_info.L[r+1] == w), w, ngl_info.ExtendNG(s1, w));
						L[d][c_w.k][c_w.r][c_w.v].insert({c_w.S, LS()}).first->second.DominateBy(l_w, true);
						mlb_log[d]->enumerated_count++;
						(*mlb_log[d]->count_by_length)[c_w.k]++;
					}
					*mlb_log[d]->extension_time += rolex_temp.Pause();
				}
			}
		}

		// Update c_dir and k_dir.
		k_dir[d]++;
		c_dir[d] = 0;
		for (int r = 0; r < ngl_info.L.size(); ++r)
			for (goc::Vertex v: vrp.D.Vertices())
				c_dir[d] += L[d][k_dir[d]][r][v].size();
	}

	// *** Second step: Merge forward and backward labels to find the best solution.
	rolex_temp.Reset().Resume();
	int k_f = k_dir[0];
	int k_b = k_dir[1];

	double best_cost = goc::INFTY;
	Core best_forward, best_backward; // forward and backward labels that merged give the best_cost.
	double best_cost_forward, best_cost_backward; // Cost when the forward and backward label merge.
	double best_time; // Time when the forward and backward merge (expressed in terms of forward instance).
	for (int r = 0; r < ngl_info_f.L.size(); ++r)
	{
		for (goc::Vertex v: vrp_f.D.Vertices())
		{
			if ((!ngl_info_f.V[r].test(v) && ngl_info_f.L[r] != v) || (r < ngl_info_f.L.size() - 1 && ngl_info_f.L[r+1] == v)) continue;
			if (rolex.Peek() > time_limit) { blb_log.status = goc::BLBStatus::TimeLimitReached; break; } // Check time limit.

			// If we are standing on a vertex from L then both labels must have the r_value in that vertex, otherwise, at one of distance.
			int r_b = ngl_info_f.L.size() - r - (ngl_info_f.L[r] != v) - 1; // r_b is the r_value for a mergeable backward label.
			for (auto &s_l_f: L[0][k_f][r][v])
			{
				for (auto& s_l_b: L[1][k_b][r_b][v])
				{
					auto &s_f = s_l_f.first, &s_b = s_l_b.first;
					auto &l_f = s_l_f.second, &l_b = s_l_b.second;

					// Check if labels can be merged.
					if ((s_f & s_b) != goc::create_bitset<MAX_N>({v})) continue; // Sf \cap Sb = {v}.
					double merge_time, cost_forward, cost_backward;
					double merged_cost = l_f.Merge(l_b, -penalties[v], &merge_time, &cost_forward, &cost_backward);
					if (merged_cost < best_cost)
					{
						best_cost = merged_cost;
						best_time = merge_time;
						best_forward = Core(k_f, r, v, s_f);
						best_backward = Core(k_b, r_b, v, s_b);
						best_cost_forward = cost_forward;
						best_cost_backward = cost_backward;
					}
				}
			}
		}
	}
	*blb_log.merge_time += rolex_temp.Pause();

	if (best_cost == goc::INFTY)
		goc:: fail("Relaxation error: Should always find a best route if problem is feasible.");

	*opt_cost = best_cost;

	// Reconstruct path by appending the forward and backward paths.
	goc::GraphPath best_route_path = reconstruct_path(vrp_f, ngl_info_f, L[0], penalties, best_forward, best_cost_forward, best_time);
	goc::GraphPath best_backward_path = reconstruct_path(vrp_b, ngl_info_b, L[1], penalties, best_backward, best_cost_backward, -best_time);
	for (int i = (int)best_backward_path.size()-2; i >= 0; --i) best_route_path.push_back(best_backward_path[i]);

	double best_route_penalty = goc::sum<goc::Vertex>(best_route_path, [&] (goc::Vertex v) { return penalties[v]; });
	*opt = goc::Route(best_route_path, 0.0, best_cost + best_route_penalty);

	// Log total execution time.
	*blb_log.time = rolex.Peek();
	*log = blb_log;

	return blb_log.status;
}
} // namespace tdtsptw

#endif //TDTSPTW_RELAXATION_SOLVER_H
