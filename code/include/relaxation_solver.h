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
// LS is the class of the LabelSequence used for the relaxation (it can be LabelSequenceTD or LabelSequenceTI).
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
	L[0][1][0][vrp_f.o][goc::create_bitset<MAX_N>({vrp_f.o})] = LS({typename LS::Label(nullptr, vrp_f.o, 0.0, vrp_f.a[vrp_f.o], vrp_f.b[vrp_f.o])});
	L[1][1][0][vrp_b.o][goc::create_bitset<MAX_N>({vrp_b.o})] = LS({typename LS::Label(nullptr, vrp_b.o, 0.0, vrp_b.a[vrp_b.o], vrp_b.b[vrp_b.o])});

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
						VertexSet S_w = s1 & ngl_info.N[w]; // Update S_w according to NG rules.
						S_w.set(w);
						Core c_w(k+1, r + (ngl_info.L[r+1] == w), w, S_w);
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

	typename LS::MergedLabel best_route;
	for (int r = 0; r < ngl_info_f.L.size(); ++r)
	{
		for (goc::Vertex v: vrp_f.D.Vertices())
		{
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
//					if (l_f.best_cost + l_b.best_cost + penalties[v] >= *opt_cost) continue;
					auto merged_route = l_f.Merge(l_b, -penalties[v]);
					if (merged_route.cost < best_route.cost) best_route = merged_route;
				}
			}
		}
	}
	*blb_log.merge_time += rolex_temp.Pause();

	if (best_route.cost == goc::INFTY)
		goc:: fail("Relaxation error: Should always find a best route if problem is feasible.");

	*opt_cost = best_route.cost;
	// Construct optimum route.
	if (*opt_cost != goc::INFTY) *opt = vrp_f.BestDurationRoute(best_route.Path());

	// Validation.
	if (*opt_cost != goc::INFTY && goc::epsilon_different(*opt_cost, opt->duration - goc::sum<goc::Vertex>(opt->path, [&] (goc::Vertex v) { return penalties[v]; })))
		goc::fail("Relaxation error: Costs do not match " + STR(*opt_cost) + " " + STR(opt->duration - goc::sum<goc::Vertex>(opt->path, [&] (goc::Vertex v) { return penalties[v]; })));

	// Log total execution time.
	*blb_log.time = rolex.Peek();
	*log = blb_log;

	return blb_log.status;
}
} // namespace tdtsptw

#endif //TDTSPTW_RELAXATION_SOLVER_H
