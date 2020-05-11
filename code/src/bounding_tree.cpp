//
// Created by Gonzalo Lera Romero on 11/05/2020.
//

#include "bounding_tree.h"

using namespace std;
using namespace goc;
using namespace nlohmann;

namespace tdtsptw
{
BoundingTree::BoundingTree(const VRPInstance* vrp, const NGLInfo* ngl_info, const vector<double>& penalties)
	: vrp(vrp), ngl_info(ngl_info), penalties(penalties)
{
	int n = vrp->D.VertexCount();
	T = vector<vector<vector<vector<pair<VertexSet, LabelSequenceTD>>>>>(
			n+1, vector<vector<vector<pair<VertexSet, LabelSequenceTD>>>>(
					ngl_info->L.size()+1, vector<vector<pair<VertexSet, LabelSequenceTD>>>(n)));

	penalty_sum = sum(penalties);
	enabled = true;
}

void BoundingTree::Disable()
{
	enabled = false;
}

void BoundingTree::Add(const Core& c, const LabelSequenceTI& L)
{
	fail("Bounding tree does not support LabelSequenceTI.");
}

void BoundingTree::Add(const Core& c, const LabelSequenceTD& L)
{
	// Compute the Core of labels that can be bounded by L.
	int n = vrp->D.VertexCount();
	int k = n - c.k + 1;
	int r = ngl_info->L.size() - c.r - 1;
	if (ngl_info->L[r] != c.v) r--;
	T[k][r][c.v].push_back({c.S, L});
}

// Sets the bounds for the labels in L and removes the labels with bound bigger than ub.
void BoundingTree::Bound(const VertexSet& S, Vertex v, LabelSequenceTD& L, double lb, double ub) const
{
	// If the bounding tree is disabled always return lb as completion bound.
	if (!enabled)
	{
		for (auto& p: L.sequence) p.completion_bound = lb;
		return;
	}

	if (L.sequence.empty()) return;
	int k = S.count();
	int r = 0;
	while (r+1 < ngl_info->L.size() && S.test(ngl_info->L[r+1])) ++r;
	for (auto& S_L: T[k][r][v])
	{
		auto& S2 = S_L.first;
		if ((S & S2) != create_bitset<MAX_N>({v})) continue; // Check that L2 is mergeable with L.
		auto& L2 = S_L.second;

		auto& s1 = L.sequence;
		auto& s2 = L2.sequence;

		// s2 comes from an opposite direction algorithm, therefore its labels are sorted backwards.
		int j = s2.size()-1;
		for (auto& l: s1)
		{
			if (l.completion_bound == -INFTY) l.completion_bound = INFTY;

			// Try to merge l with all opposite labels that come afterwards with overlap domain and the first
			// after that (for waiting purposes).
			while (j >= 0)
			{
				double late_j = -s2[j].early, early_j = -s2[j].late;
				if (epsilon_smaller(late_j, l.early)) { --j; continue; } // Move j until s2[j] is ahead of l.

				double merge_cost = INFTY;
				// If there is no overlap, then compute cost of waiting.
				if (epsilon_bigger(early_j, l.late))
				{
					double waiting_time = early_j - l.late; // necessary waiting time between l to depart at s[j].
					merge_cost = min(merge_cost, l.CostAt(l.late) + s2[j].CostAt(-early_j) + waiting_time + penalties[v]);
				}
				// If there is overlap, compute the minimum which must be in an extreme.
				else
				{
					double early_overlap = max(l.early, early_j);
					merge_cost = min(merge_cost, l.CostAt(early_overlap) + s2[j].CostAt(-early_overlap) + penalties[v]);

					double late_overlap = min(l.late, late_j);
					merge_cost = min(merge_cost, l.CostAt(late_overlap) + s2[j].CostAt(-late_overlap) + penalties[v]);
				}
				// Update completion bound.
				l.completion_bound = min(l.completion_bound, merge_cost + penalty_sum);

				if (epsilon_bigger(late_j, l.late)) break; // label s2[j] ends after the end of l, then move to next l.
				--j;
			}
			j = min((int)s2.size()-1, j+1); // We go back one label to make sure we consider them all.
		}
	}
	// Remove labels with cb >= ub.
	L.sequence.erase(remove_if(L.sequence.begin(), L.sequence.end(),
			[&] (const LabelSequenceTD::Label& l) { return epsilon_bigger(l.completion_bound, ub); }),
					L.sequence.end());
}
} // namespace tdtsptw