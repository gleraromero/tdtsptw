//
// Created by Gonzalo Lera Romero on 28/04/2020.
//

#include "label_sequence_td.h"

using namespace std;
using namespace goc;

namespace tdtsptw
{
namespace
{
// Dominates prefix from l2 by l1 including waiting times.
void dominate_label(const LabelSequenceTD::Label& l1, LabelSequenceTD::Label& l2)
{
	if (l2.early == INFTY) return; // All label is already dominated.
	if (epsilon_bigger(l1.early, l2.early)) return; // No prefix can be dominated.
	if (epsilon_bigger(l1.CostAt(l2.early), l2.CostAt(l2.early))) return; // The cost of l1 at the start of l2 is bigger, then no prefix is dominated.
	// If l2 is fully dominated, move early to INFTY.
	if (epsilon_smaller_equal(l1.CostAt(l2.late), l2.CostAt(l2.late)))
	{
		l2.early = INFTY;
		return;
	}
	// Then, l2 is partially dominated, check intersection.
	double intersection = epsilon_smaller_equal(l1.slope, l2.slope) ? INFTY : (l1.intercept - l2.intercept) / (l2.slope - l1.slope);

	// If the intersection happens after late(l1), then compute the real intersection which is between
	// the waiting time from l1 after late and l2.
	if (epsilon_bigger(intersection, l1.late))
	{
		double slope_waiting = 1.0, intercept_waiting = l1.CostAt(l1.late) - l1.late;
		intersection = epsilon_equal(l2.slope, slope_waiting) ? INFTY : (intercept_waiting - l2.intercept) / (l2.slope - slope_waiting);
	}

	// Then, l2 is dominated until the intersection.
	if (epsilon_bigger_equal(intersection, l2.late))
		l2.early = INFTY; // Move early to INFTY to signal complete dominance.
	else
		l2.early = max(intersection, l2.early); // Otherwise only dominate prefix.
}
}

LabelSequenceTD::LabelSequenceTD(const vector<Label>& labels)
		: sequence(labels)
{

}

bool LabelSequenceTD::DominateBy(const LabelSequenceTD& L2, bool include_dominating_labels)
{
	// If no dominating labels exist, then do nothing.
	if (L2.sequence.empty()) return false;

	// If this sequence has no labels, then there is nothing to dominate.
	if (sequence.empty())
	{
		if (include_dominating_labels) sequence = L2.sequence;
		return true;
	}

	auto& s1 = sequence;
	auto s2 = L2.sequence; // Copy L2 sequence in order to modify it.

	// Add fictitious labels at the end of s1 and s2 for simplicity.
	s1.emplace_back(Label::FromPoints(INFTY, INFTY, INFTY, INFTY));
	s2.emplace_back(Label::FromPoints(INFTY, INFTY, INFTY, INFTY));

	vector<Label> result_seq;
	int i = 0, j = 0;
	double t = -INFTY; // latest point covered.
	Label last_consolidated = Label::FromPoints(INFTY, INFTY, INFTY, INFTY); // Last label that was verified as not dominated.
	bool included_labels_from_L2 = false;

	// Merge labels until both have reached the fictitious label which must be last because of INFTY domain.
	while (i != s1.size() - 1 || j != s2.size() - 1)
	{
		dominate_label(last_consolidated, s1[i]);
		dominate_label(last_consolidated, s2[j]);

		// Move i and j to labels that have some domain after the last point covered (t).
		if (s1[i].early == INFTY && i + 1 < s1.size()) { ++i; continue; }
		if (s2[j].early == INFTY && j + 1 < s2.size()) { ++j; continue; }

		// Move early times of labels beyond t.
		s2[j].early = max(s2[j].early, t);
		s1[i].early = max(s1[i].early, t);

		// Dominate labels between them.
		dominate_label(s2[j], s1[i]);
		dominate_label(s1[i], s2[j]);

		// Compute label that adds that part.
		int winner_label = 0;
		if (epsilon_smaller(s1[i].early, s2[j].early)) winner_label = 1;
		else if (epsilon_smaller(s2[j].early, s1[i].early)) winner_label = 2;
		else winner_label = (epsilon_equal(s1[i].early, s1[i].late) ? 1 : 2); // If they start at the same time it is because one of them is a point that dominates the other one.

		auto& winner = winner_label == 1 ? s1[i] : s2[j];
		auto& loser = winner_label == 1 ? s2[j] : s1[i];
		t = min(min(winner.late, loser.late), loser.early); // late time of the part we must add.

		if (winner_label == 2) included_labels_from_L2 = true;

		// Add the next piece.
		if (winner_label == 1 || (winner_label == 2 && include_dominating_labels))
		{
			// If the label to be added is another part of the same label added before, then extend it.
			if (!result_seq.empty() && epsilon_equal(result_seq.back().late, winner.early) &&
				epsilon_equal(result_seq.back().slope, winner.slope) &&
				epsilon_equal(result_seq.back().intercept, winner.intercept))
				result_seq.back().late = t;
			// Otherwise, add the winner part.
			else
			{
				Label add(winner.early, t, winner.slope, winner.intercept);
				if (!result_seq.empty()) dominate_label(result_seq.back(), add);
				if (add.early < INFTY) result_seq.emplace_back(add);
			}
		}

		// Keep the latest consolidated label.
		last_consolidated = Label(winner.early, t, winner.slope, winner.intercept);
	}
	sequence = result_seq;
	if (!Validate())
	{
		fail("Wrong domination");
	}
	return included_labels_from_L2;
}

bool LabelSequenceTD::Empty() const
{
	return sequence.empty();
}

int LabelSequenceTD::Count() const
{
	return sequence.size();
}

Interval LabelSequenceTD::Domain() const
{
	if (sequence.empty()) return Interval();
	return {sequence.front().early, sequence.back().late};
}

LabelSequenceTD LabelSequenceTD::Extend(const VRPInstance& vrp, const NGLInfo& ngl_info, const Core& c, goc::Vertex w, double penalty_w, double max_dom) const
{
	LabelSequenceTD Lw;
	if (this->sequence.empty()) fail("Sequence must not be empty.");

	Vertex v = c.v;
	auto& tau_vw = vrp.tau[v][w]; // Precondition: tau is preprocessed such that it always arrives in [a_w, b_w].

	// If it is not feasible to reach w before its deadline, return empty.
	if (epsilon_bigger(this->sequence.front().early, max(dom(tau_vw)))) return Lw;

	// Can depart to arrive at b_kw or before.
	double a_kw = vrp.TWP[c.k+1][w].left;
	double b_kw = vrp.TWP[c.k+1][w].right;
	double latest_departure = min(max_dom, vrp.DepartureTime({v, w}, b_kw));
	if (latest_departure == INFTY) return Lw; // If impossible, return empty sequence.

	int j = 0;
	for (auto& l: this->sequence)
	{
		if (epsilon_bigger(l.early, latest_departure)) break;
		while (j < tau_vw.PieceCount())
		{
			if (epsilon_smaller(max(dom(tau_vw[j])), l.early)) { ++j; continue; }
			if (epsilon_bigger(min(dom(tau_vw[j])), l.late)) break;
			if (epsilon_bigger(min(dom(tau_vw[j])), latest_departure)) break;

			// Calculate overlap between l and tau_vw[j].
			double early_overlap = max(min(dom(tau_vw[j])), l.early);
			double late_overlap = min(max(dom(tau_vw[j])), l.late);
			late_overlap = min(late_overlap, latest_departure);
			late_overlap = max(early_overlap, late_overlap);

			// Calculate early and late arriving times.
			double early_lw = max(a_kw, early_overlap + tau_vw[j](early_overlap));
			double late_lw = max(a_kw, late_overlap + tau_vw[j](late_overlap));

			// Calculate cost at borders.
			double early_lw_cost = l.CostAt(early_overlap) + (early_lw - early_overlap) - penalty_w;
			double late_lw_cost = l.CostAt(late_overlap) + (late_lw - late_overlap) - penalty_w;

			// If the earliest and latest arrival times are the same, then we must consider only departing the latest.
			if (epsilon_equal(early_lw, late_lw)) early_lw_cost = late_lw_cost;

			// Label resulting from l + tau_vw[j].
			Label lw = Label::FromPoints(early_lw, late_lw, early_lw_cost, late_lw_cost);

			// Execute domination among new pieces.
			if (!Lw.sequence.empty())
			{
				dominate_label(Lw.sequence.back(), lw);
				dominate_label(lw, Lw.sequence.back());
				if (Lw.sequence.back().early == INFTY) Lw.sequence.pop_back();
			}

			if (lw.early != INFTY)
			{
				// Check if label is an extension of the previous label, if so, extend domain.
				if (!Lw.sequence.empty() &&
					epsilon_equal(Lw.sequence.back().slope, lw.slope) &&
					epsilon_equal(Lw.sequence.back().intercept, lw.intercept) &&
					epsilon_equal(Lw.sequence.back().late, lw.early))
				{
					Lw.sequence.back().late = late_lw;
				}
				else
				{
					// Label could not be squashed, then, add it to the end of the sequence.
					Lw.sequence.emplace_back(lw);
				}
			}

			// Move to next piece of tau_vw if needed.
			if (epsilon_bigger_equal(max(dom(tau_vw[j])), l.late)) break;
			++j;
		}
	}
	if (!Lw.Validate())
	{
		fail("Wrong extension");
	}
	return Lw;
}

double LabelSequenceTD::Merge(LabelSequenceTD& L, double redundant_cost, double* merge_t, double* cost_f, double* cost_b) const
{
	double best_cost = INFTY;
	if (sequence.empty() || L.sequence.empty()) return best_cost;

	auto& s1 = sequence;
	auto& s2 = L.sequence;

	// s2 comes from an opposite direction algorithm, therefore its labels are sorted backwards.
	int j = s2.size()-1;
	for (auto& l: s1)
	{
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
				merge_cost = min(merge_cost, l.CostAt(l.late) + s2[j].CostAt(-early_j) + waiting_time - redundant_cost);
			}
			// If there is overlap, compute the minimum which must be in an extreme.
			else
			{
				double early_overlap = max(l.early, early_j);
				merge_cost = min(merge_cost, l.CostAt(early_overlap) + s2[j].CostAt(-early_overlap) - redundant_cost);

				double late_overlap = min(l.late, late_j);
				merge_cost = min(merge_cost, l.CostAt(late_overlap) + s2[j].CostAt(-late_overlap) - redundant_cost);
			}

			// Update cost with merge of l and s[j].
			if (merge_cost < best_cost)
			{
				best_cost = merge_cost;
				*merge_t = min(l.late, late_j);
				*cost_f = l.CostAt(*merge_t);
				*cost_b = s2[j].CostAt(-*merge_t);
			}

			if (epsilon_bigger(late_j, l.late)) break; // label s2[j] ends after the end of l, then move to next l.
			--j;
		}
		j = min((int)s2.size()-1, j+1); // We go back one label to make sure we consider them all.
	}

	return best_cost;
}

bool LabelSequenceTD::Validate() const
{
	if (sequence.empty()) return true;

	// Check that no piece dominates the following one.
	for (int i = 0; i < sequence.size()-1; ++i)
	{
		auto p1 = sequence[i], p2 = sequence[i+1];
		dominate_label(p1, p2);
		dominate_label(p2, p1);
		if (epsilon_different(p2.early, sequence[i+1].early))
		{
			clog << sequence[i] << " dominates " << sequence[i+1] << endl;
			clog << "Non domination free sequence" << endl;
			return false;
		}
		if (epsilon_different(p1.early, sequence[i].early))
		{
			clog << sequence[i+1] << " dominates " << sequence[i] << endl;
			clog << "Non domination free sequence inverse" << endl;
			return false;
		}
	}

	// Check that is increasing.
	for (int i = 0; i < sequence.size()-1; ++i)
	{
		if (epsilon_bigger(sequence[i].late, sequence[i+1].early))
		{
			clog << sequence << endl;
			clog << i << endl;
			clog << sequence[i] << " " << sequence[i+1] << endl;
			clog << "Broken sequence because of disjointness." << endl;
			return false;
		}

		if (epsilon_equal(sequence[i].late, sequence[i+1].early) &&
			epsilon_equal(sequence[i].slope, sequence[i+1].slope) && epsilon_equal(sequence[i].intercept, sequence[i+1].intercept))
		{
			clog << sequence[i] << " vs " << sequence[i+1] << endl;
			clog << "Sequence not squashed." << endl;
			return false;
		}
	}

	for (auto& l: sequence)
	{
		if (epsilon_bigger(l.early, l.late))
		{
			clog << sequence << endl;
			clog << l << endl;
			clog << "Broken sequence because of empty label." << endl;
			return false;
		}
	}
	return true;
}

LabelSequenceTD::Label LabelSequenceTD::Initial(const goc::Interval& time_window, double initial_cost, double initial_bound)
{
	auto l = Label::FromPoints(time_window.left, time_window.right, initial_cost, initial_cost);
	l.completion_bound = initial_bound;
	return l;
}

double LabelSequenceTD::CostAt(double t) const
{
	if (sequence.empty()) return INFTY;
	if (epsilon_smaller(t, sequence.front().early)) return INFTY;
	double cost = INFTY;
	for (int i = sequence.size()-1; i >= 0; --i)
	{
		if (epsilon_smaller_equal(sequence[i].early, t))
		{
			cost = sequence[i].CostAt(t);
			break;
		}
	}
	return cost;
}

LabelSequenceTD LabelSequenceTD::WithCompletionBound(double bound) const
{
	LabelSequenceTD result;
	for (auto& l: sequence) if (epsilon_equal(floor(l.completion_bound+EPS), bound)) result.sequence.push_back(l);
	return result;
}

double LabelSequenceTD::NextBound(double bound) const
{
	double next = INFTY;
	for (auto& l: sequence) if (floor(l.completion_bound+EPS) >= bound) next = min(next, floor(l.completion_bound+EPS));
	return next;
}

void LabelSequenceTD::LowestCostArrival(double* t, double* cost) const
{
	*t = INFTY;
	*cost = INFTY;
	for (auto& l: sequence)
	{
		if (l.CostAt(l.early) < *cost)
		{
			*t = l.early;
			*cost = l.CostAt(l.early);
		}
		if (l.CostAt(l.late) < *cost)
		{
			*t = l.late;
			*cost = l.CostAt(l.late);
		}
	}
}

void LabelSequenceTD::Print(ostream& os) const
{
	os << sequence;
}

LabelSequenceTD::Label LabelSequenceTD::Label::FromPoints(double early, double late, double early_cost, double late_cost)
{
	double slope, intercept;
	// Case 1: the function is constant.
	if (epsilon_equal(early_cost, late_cost))
	{
		slope = 0;
		intercept = early_cost;
	}
	else
	{
		slope = (late_cost - early_cost) / (late - early);
		intercept = early_cost - slope * early;
	}
	return Label(early, late, slope, intercept);
}

LabelSequenceTD::Label::Label(double early, double late, double slope, double intercept)
	: early(early), late(late), slope(slope), intercept(intercept)
{

}

double LabelSequenceTD::Label::CostAt(double t) const
{
	// Check that t is inside domain or includes waiting times..
	if (epsilon_smaller(t, early)) fail("Cost of label at " + STR(t) + " is infeasible because domain is [" + STR(early) + ", " + STR(late) + "]");
	// If t is after the late time, then we must consider the cost plus the waiting time.
	if (epsilon_bigger_equal(t, late)) return slope * late + intercept  + (t - late);
	// We are inside the domain.
	return slope * t + intercept;
}

void LabelSequenceTD::Label::Print(std::ostream& os) const
{
	os << "{ domain: [" << early << ", " << late << "], image: [" << (slope * early + intercept) << ", " << (slope * late + intercept) << "], bd: " << completion_bound << " }";
}
} // namespace tdtsptw