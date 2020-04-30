//
// Created by Gonzalo Lera Romero on 28/04/2020.
//

#include "label_sequence_ti.h"

using namespace std;
using namespace goc;

namespace tdtsptw
{
namespace
{
// Dominates prefix from l2 by l1 including waiting times.
void dominate_label(const LabelSequenceTI::Label& l1, LabelSequenceTI::Label& l2)
{
	if (l2.early == INFTY) return; // All label is already dominated.
	if (epsilon_bigger(l1.cost, l2.cost)) return;
	if (epsilon_bigger(l1.early, l2.early)) return; // No prefix can be dominated.
	double dominance_end = l1.late + (l2.cost - l1.cost); // Dominate all overlapped domain and waiting time until cost is the same.
	if (epsilon_bigger_equal(dominance_end, l2.late))
		l2.early = INFTY; // Move early to INFTY to signal complete dominance.
	else
		l2.early = max(dominance_end, l2.early); // Otherwise only dominate prefix.
}
}

LabelSequenceTI::LabelSequenceTI(const vector<Label>& labels)
		: sequence(labels)
{

}

void LabelSequenceTI::DominateBy(const LabelSequenceTI& L2, bool include_dominating_labels)
{
	// If this sequence has no labels, then there is nothing to dominate.
	if (sequence.empty())
	{
		if (include_dominating_labels) sequence = L2.sequence;
		return;
	}
	// If no dominating labels exist, then do nothing.
	if (L2.sequence.empty()) return;

	auto& s1 = sequence;
	auto s2 = L2.sequence; // Copy L2 sequence in order to modify it.

	// Add fictitious labels at the end of s1 and s2 for simplicity.
	s1.emplace_back(Label(nullptr, 0, INFTY, INFTY, INFTY));
	s2.emplace_back(Label(nullptr, 0, INFTY, INFTY, INFTY));

	vector<Label> result_seq;
	int i = 0, j = 0;
	double t = -INFTY; // latest point covered.

	// Merge labels until both have reached the fictitious label which must be last because of INFTY domain.
	while (i != s1.size() - 1 || j != s2.size() - 1)
	{
		if (j > 0) dominate_label(s2[j-1], s1[i]); // Dominate by previous label because of waiting times.
		if (i > 0) dominate_label(s1[i-1], s2[j]); // Dominate by previous label because of waiting times.

		// Move i and j to labels that have some domain after the last point covered (t).
		if (epsilon_smaller_equal(s1[i].late, t) || epsilon_bigger(s1[i].early, s1[i].late)) { ++i; continue; }
		if (epsilon_smaller_equal(s2[j].late, t) || epsilon_bigger(s2[j].early, s2[j].late)) { ++j; continue; }

		// Move early times of labels beyond t.
		s1[i].early = max(s1[i].early, t);
		s2[j].early = max(s2[j].early, t);

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

		// Add the next piece.
		if (winner_label == 1 || (winner_label == 2 && include_dominating_labels))
		{
			// If the label to be added is another part of the same label added before, then extend it.
			if (!result_seq.empty() && epsilon_equal(result_seq.back().late, winner.early) && result_seq.back().prev == winner.prev)
				result_seq.back().late = t;
				// Otherwise, add the winner part.
			else
				result_seq.emplace_back(Label(winner.prev, winner.v, winner.cost, winner.early, t));
		}
	}
	sequence = result_seq;
	if (!Validate())
	{
		clog << sequence << endl;
		clog << L2.sequence << endl;
		clog << include_dominating_labels << endl;
		fail("Wrong");
	}
}

bool LabelSequenceTI::Empty() const
{
	return sequence.empty();
}

int LabelSequenceTI::Count() const
{
	return sequence.size();
}

LabelSequenceTI LabelSequenceTI::Extend(const VRPInstance& vrp, const NGLInfo& ngl_info, const Core& c, goc::Vertex w, double penalty_w) const
{
	LabelSequenceTI Lw;
	if (this->sequence.empty()) fail("Sequence must not be empty.");

	Vertex v = c.v;
	auto& tau_vw = vrp.tau[v][w]; // Precondition: tau is preprocessed such that it always arrives in [a_w, b_w].

	// If it is not feasible to reach w before its deadline, return empty.
	if (epsilon_bigger(this->sequence.front().early, max(dom(tau_vw)))) return Lw;

	int j = 0;
	for (auto& l: this->sequence)
	{
		if (epsilon_bigger(l.early, max(dom(tau_vw)))) break;
		double min_tau_lw = INFTY; // minimum travel time departing on l.
		double early_lw = INFTY; // earliest arrival to w departing on l.
		while (j < tau_vw.PieceCount())
		{
			if (epsilon_smaller(max(dom(tau_vw[j])), l.early)) { ++j; continue; }
			if (epsilon_bigger(min(dom(tau_vw[j])), l.late)) break;

			// Earliest arrival to w must be departing at l.early (FIFO).
			if (tau_vw[j].domain.Includes(l.early)) early_lw = l.early + tau_vw[j](l.early);

			// Calculate overlap between l and tau_vw[j].
			double overlap_left = max(min(dom(tau_vw[j])), l.early);
			double overlap_right = min(max(dom(tau_vw[j])), l.late);

			// Calculate minimum travel time to w in overlap zone.
			min_tau_lw = min(min_tau_lw, tau_vw[j](overlap_left));
			min_tau_lw = min(min_tau_lw, tau_vw[j](overlap_right));

			if (epsilon_bigger_equal(max(dom(tau_vw[j])), l.late)) break;
			++j;
		}

		double late_lw = min(l.late, max(dom(tau_vw))) + min_tau_lw; // latest arrival to w from l.
		double cost_lw = l.cost + min_tau_lw - penalty_w; // cost of label lw.

		Label lw = Label(&l, w, cost_lw, early_lw, late_lw);

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
				Lw.sequence.back().prev == lw.prev &&
				epsilon_equal(Lw.sequence.back().cost, lw.cost) &&
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
	}
	if (!Lw.Validate())
	{
		fail("Wrong");
	}
	return Lw;
}

LabelSequenceTI::MergedLabel LabelSequenceTI::Merge(LabelSequenceTI& L, double redundant_cost) const
{
	MergedLabel best(nullptr, nullptr, INFTY);
	if (sequence.empty() || L.sequence.empty()) return best;

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

			// Update cost with merge of l and s[j].
			double waiting_time = max(0.0, early_j - l.late); // necessary waiting time between l to depart at s[j].
			double merge_cost = l.cost + s2[j].cost - redundant_cost + waiting_time;
			if (merge_cost < best.cost)
			{
				best.cost = merge_cost;
				best.forward = &l;
				best.backward = &s2[j];
			}

			if (epsilon_bigger(late_j, l.late)) break; // label s2[j] ends after the end of l, then move to next l.
			--j;
		}
		j = min((int)s2.size()-1, j+1); // We go back one label to make sure we consider them all.
	}

	return best;
}

bool LabelSequenceTI::Validate() const
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

		if (sequence[i].prev == sequence[i+1].prev && epsilon_equal(sequence[i].late, sequence[i+1].early))
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

LabelSequenceTI::Label LabelSequenceTI::Initial(Vertex origin, const goc::Interval& time_window, double initial_cost)
{
	return Label(nullptr, origin, initial_cost, time_window.left, time_window.right);
}

LabelSequenceTI::Label::Label(const Label* prev, goc::Vertex v, double cost, double early, double late)
		: prev(prev), v(v), cost(cost), early(early), late(late)
{

}

void LabelSequenceTI::Label::Print(std::ostream& os) const
{
	os << "{ cost: " << cost << ", domain: [" << early << ", " << late << "]";
}

LabelSequenceTI::MergedLabel::MergedLabel()
		: forward(nullptr), backward(nullptr), cost(INFTY)
{

}

LabelSequenceTI::MergedLabel::MergedLabel(const Label* forward, const Label* backward, double cost)
		: forward(forward), backward(backward), cost(cost)
{

}

GraphPath LabelSequenceTI::MergedLabel::Path() const
{
	GraphPath p;
	// Add forward path.
	for (const Label* l = forward; l != nullptr; l = l->prev) p.push_back(l->v);
	p = reverse(p);
	// Add backward path.
	for (const Label* l = backward->prev; l != nullptr; l = l->prev) p.push_back(l->v);
	return p;
}
} // namespace tdtsptw