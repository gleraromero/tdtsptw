//
// Created by Gonzalo Lera Romero on 28/04/2020.
//

#ifndef TDTSPTW_LABEL_SEQUENCE_TI_H
#define TDTSPTW_LABEL_SEQUENCE_TI_H

#include <vector>
#include "core.h"
#include "ngl_info.h"
#include "vrp_instance.h"
#include "goc/goc.h"

namespace tdtsptw
{
// Specialization of class LabelSequence for Time Independent labels (i.e. representing constant functions).
class LabelSequenceTI
{
public:
	struct Label
	{
		const Label* prev; // Previous label.
		goc::Vertex v; // Last vertex.
		double cost, early, late;

		Label(const Label* prev, goc::Vertex v, double cost, double early, double late);
	};

	struct MergedLabel
	{
		const Label* forward;
		const Label* backward;
		double cost;

		MergedLabel();

		MergedLabel(const Label* forward, const Label* backward, double cost);

		goc::GraphPath Path() const;
	};

	LabelSequenceTI() = default;

	LabelSequenceTI(const std::vector<Label>& labels);

	// Dominates all the labels which have a bigger cost than those in L with the same domain.
	//  include_dominating_labels: indicates if the labels from L that dominate should be included in the sequence.
	void DominateBy(const LabelSequenceTI& L, bool include_dominating_labels=false);

	// Returns if the sequence has no labels.
	bool Empty() const;

	// Extends the sequence to another vertex w into a new sequence.
	LabelSequenceTI Extend(const VRPInstance& vrp, const NGLInfo& ngl_info, const Core& c, goc::Vertex w, double penalty_w) const;

	// Merges this sequence with the opposite direction sequence L, and returns the minimum cost of such merge.
	// 	redundant_cost: is the cost that is accounted for in both labels and therefore must be substracted.
	MergedLabel Merge(LabelSequenceTI& L, double redundant_cost) const;

private:
	std::vector<Label> sequence;
};
} // namespace tdtsptw

#endif //TDTSPTW_LABEL_SEQUENCE_TI_H
