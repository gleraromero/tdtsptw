//
// Created by Gonzalo Lera Romero on 30/04/2020.
//

#ifndef TDTSPTW_LABEL_SEQUENCE_TD_H
#define TDTSPTW_LABEL_SEQUENCE_TD_H

#include <vector>
#include "core.h"
#include "ngl_info.h"
#include "vrp_instance.h"
#include "goc/goc.h"

namespace tdtsptw
{
// Specialization of class LabelSequence for Time Dependent labels (i.e. representing linear functions).
class LabelSequenceTD
{
public:
	class Label : public goc::Printable
	{
	public:
		const Label* prev; // Previous label.
		goc::Vertex v; // Last vertex.
		double early, late; // Earliest and latest arrival time for the label.
		double slope, intercept; // slope and intercept of the function.

		Label(const Label* prev, goc::Vertex v, double early, double late, double early_cost, double late_cost);

		Label(const Label* prev, double early, double late, goc::Vertex v, double slope, double intercept);

		// Returns the cost of the function at time t.
		// If time t is before early, an exception is raised.
		// If time t is after late, the cost including waiting times is considered.
		double CostAt(double t) const;

		virtual void Print(std::ostream& os) const;
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

	LabelSequenceTD() = default;

	LabelSequenceTD(const std::vector<Label>& labels);

	// Dominates all the labels which have a bigger cost than those in L with the same domain.
	//  include_dominating_labels: indicates if the labels from L that dominate should be included in the sequence.
	void DominateBy(const LabelSequenceTD& L, bool include_dominating_labels=false);

	// Returns if the sequence has no labels.
	bool Empty() const;

	// Returns the number of labels.
	int Count() const;

	// Extends the sequence to another vertex w into a new sequence.
	LabelSequenceTD Extend(const VRPInstance& vrp, const NGLInfo& ngl_info, const Core& c, goc::Vertex w, double penalty_w) const;

	// Merges this sequence with the opposite direction sequence L, and returns the minimum cost of such merge.
	// 	redundant_cost: is the cost that is accounted for in both labels and therefore must be substracted.
	MergedLabel Merge(LabelSequenceTD& L, double redundant_cost) const;

	bool Validate() const;

	static Label Initial(goc::Vertex origin, const goc::Interval& time_window, double initial_cost);

private:
	std::vector<Label> sequence;
};
} // namespace tdtsptw

#endif //TDTSPTW_LABEL_SEQUENCE_TD_H
