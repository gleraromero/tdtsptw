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
class LabelSequenceTD : public goc::Printable
{
public:
	class Label : public goc::Printable
	{
	public:
		double early, late; // Earliest and latest arrival time for the label.
		double slope, intercept; // slope and intercept of the function.
		double completion_bound=-goc::INFTY; // completion bound of the label

		Label(double early, double late, double slope, double intercept);

		// Returns the cost of the function at time t.
		// If time t is before early, an exception is raised.
		// If time t is after late, the cost including waiting times is considered.
		double CostAt(double t) const;

		virtual void Print(std::ostream& os) const;

		static Label FromPoints(double early, double late, double early_cost, double late_cost);
	};

	LabelSequenceTD() = default;

	LabelSequenceTD(const std::vector<Label>& labels);

	// Dominates all the labels which have a bigger cost than those in L with the same domain.
	//  include_dominating_labels: indicates if the labels from L that dominate should be included in the sequence.
	// Returns if any non dominate label from L was kept.
	bool DominateBy(const LabelSequenceTD& L, bool include_dominating_labels=false);

	// Returns if the sequence has no labels.
	bool Empty() const;

	// Returns the number of labels.
	int Count() const;

	// Returns the earliest and latest arrival times.
	goc::Interval Domain() const;

	// Extends the sequence to another vertex w into a new sequence.
	LabelSequenceTD Extend(const VRPInstance& vrp, const NGLInfo& ngl_info, int k, goc::Vertex v, goc::Vertex w,
			double penalty_w, double max_dom=goc::INFTY) const;

	// Merges this sequence with the opposite direction sequence L, and returns the minimum cost of such merge.
	// 	redundant_cost: is the cost that is accounted for in both labels and therefore must be substracted.
	//  merge_t: returns the merge time when the best cost is achieved (with respect to this sequence).
	// 	cost_f: cost of forward label at time merge_t.
	// 	cost_b: cost of backward label at time merge_t.
	double Merge(LabelSequenceTD& L, double redundant_cost, double* merge_t, double* cost_f, double* cost_b) const;

	bool Validate() const;

	static Label Initial(const goc::Interval& time_window, double initial_cost, double initial_bound=-goc::INFTY);

	// Returns the cost of arriving at time t, including waiting times.
	double CostAt(double t) const;

	// Returns a new label sequence with the labels l with floor(bound(l)) == bound.
	LabelSequenceTD WithCompletionBound(double bound) const;

	// Returns the next floor(bound(l)) of a label with floor(bound(l)) >= bound.
	double NextBound(double bound) const;

	// Returns the arrival time t with lower cost.
	void LowestCostArrival(double* t, double* cost) const;

	void Print(std::ostream& os) const;

	std::vector<Label> sequence;
};
} // namespace tdtsptw

#endif //TDTSPTW_LABEL_SEQUENCE_TD_H
