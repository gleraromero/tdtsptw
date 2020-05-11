//
// Created by Gonzalo Lera Romero on 28/04/2020.
//

#ifndef TDTSPTW_BOUNDING_TREE_H
#define TDTSPTW_BOUNDING_TREE_H

#include "core.h"
#include "label_sequence_td.h"
#include "label_sequence_ti.h"
#include "vrp_instance.h"
#include "ngl_info.h"

namespace tdtsptw
{
struct BoundingTree
{
public:
	BoundingTree(const VRPInstance* vrp, const NGLInfo* ngl_info, const std::vector<double>& penalties);

	// Disables bounding tree, and then all calls to Bound will return lb.
	void Disable();

	void Add(const Core& c, const LabelSequenceTI& L);

	void Add(const Core& c, const LabelSequenceTD& L);

	// Sets the bounds for the labels in L and removes the labels with bound bigger than ub.
	void Bound(const VertexSet& S, goc::Vertex v, LabelSequenceTD& L, double lb, double ub) const;

private:
	bool enabled; // indicates if the bounding is working.
	const VRPInstance* vrp;
	const NGLInfo* ngl_info;
	std::vector<double> penalties;
	double penalty_sum;
	std::vector<std::vector<std::vector<std::vector<std::pair<VertexSet, LabelSequenceTD>>>>> T; // t[k][r][v] = {S, L}.
};
} // namespace tdtsptw

#endif //TDTSPTW_BOUNDING_TREE_H
