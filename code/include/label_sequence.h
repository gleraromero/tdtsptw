//
// Created by Gonzalo Lera Romero on 28/04/2020.
//

#ifndef TDTSPTW_LABEL_SEQUENCE_H
#define TDTSPTW_LABEL_SEQUENCE_H

#include <vector>
#include "core.h"
#include "ngl_info.h"
#include "vrp_instance.h"
#include "goc/goc.h"

namespace tdtsptw
{
// Interface that represents a sequence of ordered disjoint labels which are not dominated that share the same Core.
class LabelSequence
{
public:
	virtual void DominateBy(LabelSequence& l2) = 0;

	virtual bool Empty() const = 0;

	virtual void Extend(const LabelSequence& l, const VRPInstance& vrp, const NGLInfo& ngl_info, const Core& c, goc::Vertex w, Core* c_w, LabelSequence* l_w) const = 0;

	virtual void Fusion(LabelSequence& l) = 0;

	virtual double Merge(LabelSequence& l) const = 0;

	virtual goc::GraphPath MergePath(const LabelSequence& l) const = 0;

	double best_cost; // Best cost of the sequence.
};
} // namespace tdtsptw

#endif //TDTSPTW_LABEL_SEQUENCE_H
