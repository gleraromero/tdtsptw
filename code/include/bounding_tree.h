//
// Created by Gonzalo Lera Romero on 28/04/2020.
//

#ifndef TDTSPTW_BOUNDING_TREE_H
#define TDTSPTW_BOUNDING_TREE_H

#include "core.h"
#include "label_sequence.h"

namespace tdtsptw
{
template<typename LS>
struct BoundingTree
{
public:
	void Add(const Core& c, const LS& L)
	{

	}

	// Sets the bounds for the labels in L and removes the labels with bound bigger than ub.
	void Bound(LS& L, double ub) const
	{
		for (auto& l: L.sequence) if (l.completion_bound == -goc::INFTY) l.completion_bound = ub;
	}
};
} // namespace tdtsptw

#endif //TDTSPTW_BOUNDING_TREE_H
