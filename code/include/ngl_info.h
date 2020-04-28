//
// Created by Gonzalo Lera Romero on 24/04/2020.
//

#ifndef TDTSPTW_NGL_INFO_H
#define TDTSPTW_NGL_INFO_H

#include <vector>
#include "goc/goc.h"
#include "vrp_instance.h"

namespace tdtsptw
{
struct NGLInfo
{
	int delta; // maximum size of each neighbour.
	std::vector<VertexSet> N; // neighbours.
	goc::GraphPath L; // path for the ngL relaxation.
	std::vector<VertexSet> V; // V[r] are the vertices that may be visited after L[r] and before L[r+1].
};

// Creates NGLInfo structure with neighbours containing the _delta_ closest vertices.
// The L path consists of the longest path in the precedence digraph.
void create_default_nginfo(const VRPInstance& vrp, int delta, NGLInfo* forward, NGLInfo* backward);
} // namespace tdtsptw

#endif //TDTSPTW_NGL_INFO_H
