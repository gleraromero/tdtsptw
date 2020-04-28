//
// Created by Gonzalo Lera Romero on 28/04/2020.
//

#ifndef TDTSPTW_CORE_H
#define TDTSPTW_CORE_H

#include "goc/goc.h"
#include "vrp_instance.h"

namespace tdtsptw
{
// Represents the main components of a label as presented in the article.
struct Core
{
	int k; // number of visited vertices.
	int r; // index of the last visited vertex from L path (from ngL).
	goc::Vertex v; // last visited vertex.
	VertexSet S; // set of forbidden vertices for the next step.

	Core() = default;

	Core(int k, int r, goc::Vertex v, const VertexSet& S);
};
} // namespace tdtsptw

#endif //TDTSPTW_CORE_H
