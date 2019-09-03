//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_GRAPH_MAXFLOW_MINCUT_H
#define GOC_GRAPH_MAXFLOW_MINCUT_H

#include <functional>
#include <vector>

#include "goc/graph/digraph.h"

namespace goc
{
// Represents an S-T cut in a network. S are the vertices in the side of the source,
// T the vertices in the side of the sink.
struct STCut
{
	std::vector<int> S, T;
};

// Solves a maxflow problem on the network D with capacities c with source s and sink t.
// Parameters:
// D: digraph representing the network.
// c: function that gives the capacities c(i,j) of each arc ij.
// s: source vertex.
// t: sink vertex.
// Returns (max_flow, min_cut) of the network.
std::pair<double, STCut> maxflow_mincut(const Digraph& D, const std::function<double(Vertex i, Vertex j)>& c, int s, int t);
} // namespace goc.

#endif //GOC_GRAPH_MAXFLOW_MINCUT_H