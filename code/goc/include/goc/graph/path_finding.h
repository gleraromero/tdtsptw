//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_GRAPH_PATH_FINDING_H
#define GOC_GRAPH_PATH_FINDING_H

#include <functional>

#include "goc/graph/digraph.h"
#include "goc/graph/graph_path.h"
#include "goc/graph/vertex.h"

namespace goc
{
// Precondition: D is a DAG. s, t \in Vertices(D).
// Returns: Longest path between s and t.
GraphPath longest_path(const Digraph& D, Vertex s, Vertex t);

// Precondition: tt is a travel time function which has the FIFO property.
//	- D: digraph.
//	- s: start vertex.
//	- t0: start vertex initial time.
//	- tt(i, j, t0): travel time departing from i at t0 to j (INFTY if infeasible).
// Returns: the earliest arrival time to all vertices from s.
std::vector<double> compute_earliest_arrival_time(const Digraph& D, Vertex s, double t0, const std::function<double(Vertex, Vertex, double)>& tt);

//	- D: digraph.
//	- s: start vertex.
//  - t0: start vertex initial time.
//	- dep(i, j, tf): departing time from i to reach j at tf (INFTY if impossible).
// Returns: a vector LDT, where LDT[k] is the latest time we can depart from k to reach s.
std::vector<double> compute_latest_departure_time(const Digraph& D, Vertex s, double t0, const std::function<double(Vertex, Vertex, double)>& dep);
} // namespace goc

#endif //GOC_GRAPH_PATH_FINDING_H
