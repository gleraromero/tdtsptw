//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_GRAPH_EDGE_H
#define GOC_GRAPH_EDGE_H

#include <iostream>

#include "goc/graph/vertex.h"
#include "goc/lib/json.hpp"
#include "goc/print/printable.h"

namespace goc
{
// Represents an undirected edge in a graph.
// - (i, j) has a tail = i and a head = j.
// - Invariant: i <= j.
class Edge : public Printable
{
public:
	Vertex tail, head;
	
	// Creates an edge (tail, head).
	Edge(Vertex tail, Vertex head);
	
	// Returns: if node \in {tail, head}.
	bool IsIncident(Vertex node) const;
	
	// Returns: if the two edges share some end point.
	bool IsAdjacent(const Edge& e) const;
	
	// Prints the arc in the format: "(tail, head)".
	virtual void Print(std::ostream& os) const;
};

void to_json(nlohmann::json& j, const Edge& e);

void from_json(const nlohmann::json& j, Edge& e);

// Returns: if (tail(e1), head(e1)) < (tail(e2), head(e2)) as a pair-wise comparison.
bool operator<(const Edge& e1, const Edge& e2);

// Returns: if tail(e1)==tail(e2) and head(e1)==head(e2).
bool operator==(const Edge& e1, const Edge& e2);
} // namespace goc

// Implement hash function in order to use Edge as keys in unordered_map and unordered_set.
namespace std
{
template<>
class hash<goc::Edge> {
public:
	size_t operator()(const goc::Edge& v) const
	{
		return std::hash<int>()(v.tail) ^ std::hash<int>()(v.head);
	}
};
} // namespace std

#endif //GOC_GRAPH_EDGE_H
