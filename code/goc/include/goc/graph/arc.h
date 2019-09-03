//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_GRAPH_ARC_H
#define GOC_GRAPH_ARC_H

#include <iostream>

#include "goc/graph/vertex.h"
#include "goc/lib/json.hpp"
#include "goc/print/printable.h"

namespace goc
{
// Represents a directed arc in a digraph.
// - (i, j) has a tail = i, head = j.
class Arc : public Printable
{
public:
	Vertex tail, head;
	
	// Creates an arc (tail, head).
	Arc(Vertex tail, Vertex head);
	
	// Returns: if node \in {tail, head}.
	bool IsIncident(Vertex node) const;
	
	// Returns: if head(this) == tail(arc).
	bool IsPredecessorOf(const Arc& arc) const;
	
	// Returns: if tail(this) == head(arc).
	bool IsSuccessorOf(const Arc& arc) const;
	
	// Returns: an arc (head, tail).
	Arc Reverse() const;
	
	// Prints the arc in the format: "(tail, head)".
	virtual void Print(std::ostream& os) const;
};

void to_json(nlohmann::json& j, const Arc& e);

void from_json(const nlohmann::json& j, Arc& e);

// Returns: if (tail(e1), head(e1)) < (tail(e2), head(e2)) as a pair-wise comparison.
bool operator<(const Arc& e1, const Arc& e2);

// Returns: if tail(e1)==tail(e2) and head(e1)==head(e2).
bool operator==(const Arc& e1, const Arc& e2);
} // namespace goc

// Implement hash function in order to use Arc as keys in unordered_map and unordered_set.
namespace std
{
template<>
class hash<goc::Arc>
{
public:
	size_t operator()(const goc::Arc& v) const
	{
		return std::hash<int>()(v.tail) ^ std::hash<int>()(v.head);
	}
};
} //namespace std

#endif //GOC_GRAPH_ARC_H
