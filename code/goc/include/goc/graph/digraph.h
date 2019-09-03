//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_GRAPH_DIGRAPH_H
#define GOC_GRAPH_DIGRAPH_H

#include <iostream>
#include <vector>

#include "goc/collection/matrix.h"
#include "goc/graph/arc.h"
#include "goc/graph/vertex.h"
#include "goc/lib/json.hpp"
#include "goc/print/printable.h"

namespace goc
{
// This class represents a (simple, directed) Graph.
// - The set of vertices is numbered from 0 to n.
// - Arcs can be added and removed from the Digraph dynamically.
class Digraph : public Printable
{
public:
	// Creates a complete digraph with n vertices.
	static Digraph Complete(int n);
	
	// Creates a no-vertex digraph.
	Digraph() = default;
	
	// Creates a digraph with 'vertex_count' vertices.
	Digraph(int vertex_count);
	
	// Adds arc e to the digraph.
	// Returns: a reference to this digraph to concatenate calls.
	Digraph& AddArc(Arc e);
	
	// Adds arcs e to the digraph.
	// Returns: a reference to this digraph to concatenate calls.
	Digraph& AddArcs(const std::vector<Arc>& arcs);
	
	// Removes arc e from the digraph.
	// Returns: a reference to this digraph to concatenate calls.
	Digraph& RemoveArc(Arc e);
	
	// Removes arcs from the digraph.
	// Returns: a reference to this digraph to concatenate calls.
	Digraph& RemoveArcs(const std::vector<Arc>& arcs);
	
	// Returns: a vector with all the vertices ordered by number ascendingly.
	const std::vector<Vertex>& Vertices() const;
	
	// Returns: a vector with all the arcs.
	const std::vector<Arc>& Arcs() const;
	
	// Returns: all arcs (u, v) \in D.
	const std::vector<Arc>& InboundArcs(Vertex v) const;
	
	// Returns: all arcs (v, w) \in D.
	const std::vector<Arc>& OutboundArcs(Vertex v) const;
	
	// Returns: all vertices w such that (v, w) \in A(D).
	const std::vector<Vertex>& Successors(Vertex v) const;
	
	// Returns: all vertices u such that (u, v) \in A(D).
	const std::vector<Vertex>& Predecessors(Vertex v) const;
	
	// Returns: if arc e \in A(D).
	bool IncludesArc(const Arc& e) const;
	
	// Returns: number of vertices in the digraph.
	int VertexCount() const;
	
	// Returns: number of arcs in the digraph.
	int ArcCount() const;
	
	// Returns: this digraph reversed, meaning that an arc (j, i) \in A(reverse(D)) iif (i, j) \in A(D).
	Digraph Reverse() const;
	
	// Prints the JSON serialization of the Digraph.
	virtual void Print(std::ostream& os) const;
	
private:
	Matrix<bool> adjacency_matrix_; // adjacency_matrix[i][j] == (i, j) \in A(D).
	std::vector<std::vector<Vertex>> successor_list_; // succesor_list[i] == {j \in V(D) : (i, j) \in A(D)}.
	std::vector<std::vector<Vertex>> predecessor_list_; // predecessor_list[j] == {i \in V(D) : (i, j) \in A(D)}.
	std::vector<Vertex> vertices_; // {0, ..., vertex_count - 1}
	std::vector<Arc> arcs_; // A(D).
	std::vector<std::vector<Arc>> inbound_arcs_; // inbound_arcs[j] = {(i, j) \in A(D) }.
	std::vector<std::vector<Arc>> outbound_arcs_; // outbound_arcs[i] = {(i, j) \in A(D) }.
};

void from_json(const nlohmann::json& j, Digraph& D);

void to_json(nlohmann::json& j, const Digraph& D);
} // namespace goc

#endif //GOC_GRAPH_DIGRAPH_H
