//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_GRAPH_GRAPH_H
#define GOC_GRAPH_GRAPH_H

#include <iostream>
#include <vector>

#include "goc/collection/matrix.h"
#include "goc/graph/edge.h"
#include "goc/graph/vertex.h"
#include "goc/lib/json.hpp"
#include "goc/print/printable.h"

namespace goc
{
// This class represents a (simple, undirected) Graph.
// - The set of vertices is numbered from 0 to n.
// - Edges can be added and removed from the Graph dynamically.
class Graph : public Printable
{
public:
	// Creates a complete graph with n vertices.
	static Graph Complete(int n);
	
	// Creates a no-vertex graph.
	Graph() = default;
	
	// Creates a graph with 'vertex_count' vertices.
	Graph(int vertex_count);
	
	// Adds edge e to the graph.
	// Returns: a reference to this graph to concatenate calls.
	Graph& AddEdge(Edge e);
	
	// Adds edges to the graph.
	// Returns: a reference to this graph to concatenate calls.
	Graph& AddEdges(const std::vector<Edge>& edges);
	
	// Removes edge e from the graph.
	// Returns: a reference to this graph to concatenate calls.
	Graph& RemoveEdge(Edge e);
	
	// Returns: a vector with all the vertices ordered by number ascendingly.
	const std::vector<Vertex>& Vertices() const;
	
	// Returns: a vector with all the edges.
	const std::vector<Edge>& Edges() const;
	
	// Returns: all edges (u, v) \in E(G).
	const std::vector<Edge>& IncidentEdges(Vertex v) const;
	
	// Returns: all vertices w such that (v, w) \in E(G).
	const std::vector<Vertex>& Neighbours(Vertex v) const;
	
	// Returns: if edge e \in E(G).
	bool IncludesEdge(const Edge& e) const;
	
	// Returns: number of vertices in the graph.
	int VertexCount() const;
	
	// Returns: number of edges in the graph.
	int EdgeCount() const;
	
	// Prints the JSON serialization of the Graph.
	virtual void Print(std::ostream& os) const;

private:
	std::vector<Vertex> vertices_; // {0, ..., vertex_count - 1}
	std::vector<Edge> edges_; // E(G).
	Matrix<bool> adjacency_matrix_; // adjacency_matrix[i][j] == (i, j) \in E(G).
	std::vector<std::vector<Vertex>> adjacency_list_; // adjacency_list_[i] == {j \in V(G) : (i, j) \in E(G)}.
	std::vector<std::vector<Edge>> incident_edges_; // incident_edges_[i] = {(i, j) \in E(G) }.
};

void from_json(const nlohmann::json& j, Graph& G);

void to_json(nlohmann::json& j, const Graph& G);
} // namespace goc

#endif //GOC_GRAPH_GRAPH_H
