//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/graph/graph.h"

#include "goc/collection/collection_utils.h"

using namespace std;
using namespace nlohmann;

namespace goc
{
Graph Graph::Complete(int n)
{
	Graph G(n);
	for (int i = 0; i < n; ++i)
		for (int j = i+1; j < n; ++j)
			G.AddEdge({i,j});
	return G;
}

Graph::Graph(int vertex_count)
{
	vertices_ = range(0, vertex_count);
	edges_ = {};
	adjacency_matrix_ = Matrix<bool>(vertex_count, vertex_count, false);
	adjacency_list_.assign(vertex_count, {});
	incident_edges_.assign(vertex_count, vector<Edge>());
}

Graph& Graph::AddEdge(Edge e)
{
	if (IncludesEdge(e)) return *this;;
	edges_.push_back(e);
	adjacency_matrix_[e.tail][e.head] = adjacency_matrix_[e.head][e.tail] = true;
	incident_edges_[e.head].push_back(e);
	incident_edges_[e.tail].push_back(e);
	adjacency_list_[e.head].push_back(e.tail);
	adjacency_list_[e.tail].push_back(e.head);
	return *this;
}

Graph& Graph::AddEdges(const std::vector<Edge>& edges)
{
	for (auto& e: edges) AddEdge(e);
	return *this;
}

Graph& Graph::RemoveEdge(Edge e)
{
	if (!IncludesEdge(e)) return *this;
	edges_.erase(find(edges_.begin(), edges_.end(), e));
	adjacency_matrix_[e.tail][e.head] = adjacency_matrix_[e.head][e.tail] = false;
	incident_edges_[e.head].erase(find(incident_edges_[e.head].begin(), incident_edges_[e.head].end(), e));
	incident_edges_[e.tail].erase(find(incident_edges_[e.tail].begin(), incident_edges_[e.tail].end(), e));
	adjacency_list_[e.head].erase(find(adjacency_list_[e.head].begin(), adjacency_list_[e.head].end(), e.tail));
	adjacency_list_[e.tail].erase(find(adjacency_list_[e.tail].begin(), adjacency_list_[e.tail].end(), e.head));
	return *this;
}

const vector<Vertex>& Graph::Vertices() const
{
	return vertices_;
}

const vector<Edge>& Graph::Edges() const
{
	return edges_;
}

const vector<Edge>& Graph::IncidentEdges(Vertex v) const
{
	return incident_edges_[v];
}

const vector<Vertex>& Graph::Neighbours(Vertex v) const
{
	return adjacency_list_[v];
}

bool Graph::IncludesEdge(const Edge& e) const
{
	return adjacency_matrix_[e.tail][e.head];
}

int Graph::VertexCount() const
{
	return (int) vertices_.size();
}

int Graph::EdgeCount() const
{
	return (int) edges_.size();
}

void Graph::Print(ostream& os) const
{
	os << json(*this);
}

void from_json(const json& j, Graph& G)
{
	int n = j["vertex_count"];
	G = Graph(n);
	auto& arcs_json = j["edges"];
	for (int i = 0; i < n; ++i)
		for (int k = 0; k < n; ++k)
			if (arcs_json[i][k] == 1)
				G.AddEdge({i,k});
}

void to_json(json& j, const Graph& G)
{
	j["vertex_count"] = G.VertexCount();
	j["edge_count"] = G.EdgeCount();
	
	// Build adjacency matrix.
	Matrix<int> M(G.VertexCount(), G.VertexCount(), 0);
	for (int i = 0; i < G.VertexCount(); ++i)
		for (int j = 0; j < G.VertexCount(); ++j)
			M[i][j] = G.IncludesEdge({i,j}) ? 1 : 0;
	j["edges"] = M;
}
} // namespace goc