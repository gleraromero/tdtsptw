//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/graph/digraph.h"

#include "goc/collection/collection_utils.h"

using namespace std;
using namespace nlohmann;

namespace goc
{
Digraph Digraph::Complete(int n)
{
	Digraph D(n);
	for (int i: range(0, n))
		for (int j: range(i+1,n))
			D.AddArc({i,j});
	return D;
}

Digraph::Digraph(int vertex_count)
{
	vertices_ = range(0, vertex_count);
	inbound_arcs_.assign(vertex_count, vector<Arc>());
	outbound_arcs_.assign(vertex_count, vector<Arc>());
	successor_list_.assign(vertex_count, vector<Vertex>());
	predecessor_list_.assign(vertex_count, vector<Vertex>());
	adjacency_matrix_ = Matrix<bool>(vertex_count, vertex_count, false);
}

Digraph& Digraph::AddArc(Arc e)
{
	if (IncludesArc(e)) return *this;;
	arcs_.push_back(e);
	adjacency_matrix_[e.tail][e.head] = true;
	inbound_arcs_[e.head].push_back(e);
	outbound_arcs_[e.tail].push_back(e);
	successor_list_[e.tail].push_back(e.head);
	predecessor_list_[e.head].push_back(e.tail);
	return *this;
}

Digraph& Digraph::AddArcs(const std::vector<Arc>& arcs)
{
	for (Arc e: arcs) AddArc(e);
	return *this;
}

Digraph& Digraph::RemoveArc(Arc e)
{
	if (!IncludesArc(e)) return *this;
	arcs_.erase(find(arcs_.begin(), arcs_.end(), e));
	adjacency_matrix_[e.tail][e.head] = false;
	inbound_arcs_[e.head].erase(find(inbound_arcs_[e.head].begin(), inbound_arcs_[e.head].end(), e));
	outbound_arcs_[e.tail].erase(find(outbound_arcs_[e.tail].begin(), outbound_arcs_[e.tail].end(), e));
	successor_list_[e.tail].erase(find(successor_list_[e.tail].begin(), successor_list_[e.tail].end(), e.head));
	predecessor_list_[e.head].erase(find(predecessor_list_[e.head].begin(), predecessor_list_[e.head].end(), e.tail));
	return *this;
}

Digraph& Digraph::RemoveArcs(const vector<Arc>& arcs)
{
	for (Arc e: arcs) RemoveArc(e);
	return *this;
}

const vector<Vertex>& Digraph::Vertices() const
{
	return vertices_;
}

const vector<Arc>& Digraph::Arcs() const
{
	return arcs_;
}

const vector<Arc>& Digraph::InboundArcs(Vertex v) const
{
	return inbound_arcs_[v];
}

const vector<Arc>& Digraph::OutboundArcs(Vertex v) const
{
	return outbound_arcs_[v];
}

const vector<Vertex>& Digraph::Successors(Vertex v) const
{
	return successor_list_[v];
}

const vector<Vertex>& Digraph::Predecessors(Vertex v) const
{
	return predecessor_list_[v];
}

bool Digraph::IncludesArc(const Arc& e) const
{
	return adjacency_matrix_[e.tail][e.head];
}

int Digraph::VertexCount() const
{
	return (int) vertices_.size();
}

int Digraph::ArcCount() const
{
	return (int) arcs_.size();
}

Digraph Digraph::Reverse() const
{
	Digraph reverse_graph(VertexCount());
	for (auto& e: Arcs()) reverse_graph.AddArc(e.Reverse());
	return reverse_graph;
}

void Digraph::Print(ostream& os) const
{
	os << json(*this);
}

void from_json(const json& j, Digraph& D)
{
	int n = j["vertex_count"];
	D = Digraph(n);
	auto& arcs_json = j["arcs"];
	for (int i = 0; i < n; ++i)
		for (int k = 0; k < n; ++k)
			if (arcs_json[i][k] == 1)
				D.AddArc({i,k});
}

void to_json(json& j, const Digraph& D)
{
	j["vertex_count"] = D.VertexCount();
	j["arc_count"] = D.ArcCount();
	// Build adjacency matrix.
	Matrix<int> M(D.VertexCount(), D.VertexCount(), 0);
	for (int i = 0; i < D.VertexCount(); ++i)
		for (int j = 0; j < D.VertexCount(); ++j)
			M[i][j] = D.IncludesArc({i,j}) ? 1 : 0;
	j["arcs"] = M;
}
} // namespace goc