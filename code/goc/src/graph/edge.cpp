//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/graph/edge.h"

using namespace std;
using namespace goc;
using namespace nlohmann;

namespace goc
{
Edge::Edge(Vertex tail, Vertex head)
{
	this->tail = min(tail, head);
	this->head = max(tail, head);
}

bool Edge::IsIncident(Vertex node) const
{
	return node == head || node == tail;
}

bool Edge::IsAdjacent(const Edge& e) const
{
	return tail == e.tail || tail == e.head || head == e.tail || head == e.head;
}

void Edge::Print(ostream& os) const
{
	os << "(" << tail << "," << head << ")";
}

void to_json(json& j, const Edge& e)
{
	j = vector<json>();
	j.push_back(e.tail);
	j.push_back(e.head);
}

void from_json(const json& j, Edge& e)
{
	e = Edge(j[0], j[1]);
}

bool operator<(const Edge& e1, const Edge& e2)
{
	return make_pair(e1.tail, e1.head) < make_pair(e2.tail, e2.head);
}

bool operator==(const Edge& e1, const Edge& e2)
{
	return e1.tail == e2.tail && e1.head == e2.head;
}
} // namespace goc