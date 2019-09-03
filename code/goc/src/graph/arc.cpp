//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/graph/arc.h"

using namespace std;
using namespace goc;
using namespace nlohmann;

namespace goc
{
Arc::Arc(Vertex tail, Vertex head) : tail(tail), head(head)
{
}

bool Arc::IsIncident(Vertex node) const
{
	return node == head || node == tail;
}

bool Arc::IsPredecessorOf(const Arc& arc) const
{
	return head == arc.tail;
}

bool Arc::IsSuccessorOf(const Arc& arc) const
{
	return tail == arc.head;
}

Arc Arc::Reverse() const
{
	return {head, tail};
}

void Arc::Print(ostream& os) const
{
	os << "(" << tail << "," << head << ")";
}

void to_json(json& j, const Arc& e)
{
	j = vector<json>();
	j.push_back(e.tail);
	j.push_back(e.head);
}

void from_json(const json& j, Arc& e)
{
	e = Arc(j[0], j[1]);
}

bool operator<(const Arc& e1, const Arc& e2)
{
	return make_pair(e1.tail, e1.head) < make_pair(e2.tail, e2.head);
}

bool operator==(const Arc& e1, const Arc& e2)
{
	return e1.tail == e2.tail && e1.head == e2.head;
}
} // namespace goc