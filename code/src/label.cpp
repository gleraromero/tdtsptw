//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "label.h"

using namespace std;
using namespace goc;

namespace tdtsptw
{
Label::Label(Label* prev, goc::Vertex v, int k, VertexSet S, goc::PWLFunction D, TimeUnit t, TimeUnit cost, double P)
	: prev(prev), v(v), k(k), S(S), D(D), t(t), cost(cost), P(P)
{

}

GraphPath Label::Path() const
{
	GraphPath p = prev ? prev->Path() : GraphPath();
	p.push_back(v);
	return p;
}

void Label::Print(ostream& os) const
{
	os << "{ P: " << Path() << ", k: " << k << ", S: " << S << ", P: " << P << ", D: " << D << "}";
}
} // namespace tdtsptw