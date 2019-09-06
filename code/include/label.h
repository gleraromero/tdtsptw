//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef TDTSPTW_LABEL_H
#define TDTSPTW_LABEL_H

#include "goc/goc.h"
#include "vrp_instance.h"

namespace tdtsptw
{
class Label : public goc::Printable
{
public:
	Label* prev; // parent label.
	goc::Vertex v; // last_vertex.
	int k; // path length (number of vertices).
	VertexSet S; // forbidden vertices for extension.
	goc::PWLFunction D; // D(t) minimum duration if reaching v at t.
	TimeUnit t; // t = min(dom(D))
	TimeUnit cost; // cost = min(img(D)) - P
	double P; // sum of penalty costs
	int r; // number of vertices from the L path that where visited (NGL relaxation).
	
	Label(Label* prev, goc::Vertex v, int k, VertexSet S, goc::PWLFunction D, TimeUnit t, TimeUnit cost, double P, int r=0);
	
	goc::GraphPath Path() const;
	
	virtual void Print(std::ostream& os) const;
};
} // namespace tdtsptw

#endif //TDTSPTW_LABEL_H
