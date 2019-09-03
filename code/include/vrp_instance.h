//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef TDTSPTW_VRP_INSTANCE_H
#define TDTSPTW_VRP_INSTANCE_H

#include <vector>
#include <goc/goc.h>

// Maximum number of vertices in an instance.
#define MAX_N 100

namespace tdtsptw
{
typedef double TimeUnit; // Represents time.
typedef std::bitset<MAX_N> VertexSet;

// This class represents an instance of a vehicle routing problem.
// Considerations:
// 	- It considers two depots (origin and destination).
class VRPInstance : public goc::Printable
{
public:
	goc::Digraph D; // digraph representing the network.
	goc::Vertex o, d; // origin and destination depot.
	TimeUnit T; // end of planning horizon ([0,T]).
	std::vector<goc::Interval> tw; // time window of customers (tw[i] = time window of customer i).
	goc::Matrix<goc::PWLFunction> tau; // tau[i][j](t) = travel time of arc (i, j) if departing from i at t.
	
	// Returns: the time we finish visiting the last vertex if departing at t0.
	// If infeasible, returns INFTY.
	TimeUnit ReadyTime(const goc::GraphPath& p, TimeUnit t0=0) const;
	
	// Prints the JSON representation of the instance.
	virtual void Print(std::ostream& os) const;
};

// Serializes the instance.
void to_json(nlohmann::json& j, const VRPInstance& instance);

// Parses an instance.
void from_json(const nlohmann::json& j, VRPInstance& instance);
} // namespace tdtsptw

#endif //TDTSPTW_VRP_INSTANCE_H
