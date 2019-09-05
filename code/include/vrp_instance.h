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
	goc::Matrix<goc::PWLFunction> pretau; // pretau[i][j](t) = travel time of arc (i, j) if arriving at j at t.
	goc::Matrix<goc::PWLFunction> dep; // dep[i][j](t) = departing time of arc (i, j) if arriving to j at t.
	goc::Matrix<goc::PWLFunction> arr; // arr[i][j](t) = arrival time of arc (i, j) if departing from i at t.
	goc::Matrix<TimeUnit> LDT; // LDT[i][j] = latest time we can depart from i to reach j within its tw.
	goc::Matrix<TimeUnit> EAT; // EAT[i][j] = earliest time we can arrive to j if departing from i.
	goc::Matrix<bool> prec; // prec[i][j] = i is a predecessor of j.
	std::vector<int> prec_count; // prec_count[i] = #predecessors of i.
	
	// Returns: the time we finish visiting the last vertex if departing at t0.
	// If infeasible, returns INFTY.
	TimeUnit ReadyTime(const goc::GraphPath& p, TimeUnit t0=0) const;

	// Arrival time function.
	// @return Arrival time at end of arc e, if departing at t0, counting the time windows.
	// @details if departing at t0 is infeasible it returns INFTY.
	TimeUnit ArrivalTime(goc::Arc e, TimeUnit t0) const;

	// @return Travel time for arc e, if departing at t0, counting the time windows.
	// @details if departing at t0 is infeasible it returns INFTY.
	TimeUnit TravelTime(goc::Arc e, TimeUnit t0) const;
	
	// Prints the JSON representation of the instance.
	virtual void Print(std::ostream& os) const;
};

// Serializes the instance.
void to_json(nlohmann::json& j, const VRPInstance& instance);

// Parses an instance.
void from_json(const nlohmann::json& j, VRPInstance& instance);
} // namespace tdtsptw

#endif //TDTSPTW_VRP_INSTANCE_H
