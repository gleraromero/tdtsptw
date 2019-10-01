//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef TDTSPTW_LBL_NG_H
#define TDTSPTW_LBL_NG_H

#include <vector>
#include <tuple>

#include "goc/goc.h"

#include "vrp_instance.h"

namespace tdtsptw
{
// A structure that precomputes all the information necessary to execute an NG algorithm.
// @param delta: maximum size of the NG neighbourhoods.
// @param [out] N: NG neighbourhoods.
// @param [out] L: L path from the NGL path relaxation.
// @param [out] V: V[r] are the vertices that may be visited after L[r] and before L[r+1].
// @param [out] NGSet: NGSet[v][i] is the vertex set of the NGSet with index i of vertex v.
// @param [out] NGSub: NGSub[v][i] are the indices of the NGSets which are subsets of NGSet[v][i].
// @param [out] NGArc: NGArc[v][i] is a dictionary which given a destination vertex w, it returns the NGSet of w if departing
// 				 from v with NGSet[v][i].
struct NGStructure
{
	int delta;
	std::vector<VertexSet> N;
	goc::GraphPath L;
	std::vector<VertexSet> V;
	std::vector<std::vector<VertexSet>> NGSet;
	std::vector<std::vector<std::vector<int>>> NGSub;
	std::vector<std::vector<std::unordered_map<goc::Vertex, int>>> NGArc;
	
	// Builds the NGStructure with respect to vrp and delta.
	NGStructure(const VRPInstance& vrp, int delta);
	
	NGStructure(const VRPInstance& vrp, const std::vector<VertexSet>& N, const goc::GraphPath& L, int delta);
};

class NGLabel : public goc::Printable
{
public:
	NGLabel* prev;
	goc::Vertex v;
	int S;
	double Tdur, Thelp, lambda;
	
	NGLabel(NGLabel* prev, goc::Vertex v, int S, double Tdur, double Thelp, double lambda);
	
	goc::GraphPath Path() const;
	
	virtual void Print(std::ostream& os) const;
};

// Runs an NG labeling algorithm to find negative cost routes.
// @param vrp: VRP Instance
// @param NG: structure with NG information.
// @param lambda: penalties for vertices.
// @param UB: upper bound on the duration of the optimal route.
// @param [out] best_route: route with best cost.
// @param [out] best_cost: cost of best_route.
// @param [out] log: output log to save the execution information.
// @returns a set of negative cost routes.
std::vector<goc::Route> run_ng(const VRPInstance& vrp, const NGStructure& NG, const std::vector<double>& lambda,
							   double UB, goc::Route* best_route, double* best_cost, goc::MLBExecutionLog* log);
} // namespace tdtsptw

#endif //TDTSPTW_LBL_NG_H
