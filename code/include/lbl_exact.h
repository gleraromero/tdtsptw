//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef TDTSPTW_LBL_EXACT_H
#define TDTSPTW_LBL_EXACT_H

#include <vector>
#include <unordered_map>
#include <tuple>

#include "goc/goc.h"

#include "vrp_instance.h"
#include "pricing_problem.h"
#include "dssr.h"

namespace tdtsptw
{
// Runs an NG labeling algorithm to find negative cost routes.
// @param vrp: VRP Instance
// @param NG: structure with NG information.
// @param B: bounding structure for completion bounds.
// @param lambda: penalties for vertices.
// @param UB: upper bound on the duration of the optimal route.
// @param [out] best_route: route with best cost.
// @param [out] log: output log to save the execution information.
// @returns the optimal tour.
goc::Route run_exact(const VRPInstance& vrp, const NGStructure& NG, BoundingStructure& B,
					 const std::vector<double>& lambda, const goc::Route& UB, double lb, goc::MLBExecutionLog* log);

struct State
{
public:
	class Piece : public goc::Printable
	{
	public:
		double lb;
		goc::LinearFunction f;
		
		Piece(double lb, const goc::LinearFunction& f);
		
		// Reduces domain of p2 if any prefix or suffix is dominated (this(x) <= p2(x)).
		// If p2 is fully dominated the function returns true.
		bool Dominate(Piece& p2);
		
		virtual void Print(std::ostream& os) const;
	};
	
	// Merges pieces in F with pieces with P, preserving pieces in F when matches occur.
	// Precondition: pieces in P are disjunt and sorted by domain.
	void Merge(std::vector<Piece>& P);
	
	void DominateBy(const State& s2);
	
	std::vector<Piece> F;
};

struct Bounding
{
	Bounding(const VRPInstance& vrp, const NGStructure& NG, const std::vector<double>& lambda);
	
	void AddBound(int k, int r, goc::Vertex v, const VertexSet& S, const State& Delta);
	
	void Bound(goc::Vertex v, VertexSet S, State& Delta);
	
private:
	std::vector<std::vector<std::vector<std::vector<std::pair<VertexSet, std::vector<goc::LinearFunction>>>>>> B;
	NGStructure NG;
	VertexSet LSet;
	VRPInstance vrp;
	double Lambda;
	std::vector<double> lambda;
};

goc::Route run_dssr(const VRPInstance& vrp, NGStructure& NG, const std::vector<double>& lambda, goc::CGExecutionLog* log, double& LB);

goc::Route run_ngl(const VRPInstance& vrp, const NGStructure& NG, const std::vector<double>& lambda, goc::MLBExecutionLog* log, Bounding* B, double& LB);

goc::Route run_exact_piecewise(const VRPInstance& vrp, const goc::GraphPath& L, const std::vector<double>& lambda,
						  double LB, double UB, goc::MLBExecutionLog* log, Bounding* B);
} // namespace tdtsptw

#endif //TDTSPTW_LBL_EXACT_H
