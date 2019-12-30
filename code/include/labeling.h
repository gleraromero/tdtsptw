//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef TDTSPTW_LABELING_H
#define TDTSPTW_LABELING_H

#include <vector>
#include <unordered_map>
#include <tuple>

#include "goc/goc.h"

#include "vrp_instance.h"
#include "pricing_problem.h"

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
	goc::Matrix<goc::Interval> tw; // tw[k][v] is the smallest time window of vertex v if path length is k.

	// Builds the NGStructure with respect to vrp and delta.
	NGStructure(const VRPInstance& vrp, int delta);

	NGStructure(const VRPInstance& vrp, const std::vector<VertexSet>& N, const goc::GraphPath& L, int delta);

	void AdjustTimeWindows(const VRPInstance& vrp);
};

class NGLabel : public goc::Printable
{
public:
	NGLabel* prev;
	goc::Vertex v;
	int S;
	double Tdur, Thelp, Ttime, lambda;

	NGLabel(NGLabel* prev, goc::Vertex v, int S, double Ttime, double Tdur, double Thelp, double lambda);

	goc::GraphPath Path() const;

	virtual void Print(std::ostream& os) const;
};

double run_nglti(const VRPInstance& vrp, const NGStructure& NG, const std::vector<double>& lambda, double UB,
				 goc::Route& best_route, double& best_cost, goc::MLBExecutionLog* log);

double bidirectional_run_nglti(const VRPInstance& fvrp, const VRPInstance& bvrp, const NGStructure& fNG, const NGStructure& bNG, const std::vector<double>& lambda, double UB,
							   goc::Route& best_route, double& best_cost, goc::BLBExecutionLog* blb_log);

class State : public goc::Printable
{
public:
	class Piece : public goc::Printable
	{
	public:
		double lb;
		goc::LinearFunction f;
		Piece* prev;
		goc::Vertex v;
		
		Piece(double lb, const goc::LinearFunction& f, Piece* prev, goc::Vertex v);
		
		// Reduces domain of p2 if any prefix or suffix is dominated (this(x) <= p2(x)).
		// If p2 is fully dominated the function returns true.
		bool Dominate(Piece& p2);
		
		goc::GraphPath Path() const;
		
		virtual void Print(std::ostream& os) const;
	};
	
	// Merges pieces in F with pieces with P, preserving pieces in F when matches occur.
	// Precondition: pieces in P are disjunt and sorted by domain.
	bool Merge(std::vector<Piece>& P);
	
	void DominateBy(const State& s2);
	
	virtual void Print(std::ostream& os) const;
	
	std::vector<Piece> F;
};

struct Bounding
{
	Bounding(const VRPInstance& vrp, const NGStructure& NG, const std::vector<double>& lambda);
	
	void AddBound(int k, int r, goc::Vertex v, const VertexSet& S, const State& Delta);
	
	void Bound(goc::Vertex v, const VertexSet& S, State& Delta);
	
private:
	std::vector<std::vector<std::vector<std::vector<std::pair<VertexSet, std::vector<goc::LinearFunction>>>>>> B;
	NGStructure NG;
	VertexSet LSet;
	VRPInstance vrp;
	double Lambda;
	std::vector<double> lambda;
};

goc::Route run_dna(const VRPInstance& vrp, const VRPInstance& rvrp, NGStructure& NG, NGStructure& rNG, const std::vector<double>& lambda, goc::CGExecutionLog* log, double& LB, goc::Duration time_limit, bool bidirectional);

goc::Route run_ngltd_bidirectional(const VRPInstance& vrp, const VRPInstance& rvrp, const NGStructure& NG, const NGStructure& rNG, const std::vector<double>& lambda, goc::BLBExecutionLog* blb_log, double& LB, goc::Duration time_limit);

goc::Route run_ngltd(const VRPInstance& vrp, const NGStructure& NG, const std::vector<double>& lambda, goc::MLBExecutionLog* log, Bounding* B, double& LB, goc::Duration time_limit);

goc::Route run_exact_piecewise(const VRPInstance& vrp, const goc::GraphPath& L, const std::vector<double>& lambda,
						  double LB, double UB, goc::MLBExecutionLog* log, Bounding* B, goc::Duration time_limit);
} // namespace tdtsptw

#endif //TDTSPTW_LABELING_H
