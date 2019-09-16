//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef TDTSPTW_SPF_H
#define TDTSPTW_SPF_H

#include <vector>

#include "goc/goc.h"

#include "vrp_instance.h"

namespace tdtsptw
{
class PricingProblem;

// Represents a set-partitioning formulation for the VRP.
// min c_j y_j															(0)
// s.t.
//		sum_{j \in \Omega} a_ij y_j = 1 	\forall i \in V - {o, d}	(1)
//		y_j >= 0							\forall j \in \Omega		(2)
class SPF
{
public:
	goc::Formulation* formulation; // formulation for the restricted master problem
	std::unordered_set<goc::Variable> y; // variables associated with routes in omega.
	
	// Creates a Set-Partitioning formulation for VRP with n vertices.
	// Assumes start-depot is 0 and end depot is n-1.
	explicit SPF(int n);
	
	// Frees formulation memory.
	~SPF();
	
	// Adds the specified route to the formulation (using the duration as its cost).
	void AddRoute(const goc::Route& r);
	
	// Returns: the route associated to the variable.
	const goc::Route& RouteOf(const goc::Variable& variable) const;
	
	// Returns: the pricing problem for the duals.
	PricingProblem InterpretDuals(const std::vector<double>& duals) const;
	
	// Interprets an integer solution of the formulation and returns the routes selected.
	std::vector<goc::Route> InterpretSolution(const goc::Valuation& z) const;

private:
	int n; // number of vertices.
	std::unordered_map<goc::Variable, goc::Route> omega; // set of routes in the restricted master problem.
};
} // namespace tdtsptw

#endif //TDTSPTW_SPF_H
