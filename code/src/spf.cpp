//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "spf.h"

#include "pricing_problem.h"

using namespace std;
using namespace goc;

namespace tdtsptw
{
SPF::SPF(const Digraph& D) : n(D.VertexCount()), D(D)
{
	route_seq = 0;
	y_by_arc = Matrix<unordered_set<Variable>>(n, n);
	
	formulation = LPSolver::NewFormulation();
	formulation->Minimize(Expression()); // Minimization.
	formulation->AddConstraint(Expression().EQ(0.0)); // Start depot.
	for (int i = 1; i < n-1; ++i) formulation->AddConstraint(Expression().EQ(1.0)); // Customers.
	formulation->AddConstraint(Expression().EQ(0.0)); // End depot.
}

SPF::~SPF()
{
	delete formulation;
}

void SPF::AddRoute(const Route& r)
{
	// Add route r to Omega.
	int j = route_seq++;
	
	// Add variable y_j to the formulation.
	Variable y_j = formulation->AddVariable("y_" + STR(j), VariableDomain::Real, 0.0, INFTY);
	y.insert(y_j);
	omega[y_j] = r;
	
	// Check number of visits to each vertex.
	vector<int> a_ij(n, 0);
	for (int k = 1; k < (int)r.path.size()-1; ++k) a_ij[r.path[k]]++;
	
	// Set coefficient 1.0 in vertices visited.
	for (int k = 1; k < (int)r.path.size()-1; ++k) formulation->SetConstraintCoefficient(r.path[k], y_j, a_ij[r.path[k]]);
	
	// Update omega_by_arc structure.
	for (int k = 1; k < (int)r.path.size()-1; ++k) y_by_arc[r.path[k]][r.path[k+1]].insert(y_j);
	
	// Set duration(r) as c_j in the objective function.
	formulation->SetObjectiveCoefficient(y_j, r.duration);
}

void SPF::SetForbiddenArcs(const vector<Arc>& A)
{
	// Restore previously forbidden arcs.
	for (Arc e: forbidden_arcs) for (auto& y_j: y_by_arc[e.tail][e.head]) formulation->SetVariableBound(y_j, 0.0, INFTY);
	
	// All forbidden arcs must be set to 0.
	for (Arc e: A) for (auto& y_j : y_by_arc[e.tail][e.head]) formulation->SetVariableBound(y_j, 0.0, 0.0);
	
	forbidden_arcs = A;
}

const Route& SPF::RouteOf(const Variable& variable) const
{
	return omega.at(variable);
}

PricingProblem SPF::InterpretDuals(const vector<double>& duals) const
{
	PricingProblem pp;
	pp.penalties = vector<double>(duals.begin(), duals.begin()+n);
	return pp;
}

vector<Route> SPF::InterpretSolution(const Valuation& z) const
{
	vector<Route> solution;
	for (auto& y_route: omega)
		if (epsilon_equal(z[y_route.first], 1.0))
			solution.push_back(y_route.second);
	return solution;
}
} // namespace tdtsptw