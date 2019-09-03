//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/linear_programming/solver/bc_solver.h"

#include "goc/linear_programming/cplex/cplex_formulation.h"
#include "goc/linear_programming/cplex/cplex_solver.h"
#include "goc/time/duration.h"

using namespace std;
using namespace nlohmann;

namespace goc
{
BCSolver::BCSolver()
{
	time_limit = Duration::Max();
	config = {};
	screen_output = nullptr;
}

BCExecutionLog BCSolver::Solve(Formulation* formulation, const std::unordered_set<BCOption>& options) const
{
	return cplex::solve_bc((CplexFormulation*)formulation, screen_output, time_limit, config, initial_solutions,
						  branch_priorities, separation_strategy, options);
}

Formulation* BCSolver::NewFormulation()
{
	return new CplexFormulation();
}
} // namespace goc