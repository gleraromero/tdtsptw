//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/linear_programming/solver/lp_solver.h"

#include "goc/linear_programming/cplex/cplex_formulation.h"
#include "goc/linear_programming/cplex/cplex_solver.h"

using namespace std;
using namespace nlohmann;

namespace goc
{
LPSolver::LPSolver()
{
	// Set default values.
	time_limit = Duration::Max();
	config = {};
	screen_output = nullptr;
}

LPExecutionLog LPSolver::Solve(Formulation* formulation, const unordered_set<LPOption>& options) const
{
	return cplex::solve_lp((CplexFormulation*)formulation, screen_output, time_limit, config, options);
}

Formulation* LPSolver::NewFormulation()
{
	return new CplexFormulation();
}
} // namespace goc