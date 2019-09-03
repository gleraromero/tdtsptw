//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_LINEAR_PROGRAMMING_CPLEX_CPLEX_SOLVER_H
#define GOC_LINEAR_PROGRAMMING_CPLEX_CPLEX_SOLVER_H

#include <iostream>
#include <string>
#include <unordered_set>

#include "goc/linear_programming/cplex/cplex_formulation.h"
#include "goc/linear_programming/cuts/separation_strategy.h"
#include "goc/linear_programming/model/branch_priority.h"
#include "goc/linear_programming/solver/lp_solver.h"
#include "goc/linear_programming/solver/bc_solver.h"
#include "goc/log/bc_execution_log.h"
#include "goc/log/lp_execution_log.h"

namespace goc
{
namespace cplex
{
// Solves the formulation using the CPLEX lpopt solver.
//	formulation: lp model to be solved.
//	screen_output: stream where the cplex logs should be outputted (nullptr if no output is desired).
//	time_limit: time limit for CPLEX lpopt.
// 	config: map with the CPLEX parameter names as keys and their values.
// 	options: which options should be returned in the execution log.
// Returns: the execution log with the options specified in log_options.
LPExecutionLog solve_lp(CplexFormulation* formulation,
						std::ostream* screen_output,
						Duration time_limit,
						const nlohmann::json& config,
						const std::unordered_set<LPOption>& options);

// Solves the formulation using the CPLEX mipopt solver.
//	formulation: lp model to be solved.
//	screen_output: stream where the cplex logs should be outputted (nullptr if no output is desired).
//	time_limit: time limit for CPLEX lpopt.
// 	config: map with the CPLEX parameter names as keys and their values.
// 	initial_solutions: a sequence of initial solutions that should be added as MIPstarts to CPLEX.
// 	branch_priorities: a sequence of branch hints that should be given to CPLEX to help reduce the BB tree.
//	separation_strategy: the separation algorithm that will be called at every relaxation to add cuts.
// 	options: which options should be returned in the execution log.
// Returns: the execution log with the options specified in log_options.
BCExecutionLog solve_bc(CplexFormulation* formulation,
				std::ostream* screen_output,
				Duration time_limit,
				const nlohmann::json& config,
				const std::vector<Valuation>& initial_solutions,
				const std::vector<BranchPriority>& branch_priorities,
				const goc::SeparationStrategy& separation_strategy,
				const std::unordered_set<BCOption>& options);
} // namespace cplex
} // namespace goc

#endif //GOC_LINEAR_PROGRAMMING_CPLEX_CPLEX_SOLVER_H
