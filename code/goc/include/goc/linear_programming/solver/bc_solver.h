//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_LINEAR_PROGRAMMING_SOLVER_BC_SOLVER_H
#define GOC_LINEAR_PROGRAMMING_SOLVER_BC_SOLVER_H

#include <iostream>
#include <unordered_set>
#include <vector>

#include "goc/lib/json.hpp"
#include "goc/linear_programming/cuts/separation_strategy.h"
#include "goc/linear_programming/model/branch_priority.h"
#include "goc/linear_programming/model/formulation.h"
#include "goc/linear_programming/model/valuation.h"
#include "goc/log/bc_execution_log.h"
#include "goc/time/duration.h"

namespace goc
{
// All the log options that can be enabled/disabled.
// - ScreenOutput: 		if not included, the output will not be stored.
//						advantage: saving space.
// - RootInformation: 	if not included {root_lp_value, root_int_value, root_int_solution} will not be filled.
//						advantage: is that no callbacks need to be set to retrieve this information.
// - BestIntSolution:	if not included {best_int_solution} will not be filled.
//						advantage: if solution has many variables, getting it is linear in that size.
// - CutInformation: 	if not included {cut_count, cut_iteration_count, cut_time, cut_families, cut_family_cut_count,
//						cut_family_iteration_count, cut_family_cut_time} will not be filled.
//						advantage: saving space.
enum class BCOption {
	ScreenOutput, RootInformation, BestIntSolution, CutInformation
};

// Class representing a solver for branch and cut. Its purpose is to abstract the
// specific solver implementations from the algorithms.
class BCSolver
{
public:
	// Pointer to the stream where the output of the algorithm should go. (nullptr for no output).
	std::ostream* screen_output;
	// Maximum time to spend solving.
	Duration time_limit;
	// Json object with the configuration options to send to the solver.
	nlohmann::json config;
	// Object that indicates what families of cuts will be added and the strategy to do so.
	SeparationStrategy separation_strategy;
	// A set of initial solutions for the BC.
	std::vector<Valuation> initial_solutions;
	// Each variable might receive a number containing the priority on how important is to branch on
	// that variable. The higher the priority the earliest the variable will be selected to be branched
	std::vector<BranchPriority> branch_priorities;
	
	// Creates a default branch and cut solver. (time limit: 2 hours).
	// Currently: CPLEX.
	BCSolver();
	
	// Solves the formulation.
	// Returns: the execution log with the specified options.
	// Precondition: the formulation must have been created with the NewFormulation() method.
	BCExecutionLog Solve(Formulation* formulation, const std::unordered_set<BCOption>& options={}) const;
	
	// Returns: a formulation compatible with the solver.
	static Formulation* NewFormulation();
};
} // namespace goc

#endif //GOC_LINEAR_PROGRAMMING_SOLVER_BC_SOLVER_H
