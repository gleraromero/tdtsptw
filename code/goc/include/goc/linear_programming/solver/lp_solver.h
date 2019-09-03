//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_LINEAR_PROGRAMMING_SOLVER_LP_SOLVER_H
#define GOC_LINEAR_PROGRAMMING_SOLVER_LP_SOLVER_H

#include <iostream>
#include <unordered_set>
#include <vector>

#include "goc/lib/json.hpp"
#include "goc/linear_programming/model/formulation.h"
#include "goc/log/lp_execution_log.h"
#include "goc/time/duration.h"

namespace goc
{
// All the log options that can be enabled/disabled.
// - ScreenOutput: 	if not included, the output will not be stored.
//					advantage: saving space.
// - Duals: 		if not included {duals} will not be filled.
//					advantage: if model has many constraints, getting them is linear in that size.
// - Incumbent:		if not included {incumbent} will not be filled.
//					advantage: if solution has many variables, getting it is linear in that size.
enum class LPOption { ScreenOutput, Duals, Incumbent };

// Class representing a solver for the lp relaxation. Its purpose is to abstract the
// specific solver implementations from the algorithms.
class LPSolver
{
public:
	// Pointer to the stream where the output of the algorithm should go. (nullptr for no output).
	std::ostream* screen_output;
	// maximum time to spend solving.
	Duration time_limit;
	// json object with the configuration options to send to the solver.
	nlohmann::json config;
	
	// Creates a default lp solver. (time limit: 2 hours).
	// Currently: CPLEX.
	LPSolver();
	
	// Solves the formulation.
	// Returns: the execution log with the specified options.
	// Precondition: the formulation must have been created with the NewFormulation() method.
	LPExecutionLog Solve(Formulation* formulation, const std::unordered_set<LPOption>& options={}) const;
	
	// Returns: a formulation compatible with the solver.
	static Formulation* NewFormulation();
};
} // namespace goc

#endif //GOC_LINEAR_PROGRAMMING_SOLVER_LP_SOLVER_H
