//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_LINEAR_PROGRAMMING_SOLVER_CG_SOLVER_H
#define GOC_LINEAR_PROGRAMMING_SOLVER_CG_SOLVER_H

#include <functional>
#include <iostream>
#include <unordered_set>
#include <vector>

#include "goc/linear_programming/model/formulation.h"
#include "goc/linear_programming/solver/lp_solver.h"
#include "goc/log/cg_execution_log.h"
#include "goc/time/duration.h"

namespace goc
{
// All the log options that can be enabled/disabled.
// - IterationsInformation: if not included, {iterations} will not be stored.
//							advantage: saving the space of the iteration logs.
// - ScreenOutput: 			if not included, the output will not be stored.
//							advantage: saving space.
enum class CGOption { IterationsInformation, ScreenOutput };

// Function type for the pricing solver.
// - duals: the dual variables of the iteration.
// - incumbent_value: the value of the lp relaxation.
// - time_limit: maximum time to execute the pricing algorithm.
// - execution_log: pointer to the cg execution log to add the iteration log.
typedef std::function<void(const std::vector<double>& duals, double incumbent_value, Duration time_limit, CGExecutionLog* cg_execution_log)> PricingFunction;

// Class representing a solver for column generation. Its purpose is to abstract the
// specific solver implementations from the algorithms.
class CGSolver
{
public:
	// Pointer to the stream where the output of the algorithm should go. (nullptr for no output).
	std::ostream* screen_output;
	// Maximum time to spend solving.
	Duration time_limit;
	// Pointer to a solver for the linear relaxation of the lp.
	goc::LPSolver* lp_solver;
	// A function that receives the current lp iteration and adds variables to the lp. The column generation will
	// continue as long as the pricing function adds variables (or constraints) to the lp at a given iteration.
	PricingFunction pricing_function;
	
	// Creates a default column generation solver (no output, time_limit=2hs, lp_solver=CPLEX,
	// 	pricing_function=DONOTHING).
	CGSolver();
	
	// Solves the formulation using a column generation procedure.
	// Returns: the execution log with the specified options.
	// Precondition: the formulation must have been created with the NewFormulation() method.
	CGExecutionLog Solve(Formulation* formulation, const std::unordered_set<CGOption>& options={}) const;
	
	// Returns: a formulation compatible with the solver.
	static Formulation* NewFormulation();
};
} // namespace goc

#endif //GOC_LINEAR_PROGRAMMING_SOLVER_CG_SOLVER_H
