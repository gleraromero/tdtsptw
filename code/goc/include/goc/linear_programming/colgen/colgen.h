//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_LINEAR_PROGRAMMING_COLGEN_COLGEN_H
#define GOC_LINEAR_PROGRAMMING_COLGEN_COLGEN_H

#include <functional>
#include <iostream>
#include <unordered_set>
#include <vector>

#include "goc/linear_programming/model/formulation.h"
#include "goc/linear_programming/solver/lp_solver.h"
#include "goc/linear_programming/solver/cg_solver.h"
#include "goc/log/cg_execution_log.h"

namespace goc
{
// formulation: Master problem formulation.
// screen_output: Stream where the output to screen should be sent.
// time_limit: Maximum time that might be spent on this method.
// pricing_function: Function that given a set of dual variables and the objective value finds and adds entering variables to the master formulation base.
// lp_solver: Linear relaxation solver.
// options: Which options of the execution to keep track of.
// Returns: the execution log of the column generation with the specified options.
CGExecutionLog solve_colgen(Formulation* formulation,
				   std::ostream* screen_output,
				   Duration time_limit,
				   const PricingFunction& pricing_function,
				   LPSolver* lp_solver,
				   const std::unordered_set<CGOption>& options);
} // namespace goc

#endif //GOC_LINEAR_PROGRAMMING_COLGEN_COLGEN_H
