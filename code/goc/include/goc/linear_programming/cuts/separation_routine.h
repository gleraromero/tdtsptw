//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_LINEAR_PROGRAMMING_CUTS_SEPARATION_ROUTINE_H
#define GOC_LINEAR_PROGRAMMING_CUTS_SEPARATION_ROUTINE_H

#include <list>

#include "goc/linear_programming/model/constraint.h"
#include "goc/linear_programming/model/valuation.h"

namespace goc
{
// This is a base class for all separation routines (cuts and lazy constraints).
// All implemented routines should implement this interface in order to be used in a SeparationAlgorithm.
class SeparationRoutine
{
public:
	// Separates the current 'solution' by generating up to 'count_limit' violated cuts.
	// solution: solution to be separated from the model.
	// node_number: number of node in the BB tree (0 is root).
	// count_limit: soft limit of the number of cuts that should be generated. The routine can return more and
	//				the separation algorithm will then decide which to add.
	// node_bound: bound of the current node.
	// Returns: a sequence of cuts to add to the formulation.
	virtual std::vector<Constraint> Separate(const Valuation& solution, int node_number, int count_limit,
		double node_bound) const = 0;
};
} // namespace goc

#endif //GOC_LINEAR_PROGRAMMING_CUTS_SEPARATION_ROUTINE_H
