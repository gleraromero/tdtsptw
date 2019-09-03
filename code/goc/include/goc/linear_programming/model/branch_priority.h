//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_LINEAR_PROGRAMMING_MODEL_BRANCH_PRIORITY_H
#define GOC_LINEAR_PROGRAMMING_MODEL_BRANCH_PRIORITY_H

#include "goc/linear_programming/model/variable.h"

#include "goc/linear_programming/model/variable.h"

namespace goc
{
// In a branch and bound algorithm, there is a phase called branching where typically a variable is selected
// to be branched upon. Generally, solvers let specify priorities to the variables so that one variable is prefered
// over others in this selection process and the tree is smaller.
// This struct represents such priority and is used on the BCSolver.
struct BranchPriority
{
	enum BranchDirection { Up=0, Down=1, Any=2};
	
	Variable variable;
	int priority;
	BranchDirection direction;
	
	BranchPriority() = default;
	
	BranchPriority(const Variable& variable, int priority, BranchDirection direction)
		: variable(variable), priority(priority), direction(direction)
	{}
};
} // namespace goc

#endif //GOC_LINEAR_PROGRAMMING_MODEL_BRANCH_PRIORITY_H
