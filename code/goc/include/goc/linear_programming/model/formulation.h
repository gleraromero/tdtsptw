//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_LINEAR_PROGRAMMING_MODEL_FORMULATION_H
#define GOC_LINEAR_PROGRAMMING_MODEL_FORMULATION_H

#include <iostream>
#include <string>
#include <vector>

#include "goc/linear_programming/model/variable.h"
#include "goc/linear_programming/model/valuation.h"
#include "goc/linear_programming/model/constraint.h"
#include "goc/linear_programming/model/expression.h"
#include "goc/linear_programming/cuts/separation_routine.h"
#include "goc/math/number_utils.h"
#include "goc/print/printable.h"

namespace goc
{
// Represents the domain of the variables in a (mixed integer) linear programming model.
enum class VariableDomain { Real, Integer, Binary };

// Represents a (mixed integer) linear programming model.
// It is a protocol designed to abstract specific implementations for solvers (CPLEX, Gurobi, etc).
class Formulation : public Printable
{
public:
	enum ObjectiveSense { Minimization, Maximization };
	
	virtual ~Formulation() = default;
	
	// Adds the constraint to the formulation.
	// Observation: Constraints are normalized before added to the formulation, this is, the coefficients of each
	// variables are added up on the left side, and the constants are added up on the right side.
	// Returns: the index of the constraint added.
	virtual int AddConstraint(const Constraint& constraint) = 0;
	
	// Removes the constraint with index 'constraint_index' from the formulation.
	// If such constraint does not exists, it does nothing.
	virtual void RemoveConstraint(int constraint_index) = 0;
	
	// Adds a lazy constraint to the model. When a solution is found by the solver, it tries to add some lazy
	// constraint that separates that solution from the actual solution space.
	// Observation: if lazy_constraint == nullptr, the call to the function is ignored.
	virtual void AddLazyConstraint(SeparationRoutine* lazy_constraint) = 0;
	
	// Removes a lazy constraint from the model.
	// Observation: if lazy_constraint == nullptr, the call to the function is ignored.
	virtual void RemoveLazyConstraint(SeparationRoutine* lazy_constraint) = 0;
	
	// Adds a variable to the model.
	// Returns: a Variable with its correspondent name and index in the formulation.
	virtual Variable AddVariable(const std::string& name, VariableDomain domain=VariableDomain::Real, double lower_bound=-INFTY, double upper_bound=INFTY) = 0;
	
	// Removes the variable with the same index than the one given as a parameter.
	// Observation: The name is ignored.
	virtual void RemoveVariable(const Variable& variable) = 0;
	
	// Sets the domain of the variable in the model. Domain might be Real, Integer, or Binary.
	virtual void SetVariableDomain(const Variable& variable, VariableDomain domain) = 0;
	
	// Sets the lower and upper bounds of the variable in the formulation.
	// Observation: use -INFTY or INFTY constants to specify no bounds.
	virtual void SetVariableBound(const Variable& v, double lower_bound, double upper_bound) = 0;
	
	// Sets the lower bound of the variable in the formulation.
	// Observation: use -INFTY or INFTY constants to specify no bounds.
	virtual void SetVariableLowerBound(const Variable& v, double lower_bound) = 0;
	
	// Sets the upper bound of the variable in the formulation.
	// Observation: use -INFTY or INFTY constants to specify no bounds.
	virtual void SetVariableUpperBound(const Variable& v, double upper_bound) = 0;
	
	// Sets the objective function as a minimization function.
	virtual void Minimize(const Expression& objective_function) = 0;
	
	// Sets the objective function as a maximization function.
	virtual void Maximize(const Expression& objective_function) = 0;
	
	// Gets the constraint right hand constant.
	virtual void SetConstraintRightHandSide(int constraint_index, double value) = 0;
	
	// Sets the constraint coefficient of a certain variable.
	virtual void SetConstraintCoefficient(int constraint_index, const Variable& variable, double coefficient) = 0;
	
	// Sets the coefficient of a variable in the objective function.
	virtual void SetObjectiveCoefficient(const Variable& variable, double coefficient) = 0;
	
	// Returns: the objective function sense.
	virtual ObjectiveSense GetObjectiveSense() const = 0;
	
	// Returns: the coefficient of a variable in the objective function.
	virtual double GetObjectiveCoefficient(const Variable& variable) const = 0;
	
	// Returns: the right hand constant of constraint with ID 'constraint_id'.
	virtual double GetConstraintRightHandSide(int constraint_index) const = 0;
	
	// Returns: the coefficient of a variable in the constraint with ID 'constraint_id'.
	virtual double GetConstraintCoefficient(int constraint_index, const Variable& variable) = 0;
	
	// Returns: the domain of the variable with the specified index.
	virtual VariableDomain GetVariableDomain(const Variable& variable) const = 0;
	
	// Returns: a pair (lower, upper) with the bounds of the variables.
	// Observation: -INFTY and INFTY specifies no bounds are contemplated.
	virtual std::pair<double, double> GetVariableBound(const Variable& variable) const = 0;
	
	// Returns: the objective function expression.
	virtual Expression ObjectiveFunction() const = 0;
	
	// Returns: a sequence with all the variables.
	virtual std::vector<Variable> Variables() const = 0;
	
	// Returns: a sequence with all the constraints.
	virtual std::vector<Constraint> Constraints() const = 0;
	
	// Returns: a sequence with all the lazy constraints.
	virtual const std::vector<SeparationRoutine*>& LazyConstraints() const = 0;
	
	// Returns: the number of variables in the model.
	virtual int VariableCount() const = 0;
	
	// Returns: the number of constraints in the model.
	virtual int ConstraintCount() const = 0;
	
	// Returns: the variable at a specified index.
	virtual Variable VariableAtIndex(int variable_index) const = 0;
	
	// Returns: the objective function evaluated at valuation 'v'.
	virtual double EvaluateValuation(const Valuation& v) const = 0;
	
	// Returns: if the constraints of the model are satisfied by 'v'.
	// 	- verbose: indicates if a contraint violated should be showed.
	virtual bool IsFeasibleValuation(const Valuation& v, bool verbose=false) const = 0;
	
	// Returns: a copy in the heap of the current formulation.
	// Observation: memory should be managed by the receiver and the pointer should be freed.
	virtual Formulation* Copy() const = 0;
	
	// Prints the formulation.
	virtual void Print(std::ostream& os) const = 0;
};
} // namespace goc

#endif //GOC_LINEAR_PROGRAMMING_MODEL_FORMULATION_H
