//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_LINEAR_PROGRAMMING_MODEL_CONSTRAINT_H
#define GOC_LINEAR_PROGRAMMING_MODEL_CONSTRAINT_H

#include <iostream>

#include "goc/linear_programming/model/expression.h"
#include "goc/linear_programming/model/valuation.h"
#include "goc/print/printable.h"

namespace goc
{
// Represents a linear constraint.
// Invariant: 	it is semi-normalized, meaning that the left side is always with variables and scalar 0, and the right
//				side is always with variable coefficients 0 and the corresponding scalar.
//				Example: 2x + 5 <= 3y is not normalized. 2x - 3y <= -5 is normalized.
//				Constraints are created with class Expression LEQ, GEQ, and EQ methods.
class Constraint : public Printable
{
public:
	// Sense of the constraint (>=, <=, ==).
	enum Sense { GreaterEqual, LessEqual, Equality };
	
	// Returns: left side of the constraint (side with variables + coefficients).
	const Expression& LeftSide() const;
	
	// Returns: right side of the constriant (scalar).
	double RightSide() const;
	
	// Sense of the constraint.
	enum Sense Sense() const;
	
	// Returns: if the constraint holds for valuation v.
	bool Holds(const Valuation& v) const;
	
	// Prints the constraint.
	// Format: left (<=|>=|==) right.
	virtual void Print(std::ostream& os) const;
	
private:
	friend class Expression;
	
	// Creates a constraint left 'sense' right.
	// This constructor can only be called from Expression LEQ, GEQ, EQ methods.
	Constraint(const Expression& left, const Expression& right, enum Sense sense);
	
	Expression left_side_;
	double right_side_;
	enum Sense sense_;
};

} // namespace goc

#endif //GOC_LINEAR_PROGRAMMING_MODEL_CONSTRAINT_H
