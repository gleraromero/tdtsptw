//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_LINEAR_PROGRAMMING_MODEL_EXPRESSION_H
#define GOC_LINEAR_PROGRAMMING_MODEL_EXPRESSION_H

#include <iostream>
#include <unordered_map>
#include <vector>

#include "goc/linear_programming/model/valuation.h"
#include "goc/linear_programming/model/variable.h"
#include "goc/print/printable.h"

#define ESUM(iter_range, exp) ([&] () { goc::Expression fsum__exp; for (auto& iter_range) { fsum__exp += (exp); } return fsum__exp; })()

namespace goc
{
class Constraint;

// This class represents an expression which is a linear combination of variables multiplied by coefficients.
// Example: 2x + 3y + 5.
// Invariant: The expression is always normalized. E.g. 2x + 2x is not normalized, 4x is normalized.
// Invariant: The expression does not have 0 coefficients associated with variables.
class Expression : public Printable
{
public:
	// Creates an empty expression (represents the expression 0).
	Expression();
	
	// Creates an expression representing the scalar from the parameter.
	Expression(double scalar);
	
	// Creates an expression representing 1.0 * v.
	Expression(const Variable& v);
	
	// Adds the expression e to this expression.
	void operator+=(const Expression& e);
	
	// Subtracts the expression e to this expression.
	void operator-=(const Expression& e);
	
	// Subtracts the scalar from this expression.
	void operator-=(double scalar);
	
	// Adds the scalar from this expression.
	void operator+=(double scalar);
	
	// Multiplies the scalar by this expression.
	void operator*=(double scalar);
	
	// Divides this expression by the scalar.
	void operator/=(double scalar);
	
	// Adds 1.0 * v to this expression.
	void operator+=(const Variable& v);
	
	// Subtracts 1.0 * v to this expression.
	void operator-=(const Variable& v);
	
	// Returns: a new expression *this + e.
	Expression operator+(const Expression& e) const;
	
	// Returns: a new expression *this - e.
	Expression operator-(const Expression& e) const;
	
	// Returns: a new expression *this + scalar.
	Expression operator+(double scalar) const;
	
	// Returns: a new expression *this - scalar.
	Expression operator-(double scalar) const;
	
	// Returns: a new expression *this * scalar.
	Expression operator*(double scalar) const;
	
	// Returns: a new expression *this / scalar.
	Expression operator/(double scalar) const;
	
	// Returns: a constraint (this) <= e.
	Constraint LEQ(const Expression& e) const;
	
	// Returns: a constraint (this) >= e.
	Constraint GEQ(const Expression& e) const;
	
	// Returns: a constraint (this) == e.
	Constraint EQ(const Expression& e) const;
	
	// Sets the coefficient of the variable to 'coefficient'.
	// Observation: if coefficient is 0, the variable is removed from the expression.
	void SetVariableCoefficient(const Variable& variable, double coefficient);
	
	// Sets the scalar value to scalar.
	void SetScalar(double scalar);
	
	// Returns: the coefficient associated with the variable.
	// Observation: if the variable is not present a coefficient of 0 is returned.
	double VariableCoefficient(const Variable& variable) const;
	
	// Returns: the scalar value.
	double Scalar() const;
	
	// Returns: number of variable terms with non zero coefficient.
	int NonZeroVariableTermCount() const;
	
	// Returns: a sequence with the terms involving variables paired with their coefficients (!= 0).
	std::vector<std::pair<Variable, double>> Terms() const;
	
	// Returns: the expression evaluated on the valuation v.
	double Value(const Valuation& v) const;
	
	// Prints the expression to the stream.
	// Format: c1x1 + c2x2 + ... + cnxn + scalar. (+ scalar if scalar != 0 or n==0).
	virtual void Print(std::ostream& os) const;
	
private:
	std::unordered_map<Variable, double> coefficients_; // Stores the coefficients of each variable.
	double scalar_; // Stores the scalar of the expression.
};

// Returns: a new expression representing scalar * v.
Expression operator*(double scalar, const Variable& v);

// Returns: a new expression representing 1.0 * v + e.
Expression operator+(const Variable& v, const Expression& e);

// Returns: a new expression representing 1.0 * v - e.
Expression operator-(const Variable& v, const Expression& e);

// Returns: a new expression representing 1.0 * v1 + 1.0 * v2.
Expression operator+(const Variable& v1, const Variable& v2);

// Returns: a new expression representing 1.0 * v1 - 1.0 * v2.
Expression operator-(const Variable& v1, const Variable& v2);

// Returns: a new expression scalar - v.
Expression operator-(double scalar, const Variable& v);

// Returns: a new expression representing scalar * e.
Expression operator*(double scalar, const Expression& e);
} // namespace goc

#endif //GOC_LINEAR_PROGRAMMING_MODEL_EXPRESSION_H
