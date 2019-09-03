//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/linear_programming/model/expression.h"

#include <math.h>

#include "goc/linear_programming/model/valuation.h"
#include "goc/linear_programming/model/constraint.h"
#include "goc/collection/collection_utils.h"
#include "goc/math/number_utils.h"
#include "goc/exception/exception_utils.h"

using namespace std;

namespace goc
{
Expression::Expression() : Expression(0.0)
{ }

Expression::Expression(double scalar) : scalar_(scalar)
{ }

Expression::Expression(const Variable& v) : Expression()
{
	SetVariableCoefficient(v, 1.0);
}

void Expression::operator+=(const Expression& e)
{
	scalar_ += e.scalar_;
	for (auto& variable_coefficient: e.coefficients_)
	{
		double new_coefficient = VariableCoefficient(variable_coefficient.first)+variable_coefficient.second;
		if (epsilon_equal(new_coefficient, 0.0)) coefficients_.erase(variable_coefficient.first);
		else coefficients_[variable_coefficient.first] = new_coefficient;
	}
}

void Expression::operator-=(const Expression& e)
{
	scalar_ -= e.scalar_;
	for (auto& variable_coefficient: e.coefficients_)
	{
		double new_coefficient = VariableCoefficient(variable_coefficient.first)-variable_coefficient.second;
		if (epsilon_equal(new_coefficient, 0.0)) coefficients_.erase(variable_coefficient.first);
		else coefficients_[variable_coefficient.first] = new_coefficient;
	}
}

void Expression::operator-=(double scalar)
{
	scalar_ -= scalar;
}

void Expression::operator+=(double scalar)
{
	scalar_ += scalar;
}

void Expression::operator*=(double scalar)
{
	scalar_ *= scalar;
	if (epsilon_equal(scalar, 0.0)) { coefficients_.clear(); return; }
	for (auto& variable_coefficient: coefficients_) variable_coefficient.second *= scalar;
}

void Expression::operator/=(double scalar)
{
	if (epsilon_equal(scalar, 0.0)) fail("Expression: Division by 0.");
	scalar_ /= scalar;
	for (auto& variable_coefficient: coefficients_) variable_coefficient.second /= scalar;
}

void Expression::operator+=(const Variable& v)
{
	*this += 1.0*v;
}

void Expression::operator-=(const Variable& v)
{
	*this -= 1.0*v;
}

Expression Expression::operator+(const Expression& e) const
{
	Expression new_expression = *this;
	new_expression += e;
	return new_expression;
}

Expression Expression::operator-(const Expression& e) const
{
	Expression new_expression = *this;
	new_expression -= e;
	return new_expression;
}

Expression Expression::operator+(double scalar) const
{
	Expression new_expression = *this;
	new_expression += scalar;
	return new_expression;
}

Expression Expression::operator-(double scalar) const
{
	Expression new_expression = *this;
	new_expression -= scalar;
	return new_expression;
}

Expression Expression::operator*(double scalar) const
{
	Expression new_expression = *this;
	new_expression *= scalar;
	return new_expression;
}

Expression Expression::operator/(double scalar) const
{
	Expression new_expression = *this;
	new_expression /= scalar;
	return new_expression;
}

Constraint Expression::LEQ(const Expression& e) const
{
	return Constraint(*this, e, Constraint::LessEqual);
}

Constraint Expression::GEQ(const Expression& e) const
{
	return Constraint(*this, e, Constraint::GreaterEqual);
}

Constraint Expression::EQ(const Expression& e) const
{
	return Constraint(*this, e, Constraint::Equality);
}

void Expression::SetVariableCoefficient(const Variable& variable, double coefficient)
{
	if (epsilon_equal(coefficient, 0.0)) coefficients_.erase(variable);
	else coefficients_[variable] = coefficient;
}

// Sets the scalar value to scalar.
void Expression::SetScalar(double scalar)
{
	scalar_ = scalar;
}

// Returns: the coefficient associated with the variable.
// Observation: if the variable is not present a coefficient of 0 is returned.
double Expression::VariableCoefficient(const Variable& variable) const
{
	return includes_key(coefficients_, variable) ? coefficients_.at(variable) : 0.0;
}

// Returns: the scalar value.
double Expression::Scalar() const
{
	return scalar_;
}

int Expression::NonZeroVariableTermCount() const
{
	return coefficients_.size();
}

vector<pair<Variable, double>> Expression::Terms() const
{
	vector<pair<Variable, double>> terms;
	for (auto& variable_coefficient: coefficients_) terms.push_back(variable_coefficient);
	return terms;
}

// Returns: the expression evaluated on the valuation v.
double Expression::Value(const Valuation& v) const
{
	double value = scalar_;
	for (auto& variable_coefficient: coefficients_) value += v[variable_coefficient.first] * variable_coefficient.second;
	return value;
}

void Expression::Print(ostream& os) const
{
	if (coefficients_.empty()) { os << scalar_; return; }
	bool is_first = true;
	for (auto& variable_coefficient: coefficients_)
	{
		if (!is_first) os << " + ";
		is_first = false;
		if (epsilon_equal(variable_coefficient.second, 1.0)) os << variable_coefficient.first;
		else os << variable_coefficient.second << " " << variable_coefficient.first;
	}
	if (epsilon_different(scalar_, 0.0)) os << " + " << scalar_;
}

Expression operator*(double scalar, const Variable& v)
{
	auto e = Expression();
	e.SetVariableCoefficient(v, scalar);
	return e;
}

Expression operator+(const Variable& v, const Expression& e)
{
	return e + v;
}

Expression operator-(const Variable& v, const Expression& e)
{
	return v + e * -1;
}

Expression operator+(const Variable& v1, const Variable& v2)
{
	return 1.0*v1+1.0*v2;
}

Expression operator-(const Variable& v1, const Variable& v2)
{
	return 1.0*v1-1.0*v2;
}

Expression operator-(double scalar, const Variable& v)
{
	return -1 * v + scalar;
}

Expression operator*(double scalar, const Expression& e)
{
	return e * scalar;
}
} // namespace goc