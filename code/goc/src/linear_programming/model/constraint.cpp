//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/linear_programming/model/constraint.h"

#include "goc/math/number_utils.h"
#include "goc/exception/exception_utils.h"

using namespace std;

namespace goc
{
const Expression& Constraint::LeftSide() const
{
	return left_side_;
}

double Constraint::RightSide() const
{
	return right_side_;
}

enum Constraint::Sense Constraint::Sense() const
{
	return sense_;
}

bool Constraint::Holds(const Valuation& v) const
{
	double left_value = left_side_.Value(v);
	double right_value = right_side_;
	switch (Sense())
	{
		case LessEqual: return epsilon_smaller_equal(left_value, right_value);
		case GreaterEqual: return epsilon_bigger_equal(left_value, right_value);
		case Equality: return epsilon_equal(left_value, right_value);
	}
	fail("Unrecognized constraint sense.");
	return false;
}

void Constraint::Print(ostream& os) const
{
	os << left_side_;
	switch (Sense())
	{
		case LessEqual: os << " ≤ "; break;
		case GreaterEqual: os << " ≥ "; break;
		case Equality: os << " = "; break;
	}
	os << right_side_;
}

Constraint::Constraint(const Expression& left, const Expression& right, enum Sense sense)
{
	left_side_ = left - right;
	right_side_ = -left_side_.Scalar();
	left_side_.SetScalar(0.0);
	sense_ = sense;
}
} // namespace goc