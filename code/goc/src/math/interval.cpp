//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/math/interval.h"

#include "goc/math/number_utils.h"

using namespace std;
using namespace nlohmann;

namespace goc
{
Interval::Interval()
{
	left = INFTY;
	right = -INFTY;
}

Interval::Interval(double left, double right) : left(left), right(right)
{
}

bool Interval::IsIncludedIn(const Interval& r) const
{
	return epsilon_bigger_equal(left, r.left) && epsilon_smaller_equal(right, r.right);
}

bool Interval::Intersects(const Interval& r) const
{
	return !(epsilon_smaller(r.right, left) || epsilon_bigger(r.left, right));
}

void Interval::Print(std::ostream& os) const
{
	os << "[" << left << ", " << right << "]";
}

bool Interval::operator==(const Interval& i) const
{
	if (Empty()) return i.Empty();
	return epsilon_equal(left, i.left) && epsilon_equal(right, i.right);
}

bool Interval::operator!=(const Interval& i) const
{
	return !(*this == i);
}

void from_json(const json& j, Interval& i)
{
	i.left = j[0];
	i.right = j[1];
}

void to_json(json& j, const Interval& i)
{
	j = vector<json>();
	j.push_back(i.left);
	j.push_back(i.right);
}
} // namespace goc