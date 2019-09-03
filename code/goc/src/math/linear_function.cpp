//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/math/linear_function.h"

#include <vector>

#include "goc/exception/exception_utils.h"
#include "goc/math/number_utils.h"
#include "goc/string/string_utils.h"

using namespace std;
using namespace nlohmann;

namespace goc
{
LinearFunction::LinearFunction(const Point2D& p1, const Point2D& p2)
	: domain(p1.x, p2.x), image(min(p1.y, p2.y), max(p1.y, p2.y))
{
	slope = epsilon_equal(p2.x, p1.x) ? 0.0 : (p2.y - p1.y) / (p2.x - p1.x);
	if (epsilon_equal(slope, 0.0)) slope = 0.0;
	intercept = p2.y - slope * p2.x;
}

double LinearFunction::Value(double x) const
{
	return slope * x + intercept;
}

double LinearFunction::operator()(double x) const
{
	return Value(x);
}

double LinearFunction::PreValue(double y) const
{
	if (!image.Includes(y)) fail(STR(y) + " is not in the image " + STR(image));
	if (epsilon_equal(slope, 0.0)) return domain.right;
	return (y - intercept) / slope;
}

bool LinearFunction::Intersects(const LinearFunction& f) const
{
	// If both pieces have the same slope, check if they have the same intercept.
	if (epsilon_equal(slope, f.slope)) return epsilon_equal(intercept, f.intercept);
	// Otherwise check if intersection is inside both functions domains.
	double intersection = Intersection(f);
	return domain.Includes(intersection) && f.domain.Includes(intersection);
}

double LinearFunction::Intersection(const LinearFunction& l) const
{
	return (l.intercept - intercept) / (slope - l.slope);
}

LinearFunction LinearFunction::Inverse() const
{
	return LinearFunction({min(image), PreValue(min(image))}, {max(image), PreValue(max(image))});
}

LinearFunction LinearFunction::RestrictDomain(const Interval& domain) const
{
	double left = max(this->domain.left, domain.left);
	double right = min(this->domain.right, domain.right);
	return LinearFunction({left, Value(left)}, {right, Value(right)});
}

LinearFunction LinearFunction::RestrictImage(const Interval& image) const
{
	if (epsilon_equal(slope, 0.0))
	{
		if (image.Includes(intercept)) return *this;
		else fail("Linear function is empty");
	}
	double left = max(this->image.left, image.left);
	double right = min(this->image.right, image.right);
	double pre_left = PreValue(left), pre_right = PreValue(right);
	if (pre_left > pre_right) swap(pre_left, pre_right);
	return LinearFunction({pre_left, left}, {pre_right, right});
}

void LinearFunction::Print(ostream& os) const
{
	os << "{" << Point2D(domain.left, Value(domain.left)) << "->" << Point2D(domain.right, Value(domain.right)) << "}";
}

bool LinearFunction::operator==(const LinearFunction& f) const
{
	return domain == f.domain && image == f.image && epsilon_equal(slope, f.slope) && epsilon_equal(intercept, f.intercept);
}

bool LinearFunction::operator!=(const LinearFunction& f) const
{
	return !(*this == f);
}

void from_json(const json& j, LinearFunction& f)
{
	f = LinearFunction(j[0], j[1]);
}

void to_json(json& j, const LinearFunction& f)
{
	j = vector<Point2D>();
	j.push_back(Point2D(f.domain.left, f.Value(f.domain.left)));
	j.push_back(Point2D(f.domain.right, f.Value(f.domain.right)));
}

LinearFunction operator+(const LinearFunction& f, const LinearFunction& g)
{
	if (!f.domain.Intersects(g.domain)) return LinearFunction(Point2D(0.0, 0.0), Point2D(-1.0, 0.0));
	double l = max(f.domain.left, g.domain.left);
	double r = min(f.domain.right, g.domain.right);
	return LinearFunction(Point2D(l, f.Value(l)+g.Value(l)), Point2D(r, f.Value(r)+g.Value(r)));
}

// Returns: h(x) = f(x)*g(x).
// Precondition: dom(f) == dom(g).
LinearFunction operator*(const LinearFunction& f, const LinearFunction& g)
{
	if (!f.domain.Intersects(g.domain)) return LinearFunction(Point2D(0.0, 0.0), Point2D(-1.0, 0.0));
	double l = max(f.domain.left, g.domain.left);
	double r = min(f.domain.right, g.domain.right);
	return LinearFunction(Point2D(l, f.Value(l)*g.Value(l)), Point2D(r, f.Value(r)*g.Value(r)));
}

// Returns: h(x) = f(x)+a.
LinearFunction operator+(const LinearFunction& f, double a)
{
	return LinearFunction(Point2D(f.domain.left, f.Value(f.domain.left)+a), Point2D(f.domain.right, f.Value(f.domain.right)+a));
}

// Returns: h(x) = f(x)*a.
LinearFunction operator*(const LinearFunction& f, double a)
{
	return LinearFunction(Point2D(f.domain.left, f.Value(f.domain.left)*a), Point2D(f.domain.right, f.Value(f.domain.right)*a));
}
} // namespace goc