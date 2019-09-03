//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_MATH_LINEAR_FUNCTION_H
#define GOC_MATH_LINEAR_FUNCTION_H

#include <iostream>

#include "goc/math/interval.h"
#include "goc/math/point_2d.h"
#include "goc/lib/json.hpp"
#include "goc/print/printable.h"

namespace goc
{
// This class represents a segment or bounded linear function.
// An instance can be represented by two points (p1, p2) and the function
// connecting those points.
// Invariant: p1.x <= p2.x.
class LinearFunction : public Printable
{
public:
	Interval domain, image;
	double slope, intercept;
	
	LinearFunction() = default;
	
	// Creates a linear function from p1 to p2.
	// Precondition: p1.x <= p2.x.
	LinearFunction(const Point2D& p1, const Point2D& p2);
	
	// Returns: the value of the function f(x).
	// Precondition: x \in domain.
	double Value(double x) const;
	
	// Returns: the value of the function f(x).
	// Precondition: x \in domain.
	double operator()(double x) const;
	
	// Returns: the last x such that f(x) = y.
	// Observation: if slope == 0 the function is not inversible therefore many f(x) = y might exist.
	// Precondition: y \in image.
	double PreValue(double y) const;
	
	// Returns: if this function intersects the function f inside both functions domains.
	bool Intersects(const LinearFunction& f) const;
	
	// Returns: the coordinate x such that this(x) == f(x).
	// Precondition: Intersects(*this, f).
	double Intersection(const LinearFunction& f) const;
	
	// Returns: the inverse function of this function (f) if exists, otherwise returns g(y) = max{x : f(x)=y }.
	LinearFunction Inverse() const;
	
	// Restricts the domain.
	// Returns: the restricted function.
	LinearFunction RestrictDomain(const Interval& domain) const;
	
	// Restricts the image.
	// Returns: the restricted function.
	LinearFunction RestrictImage(const Interval& image) const;
	
	// Prints the function.
	// Format: {(x1, y1)->(x2,y2)}.
	void Print(std::ostream& os) const;
	
	// Returns: if two functions have the same p1 and p2.
	bool operator==(const LinearFunction& f) const;
	
	// Returns: if two functions have different p1 and p2.
	bool operator!=(const LinearFunction& f) const;
};

// JSON format: [p1, p2].
void from_json(const nlohmann::json& j, LinearFunction& f);

void to_json(nlohmann::json& j, const LinearFunction& f);

// Returns: f.domain
inline Interval dom(const LinearFunction& f) { return f.domain; }

// Returns: f.image
inline Interval img(const LinearFunction& f) { return f.image; }

// Returns: h(x) = f(x)+g(x).
// Obs: dom(h) = dom(f) \cap dom(g).
LinearFunction operator+(const LinearFunction& f, const LinearFunction& g);

// Returns: h(x) = f(x)*g(x).
// Obs: dom(h) = dom(f) \cap dom(g).
LinearFunction operator*(const LinearFunction& f, const LinearFunction& g);

// Returns: h(x) = f(x)+a.
LinearFunction operator+(const LinearFunction& f, double a);

// Returns: h(x) = f(x)*a.
LinearFunction operator*(const LinearFunction& f, double a);
} // namespace goc
#endif //GOC_MATH_LINEAR_FUNCTION_H
