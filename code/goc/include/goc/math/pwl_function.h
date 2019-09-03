//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_MATH_PWL_FUNCTION_H
#define GOC_MATH_PWL_FUNCTION_H

#include <iostream>
#include <vector>

#include "goc/lib/json.hpp"
#include "goc/math/interval.h"
#include "goc/math/linear_function.h"
#include "goc/print/printable.h"

namespace goc
{
// This class represents a piecewise linear function. It has a sequence of linear functions with bounded domains.
// Invariant: the linear functions are non overlapping and are increasing in domain.
// Invariant: the function is stored normalized. A function is normalized iif no two consecutive pieces have the same
// 			  slope, intercept, and share the end and beginning of their domains.
// Example: [p1={(1,2),(2,3)},p2={(2,3),(3,4)}] is not normalized. [p1={(1,2),(3,4)}] is normalized.
class PWLFunction : public Printable
{
public:
	// Returns: f(x)=a with the specific domain.
	static PWLFunction ConstantFunction(double a, Interval domain);
	
	// Returns: f(x)=x with the specific domain.
	static PWLFunction IdentityFunction(Interval domain);
	
	// Creates an empty piecewise linear function.
	PWLFunction();
	
	// Creates a piecewise linear function with the specified pieces.
	PWLFunction(const std::vector<LinearFunction>& pieces);
	
	// Adds the piece at the end of the function.
	// Keeps the normalization invariant automatically.
	void AddPiece(const LinearFunction& piece);
	
	// Removes the last piece from the function.
	// Precondition: PieceCount() > 0.
	void PopPiece();
	
	// Returns: if the function has no pieces.
	bool Empty() const;
	
	// Returns: the numer of pieces of the function.
	int PieceCount() const;
	
	// Returns: a vector with the function pieces.
	const std::vector<LinearFunction>& Pieces() const;
	
	// Returns: the i-th piece of the function.
	// Precondition: i > PieceCount().
	const LinearFunction& Piece(int i) const;
	
	// Returns: the i-th piece of the function.
	// Precondition: i > PieceCount().
	const LinearFunction& operator[](int i) const;
	
	// Returns: the first piece of the function.
	// Precondition: !Empty().
	const LinearFunction& FirstPiece() const;
	
	// Returns: the last piece of the function.
	// Precondition: !Empty().
	const LinearFunction& LastPiece() const;
	
	// Returns: the last piece index that includes x in its domain.
	// Precondition: x \in dom(p) for any piece p.
	int PieceIncluding(double x) const;
	
	// Returns: the smallest interval [m, M] that includes all pieces domains.
	// Observation: if Empty() then returns [INFTY, -INFTY].
	Interval Domain() const;
	
	// Returns: the smallest interval [m, M] that includes all pieces images.
	// Observation: if Empty() then returns [INFTY, -INFTY].
	Interval Image() const;
	
	// Returns: the evaluation of the piece that includes x in its domain.
	// Exception: if no piece includes x in its domain, it throws an exception.
	double Value(double x) const;
	
	// Returns: the evaluation of the piece that includes x in its domain.
	// Exception: if no piece includes x in its domain, it throws an exception.
	double operator()(double x) const;
	
	// Returns: the last x such that f(x) = y. Notice that if the function is not bijective it may contain multiple
	// x such that f(x) = y.
	// Exception: if no f(x) = y, then it throws an exception.
	double PreValue(double y) const;
	
	// Returns: the composition of this function (f) and g, i.e. fog(x) == f(g(x)).
	// Observation: the domain of the new function are those x such that g(x) \in dom(f).
	PWLFunction Compose(const PWLFunction& g) const;
	
	// Returns: the inverse of this function (f) if is inversible, otherwise returns g(y) = max{x : f(x) = y}.
	PWLFunction Inverse() const;
	
	// Restricts the domain to only the pieces included in the parameter.
	// Returns: the restricted function.
	PWLFunction RestrictDomain(const Interval& domain) const;
	
	// Restricts the image to only the pieces included in the parameter.
	// Returns: the restricted function.
	PWLFunction RestrictImage(const Interval& image) const;
	
	// Prints the function.
	// Format: [p1, p2, ..., pn].
	virtual void Print(std::ostream& os) const;
	
	// Returns: if all the pieces of both functions are the same.
	bool operator==(const PWLFunction& f) const;
	
	// Returns: if any piece of both functions is different.
	bool operator!=(const PWLFunction& f) const;
private:
	// Updates the image_ attribute to keep it updated after a Pop() operation.
	void UpdateImage();
	
	std::vector<LinearFunction> pieces_;
	Interval domain_, image_;
};

// JSON format: [p1, p2, ..., pn].
void from_json(const nlohmann::json& j, PWLFunction& f);

void to_json(nlohmann::json& j, const PWLFunction& f);

// Returns: the function h(x) = f(x)+g(x).
// Observation: only returns h(x) for x \in dom(f) \cap dom(g).
PWLFunction operator+(const PWLFunction& f, const PWLFunction& g);

// Returns: the function h(x) = f(x)-g(x).
// Observation: only returns h(x) for x \in dom(f) \cap dom(g).
PWLFunction operator-(const PWLFunction& f, const PWLFunction& g);

// Returns: the function h(x) = f(x)*g(x).
// Observation: only returns h(x) for x \in dom(f) \cap dom(g).
PWLFunction operator*(const PWLFunction& f, const PWLFunction& g);

// Returns: the function h(x) = f(x)+a.
PWLFunction operator+(const PWLFunction& f, double a);
PWLFunction operator+(double a, const PWLFunction& f);

// Returns: the function h(x) = f(x)-a.
PWLFunction operator-(const PWLFunction& f, double a);

// Returns: the function h(x) = a-f(x).
PWLFunction operator-(double a, const PWLFunction& f);

// Returns: the function h(x) = f(x)*a.
PWLFunction operator*(const PWLFunction& f, double a);
PWLFunction operator*(double a, const PWLFunction& f);

// Returns: h(x) = max(f(x), g(x)).
// Obs: If x \in dom(f), but x \not\in dom(g), then h(x) = f(x). Analogously, for the opposite case.
goc::PWLFunction Max(const goc::PWLFunction& f, const goc::PWLFunction& g);

// Returns: h(x) = max(f(x), a).
goc::PWLFunction Max(const goc::PWLFunction& f, double a);
goc::PWLFunction Max(double a, const goc::PWLFunction& f);

// Returns: h(x) = min(f(x), g(x)).
// Obs: If x \in dom(f), but x \not\in dom(g), then h(x) = f(x). Analogously, for the opposite case.
goc::PWLFunction Min(const goc::PWLFunction& f, const goc::PWLFunction& g);

// Returns: h(x) = min(f(x), a).
goc::PWLFunction Min(const goc::PWLFunction& f, double a);
goc::PWLFunction Min(double a, const goc::PWLFunction& f);

// Returns: f.domain
inline Interval dom(const PWLFunction& f) { return f.Domain(); }

// Returns: f.image
inline Interval img(const PWLFunction& f) { return f.Image(); }
} // namespace goc

#endif //GOC_MATH_PWL_FUNCTION_H
