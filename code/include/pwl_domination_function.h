//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef TDTSPTW_PWL_DOMINATION_FUNCTION_H
#define TDTSPTW_PWL_DOMINATION_FUNCTION_H

#include <vector>

#include "goc/goc.h"

namespace tdtsptw
{
// This class represents a piecewise linear function with domination functionality.
// Basically, domination between PWL functions depends on which one is smaller than the other.
// We say that f1(t) is dominated by f2(t) in t \in dom(f1) if f1(t) >= f2(t).
// Invariant:
// 		- The domain of the pieces is disjunt and ascending.
class PWLDominationFunction
{
public:
	// Creates a PWLDominationFunction with f's pieces.
	PWLDominationFunction(const goc::PWLFunction& f);
	
	// Converts this function to a PWLFunction.
	explicit operator goc::PWLFunction() const;
	
	// @return Returns: The smallest interval [min_x, max_x] such that dom(f) \subseteq [min_x, max_x].
	goc::Interval Domain() const;
	
	// Returns: if no pieces are left.
	bool Empty() const;
	
	// Leaves in this function (f1) only x in dom(f1) | f1(x) < f2(x)+delta.
	// Returns: if this function (f1) is fully dominated (has no pieces left).
	bool DominatePieces(const goc::PWLFunction& f2, double delta = 0);
	
	// Returns: if this function (f1), forall x in dom(f1), f1(x) >= f2(x)+delta.
	// Precondition: f1, f2 are continuous.
	bool IsAlwaysDominated(const goc::PWLFunction& f2, double delta = 0);

private:
	// Removes piece i. (leaves domain_ and image_ inconsistent).
	// prev_i: index of the previous element of i.
	// Returns: next piece of i.
	int ErasePiece(int i, int prev_i);
	
	// Adds the piece after the piece i. (leaves domain_ and image_ inconsistent).
	// Returns: new piece index.
	int AddPieceAfter(int i, const goc::LinearFunction& piece);
	
	std::vector<goc::LinearFunction> pieces_;
	std::vector<int> next_;
	int first_, last_, size_;
	goc::Interval domain_; // minimum and maximum values of t where f(t) is in domain.
};
} // namespace tdtsptw

#endif //TDTSPTW_PWL_DOMINATION_FUNCTION_H
