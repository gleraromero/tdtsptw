//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_LINEAR_PROGRAMMING_MODEL_VALUATION_H
#define GOC_LINEAR_PROGRAMMING_MODEL_VALUATION_H

#include <iostream>
#include <unordered_map>

#include "goc/lib/json.hpp"
#include "goc/linear_programming/model/variable.h"
#include "goc/print/printable.h"

namespace goc
{
// A valuation is an assignment of values to variables.
// - Invariant: all values stored must not be zero.
// - Observation: if a variable is not present in the valuation its value is 0.
class Valuation : public Printable
{
public:
	// Creates a valuation with all values equal to 0.
	Valuation() = default;
	
	// Sets the value of the varaible v.
	void SetValue(const Variable& v, double value);
	
	// Returns: the value of the variable v.
	// Observation: if the variable v is not in the inner structure it returns 0.0.
	double operator[](const Variable& v) const;
	
	// Returns: the value of the variable v.
	// Observation: if the variable v is not in the inner structure it returns 0.0.
	double at(const Variable& v) const;
	
	// Returns: if all the values associated with variables are integer.
	bool IsInteger() const;
	
	// Prints the valuation.
	// Format: the output is a json { "var_name1": value1, ..., "var_namen": valuen }.
	virtual void Print(std::ostream& os) const;
	
	// Iterators to traverse the non-zero valued variables in the valuation as pairs (variable, value).
	std::unordered_map<Variable, double>::const_iterator begin() const;
	std::unordered_map<Variable, double>::const_iterator end() const;
	
private:
	std::unordered_map<Variable, double> values_;
};

void to_json(nlohmann::json& j, const Valuation& v);
} // namespace goc

#endif //GOC_LINEAR_PROGRAMMING_MODEL_VALUATION_H
