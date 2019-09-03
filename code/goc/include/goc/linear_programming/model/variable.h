//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_LINEAR_PROGRAMMING_MODEL_VARIABLE_H
#define GOC_LINEAR_PROGRAMMING_MODEL_VARIABLE_H

#include <iostream>
#include <string>

#include "goc/print/printable.h"

namespace goc
{
// This class represents a variable inside a formulation.
// We can identify a variable by the index of the column in the formulation.
class Variable : public Printable
{
public:
	std::string name; // name of the variable.
	
	// Default constructor.
	Variable();
	
	// Returns: the index of the variable in the model (i.e. column index).
	// Observation: It can change if some variables are deleted.
	int Index() const;
	
	// Prints the name of the variable.
	virtual void Print(std::ostream& os) const;
	
	// Compares the index of the variables.
	bool operator<(const Variable& v) const;
	
	// Compares the index of the variables.
	bool operator==(const Variable& v) const;
	
	// Compares the index of the variables.
	bool operator!=(const Variable& v) const;
	
private:
	friend class CplexFormulation;
	friend class std::hash<Variable>;
	
	Variable(const std::string& name, int* index);
	int* index_; // a pointer to an integer that keeps the index of this variable updates with respect to the model.
};
} // namespace goc

// Implement hash function in order to use Variable as keys in unordered_map and unordered_set.
namespace std
{
template<>
class hash<goc::Variable> {
public:
	size_t operator()(const goc::Variable& v) const
	{
		return std::hash<int*>()(v.index_);
	}
};
} // namespace std

#endif //GOC_LINEAR_PROGRAMMING_MODEL_VARIABLE_H
