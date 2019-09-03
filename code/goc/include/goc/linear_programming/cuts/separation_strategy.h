//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_LINEAR_PROGRAMMING_CUTS_SEPARATION_STRATEGY_H
#define GOC_LINEAR_PROGRAMMING_CUTS_SEPARATION_STRATEGY_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "goc/lib/json.hpp"
#include "goc/linear_programming/cuts/separation_routine.h"
#include "goc/print/printable.h"

namespace goc
{
// This class represents the strategies for cut algorithms.
// It allows adding which families to separate and their routines.
// Also for each family we can specify
//	- cut_limit: maximum number of cuts to add.
//	- iteration_limit: maximum number of iterations to run.
//	- cuts_per_iteration: maximum number of cuts to add per iteration.
//	- node_limit: maximum number of nodes where this family is separated.
//	- dependencies: families which must not find a cut in order to search for this family.
// Example:
// {
// 		"families": ["sec", "clique"],
// 		"dependencies": "clique<sec",
// 		"cut_limit": "sec:100,clique:1000",
// 		"iteration_limit": "*:50"
//	}
// This json serialization of a strategy specifies that there are two families "sec" and "clique".
// 	i) clique will only be separated if no secs are found on that round.
//	ii) only 100 secs will be added to the model and 1000 clique.
//	iii) up to 50 separation iterations will be executed for each family.
class SeparationStrategy : public Printable
{
public:
	// Strategy with no limits and no families.
	SeparationStrategy();
	
	// Add a cut family to the separation strategy.
	void AddFamily(const std::string& family);
	
	// Defines the separation routine that is executed for 'family'.
	void SetSeparationRoutine(const std::string& family, const SeparationRoutine* routine);
	
	// Returns: the separation routine defined for the family.
	const SeparationRoutine* SeparationRoutineFor(const std::string& family) const;
	
	// Adds a dependency of the type (f1, f2), this implies that f1 will only be separated if no f2 cuts were found.
	void AddDependency(const std::string& f1, const std::string& f2);
	
	// Returns: a sequence with all the cut families considered.
	const std::vector<std::string>& Families() const;
	
	// Returns: if f1 has a dependency on f2. Meaning that f1 will only be separated if no f2 cuts were found.
	bool HasDependency(const std::string& f1, const std::string& f2) const;
	
	// Returns: dependencies of f1.
	const std::vector<std::string>& Dependencies(const std::string& f1) const;
	
	// Prints the strategy.
	// Format: a json object with keys {families, cut_limit, iteration_limit, cuts_per_iteration, node_limit, dependencies}.
	// {cut_limit, iteration_limit, cuts_per_iteration, node_limit} are strings "f1:num1, f2:num2, ..." where fi is the
	// 	family and numi is the limit. If fi == "*" then all families are affected.
	// 	dependencies is a string "f1<f2,f3<f4,...,fk<fk+1" with all the dependencies. fi<fj is fi depends on fj.
	virtual void Print(std::ostream& os) const;
	
	std::unordered_map<std::string, int> cut_limit; // maximum number of cuts to add in total for a family.
	std::unordered_map<std::string, int> iteration_limit; // maximum number of iterations to run for a family.
	std::unordered_map<std::string, int> cuts_per_iteration; // maximum number of cuts to add per iteration for a family.
	std::unordered_map<std::string, int> node_limit; // maximum number of nodes where cuts will be searched for a family.
	std::unordered_map<std::string, double> improvement; // cuts will be searched if the objective value improved at least this value in the last iteration.
	
private:
	std::unordered_map<std::string, const SeparationRoutine*> routines_; // one separation routine associated to each cut family.
	std::vector<std::string> families_; // family of cuts that should be separated.
	std::unordered_map<std::string, std::vector<std::string>> dependencies_; // dependencies[f1] = set of all families
																			 // that should not find cuts in order to
																			 // separate f1.
};

// Sets the json object j with the serialized representation explained in the Print method of s.
void to_json(nlohmann::json& j, const SeparationStrategy& s);

// Parses back from the json representation to a strategy.
void from_json(const nlohmann::json& j, SeparationStrategy& s);
} // namespace goc

#endif //GOC_LINEAR_PROGRAMMING_CUTS_SEPARATION_STRATEGY_H
