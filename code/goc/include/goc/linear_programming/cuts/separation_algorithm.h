//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_LINEAR_PROGRAMMING_CUTS_SEPARATION_ALGORITHM_H
#define GOC_LINEAR_PROGRAMMING_CUTS_SEPARATION_ALGORITHM_H

#include <mutex>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "goc/linear_programming/cuts/separation_strategy.h"
#include "goc/linear_programming/model/constraint.h"
#include "goc/time/duration.h"

namespace goc
{
// This class is responsible for executing a separation strategy on a branch and cut algorithm.
// It also keeps track of the cuts added and time spent for statistics.
class SeparationAlgorithm
{
public:
	// separation_strategy: Strategy of the cutting planes algorithm (node limit, iterations, etc.).
	SeparationAlgorithm(const SeparationStrategy& separation_strategy);
	
	// Separates the solution using the separation strategy.
	// solution: current fractional solution that has to be cut.
	// node_number: node in the BB tree enumeration (0 is root).
	// node_bound: bound of the current node in the BB tree.
	// Returns: the violated constraints found.
	std::vector<Constraint> Separate(const Valuation& solution, int node_number, double node_bound) const;
	
	// Returns: true if there is at least one cut family which can be executed for separation according to the
	// separation strategy. If all families stopped separating because of their limits, then it returns false.
	bool IsEnabled() const;
	
	// Returns: the total number of cuts added.
	int CutsAdded() const;
	
	// Returns: the total number of cuts from the family added.
	int CutsAdded(const std::string& family) const;
	
	// Returns: the total number of cut iterations performed.
	int IterationCount() const;
	
	// Returns: the total number of cut iterations performed from the family.
	int IterationCount(const std::string& family) const;
	
	// Returns: the total time spent on separation routines.
	Duration SeparationTime() const;
	
	// Returns: the total time spent on separation routines from the family.
	Duration SeparationTime(const std::string& family) const;
	
	// Returns: the separation strategy being used.
	const SeparationStrategy& Strategy() const;
	
	// Disables the cutting plane algorithm. This will prevent future cuts to be added.
	void Disable() const;
	
private:
	// Limit of cuts for a given family in the current iteration.
	int CutLimitForThisIteration(const std::string& family, double node_bound) const;
	
	mutable std::mutex lock_; 	// This lock is important for parallel branch and bound execution.
							 	// It makes sure only one separation routine is called at each time.
							 	
	SeparationStrategy strategy_; // Strategy to be used.
	std::vector<std::string> families_ordered_by_dependencies_; // Families in topological order by '<': "depends on".
	
	// Keep track of what happened so far.
	mutable double last_objective_; // last objective value recorded.
	mutable std::unordered_map<std::string, int> cuts_added_, iteration_count_; // cuts added and number of iterations run.
	mutable std::unordered_map<std::string, Duration> separation_time_; // Time spent on separation of each family.
	mutable std::unordered_set<std::string> disabled_families_; // families that reached some of the limits and will never separate a cut in the future.
	mutable bool is_disabled_; // True if the algorithm is disabled, false otherwise.
};
} // namespace goc

#endif //GOC_LINEAR_PROGRAMMING_CUTS_SEPARATION_ALGORITHM_H
