//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/linear_programming/cuts/separation_algorithm.h"

#include <set>
#include <climits>

#include "goc/time/stopwatch.h"
#include "goc/collection/collection_utils.h"
#include "goc/math/number_utils.h"

using namespace std;

namespace goc
{
SeparationAlgorithm::SeparationAlgorithm(const SeparationStrategy& separation_strategy)
	: strategy_(separation_strategy), is_disabled_(false)
{
	// Calculate topological order of families.
	families_ordered_by_dependencies_ = strategy_.Families();
	sort(families_ordered_by_dependencies_.begin(), families_ordered_by_dependencies_.end(),
		[&] (const string& f1, const string& f2) { return strategy_.HasDependency(f2, f1); });
	
	// Init tracking limit structure.
	for (auto& family: families_ordered_by_dependencies_)
	{
		cuts_added_[family] = iteration_count_[family] = 0;
		separation_time_[family] = 0.0_sec;
	}
	last_objective_ = INFTY;
}

vector<Constraint> SeparationAlgorithm::Separate(const Valuation& solution, int node_number, double node_bound) const
{
	vector<Constraint> cuts;
	if (is_disabled_) return cuts;
	if (disabled_families_.size() == strategy_.Families().size()) return cuts;
	
	// Adquire lock.
	lock_.lock();
	
	// Keep track of cut families that found violated cuts for dependencies purposes.
	unordered_set<string> families_with_cuts;
	
	// Try to find violated cuts for all families.
	for (auto& family: families_ordered_by_dependencies_)
	{
		if (includes(disabled_families_, family)) continue;
		
		// Check if family should be disabled.
		if (strategy_.node_limit.at(family) <= node_number) disabled_families_.insert(family);
		else if (strategy_.cut_limit.at(family) <= cuts_added_[family]) disabled_families_.insert(family);
		else if (strategy_.iteration_limit.at(family) <= iteration_count_.at(family)) disabled_families_.insert(family);
		
		if (includes(disabled_families_, family)) continue;
		
		// Check if all the dependecies of the family of inequalities have failed to find cuts. If so, then we proceed
		// to find cuts.
		bool all_dependencies_failed = true;
		for (auto& dependency: strategy_.Dependencies(family))
		{
			if (includes(families_with_cuts, dependency))
			{
				all_dependencies_failed = false;
				break;
			}
		}
		if (!all_dependencies_failed) continue;
		
		// Check what is the max amount of cuts that can be added in this iteration.
		int cut_limit = CutLimitForThisIteration(family, node_bound);
		if (cut_limit == 0) continue;
		
		// Execute the separation routine.
		Stopwatch rolex(true);
		vector<Constraint> family_cuts = strategy_.SeparationRoutineFor(family)->Separate(solution, node_number, cut_limit, node_bound);
		rolex.Pause();
		for (int i = 0; i < min((int)family_cuts.size(), cut_limit); ++i) cuts.push_back(family_cuts[i]);
		if (!family_cuts.empty()) families_with_cuts.insert(family);
		
		// Keep track of statistics.
		cuts_added_[family] += min((int)family_cuts.size(), cut_limit);
		iteration_count_[family]++;
		separation_time_[family] += rolex.Peek();
	}
	last_objective_ = node_bound;
	lock_.unlock();
	return cuts;
}

bool SeparationAlgorithm::IsEnabled() const
{
	return !is_disabled_ && disabled_families_.size() < strategy_.Families().size();
}

int SeparationAlgorithm::CutsAdded() const
{
	int count = 0;
	for (auto& it: cuts_added_) count += it.second;
	return count;
}

int SeparationAlgorithm::CutsAdded(const string& family) const
{
	return cuts_added_.at(family);
}

int SeparationAlgorithm::IterationCount() const
{
	int count = 0;
	for (auto& it: iteration_count_) count += it.second;
	return count;
}

int SeparationAlgorithm::IterationCount(const string& family) const
{
	return iteration_count_.at(family);
}

Duration SeparationAlgorithm::SeparationTime() const
{
	Duration count;
	for (auto& it: separation_time_) count += it.second;
	return count;
}

Duration SeparationAlgorithm::SeparationTime(const string& family) const
{
	return separation_time_.at(family);
}

const SeparationStrategy& SeparationAlgorithm::Strategy() const
{
	return strategy_;
}

void SeparationAlgorithm::Disable() const
{
	is_disabled_ = true;
}

int SeparationAlgorithm::CutLimitForThisIteration(const string& family, double node_bound) const
{
	int limit = INT_MAX;
	limit = min(limit, strategy_.cut_limit.at(family) - cuts_added_.at(family));
	limit = min(limit, strategy_.iteration_limit.at(family));
	limit = min(limit, strategy_.cuts_per_iteration.at(family));
	if (fabs(node_bound-last_objective_) < strategy_.improvement.at(family)) limit = 0;
	return max(limit, 0);
}
} // namespace goc