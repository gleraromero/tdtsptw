//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/linear_programming/cuts/separation_strategy.h"

#include <climits>
#include <stdlib.h>

#include "goc/collection/collection_utils.h"
#include "goc/json/json_utils.h"
#include "goc/math/number_utils.h"
#include "goc/string/string_utils.h"

using namespace std;
using namespace nlohmann;

namespace goc
{

SeparationStrategy::SeparationStrategy()
{
}

void SeparationStrategy::AddFamily(const string& family)
{
	families_.push_back(family);
	dependencies_[family] = {};
	cut_limit[family] = INT_MAX;
	iteration_limit[family] = INT_MAX;
	cuts_per_iteration[family] = INT_MAX;
	node_limit[family] = INT_MAX;
	improvement[family] = 0.0;
}

void SeparationStrategy::SetSeparationRoutine(const string& family, const SeparationRoutine* routine)
{
	routines_[family] = routine;
}

const SeparationRoutine* SeparationStrategy::SeparationRoutineFor(const string& family) const
{
	return routines_.at(family);
}

void SeparationStrategy::AddDependency(const string& f1, const string& f2)
{
	dependencies_[f1].push_back(f2);
}

const std::vector<string>& SeparationStrategy::Families() const
{
	return families_;
}

bool SeparationStrategy::HasDependency(const string& f1, const string& f2) const
{
	return includes(dependencies_.at(f1), f2);
}

const vector<string>& SeparationStrategy::Dependencies(const string& f1) const
{
	return dependencies_.at(f1);
}

void SeparationStrategy::Print(ostream& os) const
{
	os << json(*this);
}

void to_json(json& j, const SeparationStrategy& s)
{
	if (!s.Families().empty()) j["families"] = s.Families();

	string dependencies_str;
	for (auto& family: s.Families())
	{
		for (auto& dependency: s.Dependencies(family))
		{
			if (dependencies_str != "") dependencies_str += ",";
			dependencies_str += dependency + "<" + family;
		}
	}
	if (!dependencies_str.empty()) j["dependencies"] = dependencies_str;
	
	string cut_limit_str;
	for (auto& family: s.Families())
	{
		if (s.cut_limit.at(family) == INT_MAX) continue;
		if (cut_limit_str != "") cut_limit_str += ",";
		cut_limit_str += family + ":" + STR(s.cut_limit.at(family));
	}
	if (!cut_limit_str.empty()) j["cut_limit"] = cut_limit_str;
	
	string iteration_limit_str;
	for (auto& family: s.Families())
	{
		if (s.iteration_limit.at(family) == INT_MAX) continue;
		if (iteration_limit_str != "") iteration_limit_str += ",";
		iteration_limit_str += family + ":" + STR(s.iteration_limit.at(family));
	}
	if (!iteration_limit_str.empty()) j["iteration_limit"] = iteration_limit_str;
	
	string cuts_per_iteration_str;
	for (auto& family: s.Families())
	{
		if (s.cuts_per_iteration.at(family) == INT_MAX) continue;
		if (cuts_per_iteration_str != "") cuts_per_iteration_str += ",";
		cuts_per_iteration_str += family + ":" + STR(s.cuts_per_iteration.at(family));
	}
	if (!cuts_per_iteration_str.empty()) j["cuts_per_iteration"] = cuts_per_iteration_str;
	
	string node_limit_str = "";
	for (auto& family: s.Families())
	{
		if (s.node_limit.at(family) == INT_MAX) continue;
		if (node_limit_str != "") node_limit_str += ",";
		node_limit_str += family + ":" + STR(s.node_limit.at(family));
	}
	if (!node_limit_str.empty())j["node_limit"] = node_limit_str;
	
	string improvement_str = "";
	for (auto& family: s.Families())
	{
		if (s.improvement.at(family) == 0.0) continue;
		if (improvement_str != "") improvement_str += ",";
		improvement_str += family + ":" + STR(s.improvement.at(family));
	}
	if (!improvement_str.empty()) j["improvement"] = improvement_str;
}

void from_json(const json& j, SeparationStrategy& s)
{
	if (has_key(j, "families")) for (auto& family: j["families"]) s.AddFamily(family);
	// Parse cut_limit.
	// Example: {"cut_limit": "f1:50, f2: 60"}.
	// Example: {"cut_limit": "*:100"} <- * means all cuts.
	if (has_key(j, "cut_limit") && j["cut_limit"] != "")
	{
		auto limits = split(j["cut_limit"], ',');
		for (auto& limit: limits)
		{
			auto split_limit = split(limit, ':');
			string f = trim(split_limit[0]);
			int l = atoi(trim(split_limit[1]).c_str());
			if (f == "*")
			{
				for (auto& family: s.Families())
					s.cut_limit[family] = l;
			}
			else
			{
				s.cut_limit[f] = l;
			}
		}
	}
	// Parse iteration_limit.
	// Example: {"iteration_limit": "f1:50, f2: 60"}.
	// Example: {"iteration_limit": "*:100"} <- * means all cuts.
	if (has_key(j, "iteration_limit") && j["iteration_limit"] != "")
	{
		auto limits = split(j["iteration_limit"], ',');
		for (auto& limit: limits)
		{
			auto split_limit = split(limit, ':');
			string f = trim(split_limit[0]);
			int l = atoi(trim(split_limit[1]).c_str());
			if (f == "*")
			{
				for (auto& family: s.Families())
					s.iteration_limit[family] = l;
			}
			else
			{
				s.iteration_limit[f] = l;
			}
		}
	}
	// Parse cuts_per_iteration.
	// Example: {"cuts_per_iteration": "f1:50, f2: 60"}.
	// Example: {"cuts_per_iteration": "*:100"} <- * means all cuts.
	if (has_key(j, "cuts_per_iteration") && j["cuts_per_iteration"] != "")
	{
		auto limits = split(j["cuts_per_iteration"], ',');
		for (auto& limit: limits)
		{
			auto split_limit = split(limit, ':');
			string f = trim(split_limit[0]);
			int l = atoi(trim(split_limit[1]).c_str());
			if (f == "*")
			{
				for (auto& family: s.Families())
					s.cuts_per_iteration[family] = l;
			}
			else
			{
				s.cuts_per_iteration[f] = l;
			}
		}
	}
	// Parse node_limit.
	// Example: {"node_limit": "f1:50, f2: 60"}.
	// Example: {"node_limit": "*:100"} <- * means all cuts.
	if (has_key(j, "node_limit") && j["node_limit"] != "")
	{
		auto limits = split(j["node_limit"], ',');
		for (auto& limit: limits)
		{
			auto split_limit = split(limit, ':');
			string f = trim(split_limit[0]);
			int l = atoi(trim(split_limit[1]).c_str());
			if (f == "*")
			{
				for (auto& family: s.Families())
					s.node_limit[family] = l;
			}
			else
			{
				s.node_limit[f] = l;
			}
		}
	}
	// Parse improvement.
	// Example: {"improvement": "f1:0.1, f2: 0.0"}.
	// Example: {"improvement": "*:0.1"} <- * means all cuts.
	if (has_key(j, "improvement") && j["improvement"] != "")
	{
		auto limits = split(j["improvement"], ',');
		for (auto& limit: limits)
		{
			auto split_limit = split(limit, ':');
			string f = trim(split_limit[0]);
			double l = atof(trim(split_limit[1]).c_str());
			if (f == "*")
			{
				for (auto& family: s.Families())
					s.improvement[family] = l;
			}
			else
			{
				s.improvement[f] = l;
			}
		}
	}
	
	// Parse dependencies.
	// Example: {"dependencies": "f1<f2, f3 < f4"}.
	if (has_key(j, "dependencies") && j["dependencies"] != "")
	{
		// Dependencies is a string with clauses f1<f2 (meaning family f2 depends on family f1), separated by commas.
		auto dependencies = split(j["dependencies"], ',');
		for (auto& dependency: dependencies)
		{
			auto splitted_dependency = split(dependency, '<');
			string f1 = trim(splitted_dependency[0]);
			string f2 = trim(splitted_dependency[1]);
			
			s.AddDependency(f2, f1);
		}
	}
}
} // namespace goc