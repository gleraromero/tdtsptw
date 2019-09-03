//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/linear_programming/model/valuation.h"

#include "goc/lib/json.hpp"
#include "goc/math/number_utils.h"
#include "goc/collection/collection_utils.h"

using namespace std;
using namespace nlohmann;

namespace goc
{
void Valuation::SetValue(const Variable& v, double value)
{
	if (epsilon_equal(value, 0.0)) { values_.erase(v); return; }
	values_[v] = value;
}

double Valuation::operator[](const Variable& v) const
{
	return includes_key(values_, v) ? values_.at(v) : 0.0;
}

double Valuation::at(const Variable& v) const
{
	return (*this)[v];
}

bool Valuation::IsInteger() const
{
	for (auto& variable_value: values_)
		if (epsilon_different(variable_value.second, round(variable_value.second)))
			return false;
		
	return true;
}

void Valuation::Print(ostream& os) const
{
	os << json(*this);
}

unordered_map<Variable, double>::const_iterator Valuation::begin() const
{
	return values_.begin();
}

unordered_map<Variable, double>::const_iterator Valuation::end() const
{
	return values_.end();
}

void to_json(json& j, const Valuation& v)
{
	j = map<string, double>();
	for (auto& variable_value: v) j[variable_value.first.name] = variable_value.second;
}
} // namespace goc