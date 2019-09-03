//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/log/log.h"

#include <unordered_map>
#include <string>

using namespace std;

namespace goc
{
void Log::Print(std::ostream& os) const
{
	os << ToJSON();
}

void to_json(nlohmann::json& j, const Log& l)
{
	j = l.ToJSON();
}
} // namespace goc