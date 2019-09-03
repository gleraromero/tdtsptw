//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/json/json_utils.h"

using namespace std;
using namespace nlohmann;

namespace goc
{
bool has_key(const json& object, const string& key)
{
	return object.find(key) != object.end();
}

const json& value_or_default(const json& object, const string& key, const json& def)
{
	return has_key(object, key) ? object[key] : def;
}
} // namespace goc