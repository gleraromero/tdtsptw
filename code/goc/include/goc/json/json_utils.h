//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_JSON_JSON_UTILS_H
#define GOC_JSON_JSON_UTILS_H

#include <string>

#include "goc/lib/json.hpp"

namespace goc
{
// Returns: if the object has the specified key defined.
bool has_key(const nlohmann::json& object, const std::string& key);

// Returns: The object value for the key if it is defined, otherwise returns def.
const nlohmann::json& value_or_default(const nlohmann::json& object, const std::string& key, const nlohmann::json& def);
} // namespace goc

// Add implementations for the to_json of common objects.
namespace std
{
// to_json implementation of a generic vector of elements.
template<typename T>
void to_json(nlohmann::json& j, const std::vector<T>& v)
{
	j = vector<nlohmann::json>();
	for (auto& e: v) j.push_back(e);
}
} // namespace std

#endif //GOC_JSON_JSON_UTILS_H
