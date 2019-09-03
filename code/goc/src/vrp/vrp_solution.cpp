//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/vrp/vrp_solution.h"

using namespace std;
using namespace nlohmann;

namespace goc
{
VRPSolution::VRPSolution(double value, const vector<Route>& routes)
	: value(value), routes(routes)
{

}

void VRPSolution::Print(ostream& os) const
{
	os << json(*this);
}

void to_json(json& j, const VRPSolution& solution)
{
	j["kd_type"] = "vrp_solution";
	j["value"] = solution.value;
	j["routes"] = solution.routes;
}

void from_json(const json& j, VRPSolution& solution)
{
	solution.value = j["value"];
	for (Route r: j["routes"]) solution.routes.push_back(r);
}

} // namespace goc