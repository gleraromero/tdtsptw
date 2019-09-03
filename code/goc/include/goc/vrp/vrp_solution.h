//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_VRP_VRP_SOLUTION_H
#define GOC_VRP_VRP_SOLUTION_H

#include <vector>
#include <iostream>

#include "goc/lib/json.hpp"
#include "goc/print/printable.h"
#include "goc/vrp/route.h"

namespace goc
{
// Represents a solution to a Vehicle Routing Problem.
// This solution has a value, and a set of routes.
// - It knows how te serialize itself in JSON to be compatible with Kaleidoscope kd_type "vrp_solution".
class VRPSolution : public Printable
{
public:
	double value; // Value associated with the solution.
	std::vector<Route> routes; // Solution routes.
	
	VRPSolution() = default;
	
	// Creates the solution with the specified parameters.
	VRPSolution(double value, const std::vector<Route>& routes);
	
	// Prints the JSON representation of the solution.
	virtual void Print(std::ostream& os) const;
};

// Serializes the solution.
void to_json(nlohmann::json& j, const VRPSolution& solution);

// Parses an solution.
void from_json(const nlohmann::json& j, VRPSolution& solution);

} // namespace goc

#endif //GOC_VRP_VRP_SOLUTION_H
