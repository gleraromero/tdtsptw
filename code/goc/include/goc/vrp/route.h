//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_VRP_ROUTE_H
#define GOC_VRP_ROUTE_H

#include <iostream>

#include "goc/graph/graph_path.h"
#include "goc/lib/json.hpp"
#include "goc/print/printable.h"

namespace goc
{
// In a vehicle routing problem, a Route is a path traversed at a certain time.
// To define a route we need the path and the departing time t0. We keep also the duration for more information.
class Route : public Printable
{
public:
	GraphPath path; // Path traversed by the vehicle.
	double t0; // Departure time.
	double duration; // Route duration departing at t0.
	
	// Initializes the empty route with t0=0, duration=0.
	Route();
	
	// Initializes route with specified parameters.
	Route(const GraphPath& path, double t0, double duration);
	
	// Prints the JSON representation of the Route.
	virtual void Print(std::ostream& os) const;
};

// Format: {"path": ..., "t0": ..., "duration": ...}
void to_json(nlohmann::json& j, const Route& r);
void from_json(const nlohmann::json& j, Route& r);

} // namespace goc

#endif //GOC_VRP_ROUTE_H
