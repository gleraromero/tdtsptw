//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "pricing_problem.h"

using namespace std;
using namespace goc;
using namespace nlohmann;

namespace tdtsptw
{
double PricingProblem::PathReducedCost(const goc::GraphPath& p, double c_p) const
{
	double rc = c_p;
	for (Vertex i: p) rc -= penalties[i];
	return rc;
}

void PricingProblem::Print(std::ostream& os) const
{
	os << json(*this);
}

void to_json(json& j, const PricingProblem& p)
{
	j["penalties"] = p.penalties;
}
} // namespace tdtsptw