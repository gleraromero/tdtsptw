//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef TDTSPTW_PRICING_PROBLEM_H
#define TDTSPTW_PRICING_PROBLEM_H

#include <vector>
#include <iostream>

#include "goc/goc.h"

#include "spf.h"

namespace tdtsptw
{
class PricingProblem : public goc::Printable
{
public:
	std::vector<double> penalties;
	
	// Returns: the reduced cost of graph path p with cost c_p.
	double PathReducedCost(const goc::GraphPath& p, double c_p) const;
	
	virtual void Print(std::ostream& os) const;
};

void to_json(nlohmann::json& j, const PricingProblem& p);

} // namespace tdtsptw

#endif //TDTSPTW_PRICING_PROBLEM_H
