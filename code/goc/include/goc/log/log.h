//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_LOG_LOG_H
#define GOC_LOG_LOG_H

#include <iostream>

#include "goc/lib/json.hpp"
#include "goc/print/printable.h"

namespace goc
{
// Represent the execution log of an algorithm. It must know how to serialize itself as a JSON object so to add it to
// the program output.
class Log : public Printable
{
public:
	// Maps the execution log to a json object.
	virtual nlohmann::json ToJSON() const = 0;
	
	virtual void Print(std::ostream& os) const;
};

void to_json(nlohmann::json& j, const Log& l);
} // namespace goc

#endif //GOC_LOG_LOG_H
