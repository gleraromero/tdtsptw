//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_PRINT_PRINTABLE_H
#define GOC_PRINT_PRINTABLE_H

#include <iostream>

namespace goc
{
// Represents the interface for all objects that can be printed. It automatically implements the operator<< for objects
// that implement this interface.
class Printable
{
public:
	virtual void Print(std::ostream& os) const = 0;
};

std::ostream& operator<<(std::ostream& os, const Printable& object);
} // namespace goc

#endif //GOC_PRINT_PRINTABLE_H
