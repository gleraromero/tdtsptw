//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/print/printable.h"

using namespace std;

namespace goc
{
std::ostream& operator<<(std::ostream& os, const Printable& object)
{
	object.Print(os);
	return os;
}
} // namespace goc