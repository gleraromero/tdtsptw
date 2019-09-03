//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_EXCEPTION_EXCEPTION_UTILS_H
#define GOC_EXCEPTION_EXCEPTION_UTILS_H

#include <exception>
#include <iostream>
#include <string>

namespace goc
{
// A call to this function will throw an exception with the given 'message'.
// The purpose of this function is to standardize the way we fail in the program.
inline void fail(const std::string& message, int exit_code = 1)
{
	std::cerr << "ERROR: " << message << std::endl;
	throw std::runtime_error(message);
	exit(exit_code);
}
}

#endif // GOC_EXCEPTION_EXCEPTION_UTILS_H