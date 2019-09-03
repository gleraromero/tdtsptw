//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_STRING_STRING_UTILS_H
#define GOC_STRING_STRING_UTILS_H

#include <sstream>
#include <string>
#include <vector>

// This is a macro that converts any object into a string by using its operator<<(ostream&).
// Example: STR(15) == "15".
// Example: STR(vector<int>({3, 5, 7})) == "[3, 5, 7]".
#define STR(x) ([&]() { std::ostringstream qwert; qwert << (x); return qwert.str();})()

namespace goc
{
// Returns: If the string 'container' includes the string 'contained'.
bool string_contains(const std::string& container, const std::string& contained);

// Returns: string s without the appearances of c.
std::string remove(const std::string& s, char c);

// Returns: A list of strings which result from splitting the string s by the delimiter.
// Example: split("Hello how are you", ' ') == ["Hello", "how", "are", "you"].
std::vector<std::string> split(const std::string& s, char delimiter=',');

// Removes the whitespaces at the beginning and end of s.
// Example: "    hello  " == "hello".
std::string trim(const std::string& s);
} // namespace goc

#endif // GOC_STRING_STRING_UTILS_H