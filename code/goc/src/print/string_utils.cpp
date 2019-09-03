//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/string/string_utils.h"

#include <sstream>

using namespace std;

namespace goc
{
bool string_contains(const string& container, const string& contained)
{
	return container.find(contained) != string::npos;
}

string remove(const string& s, char c)
{
	string r = "";
	for (char d: s) if (d != c) r += d;
	return r;
}

vector<string> split(const string& s, char delimiter)
{
	stringstream ss(s);
	vector<string> result;
	while( ss.good() )
	{
		string substr;
		getline(ss, substr, delimiter);
		result.push_back(substr);
	}
	return result;
}

string trim(const string& s)
{
	string t;
	for (char c: s) if (c != ' ') t.append(1u, c);
	return t;
}
} // namespace goc