//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/collection/collection_utils.h"

#include <numeric>

using namespace std;

namespace goc
{
vector<int> range(int left, int right)
{
	vector<int> v(right-left);
	iota(v.begin(), v.end(), left);
	return v;
}
} // namespace goc