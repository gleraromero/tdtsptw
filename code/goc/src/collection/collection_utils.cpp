//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/collection/collection_utils.h"

using namespace std;

namespace goc
{
vector<int> range(int left, int right)
{
	vector<int> v;
	for (int i = left; i < right; ++i) v.push_back(i);
	return v;
}
} // namespace goc