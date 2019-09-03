//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/graph/graph_path.h"

#include <unordered_set>

using namespace std;

namespace goc
{
bool has_cycle(GraphPath p, int max_size)
{
	unordered_set<int> V;
	max_size = min(max_size, (int)p.size());
	
	// Add to V all vertices in (p[0], ..., p[max_size-1]).
	// If any vertex is repeated, then it has a cycle.
	for (int i = 0; i < max_size; ++i)
	{
		if (V.find(p[i]) != V.end()) return true;
		V.insert(p[i]);
	}
	
	// Now move the window (p[i-max_size+1], ..., p[i]) until the end i==|p|-1.
	for (int i = max_size; i < p.size(); ++i)
	{
		V.erase(p[i-max_size]);
		if (V.find(p[i]) != V.end()) return true;
		V.insert(p[i]);
	}
	
	return false;
}
} // namespace goc
