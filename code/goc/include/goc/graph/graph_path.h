//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_GRAPH_GRAPH_PATH_H
#define GOC_GRAPH_GRAPH_PATH_H

#include <limits.h>
#include <vector>

#include "goc/graph/vertex.h"

namespace goc
{
// Represents a path in a graph or digraph.
typedef std::vector<Vertex> GraphPath;

// Returns: if the path 'p' contains a cycle of size 'max_size' vertices or less.
bool has_cycle(GraphPath p, int max_size=INT_MAX);
} // namespace goc

#endif //GOC_GRAPH_GRAPH_PATH_H
