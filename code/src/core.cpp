//
// Created by Gonzalo Lera Romero on 28/04/2020.
//

#include "core.h"

using namespace std;
using namespace goc;

namespace tdtsptw
{
Core::Core(int k, int r, goc::Vertex v, const VertexSet& S)
	: k(k), r(r), v(v), S(S)
{}
} // namespace tdtsptw