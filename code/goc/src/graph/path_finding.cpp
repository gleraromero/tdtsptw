//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/graph/path_finding.h"

#include <limits.h>

#include "goc/collection/collection_utils.h"
#include "goc/math/number_utils.h"

using namespace std;

namespace goc
{
GraphPath longest_path(const Digraph& D, Vertex s, Vertex t)
{
	// Calculate topological order.
	vector<Vertex> topo = range(0, D.VertexCount());
	sort(topo.begin(), topo.end(), [&] (Vertex v1, Vertex v2) { return D.IncludesArc({v2, v1}); });
	
	// Calculate max distance and parents.
	vector<int> max_dist(D.VertexCount(), -INT_MAX);
	vector<Vertex> parent(D.VertexCount(), -1);
	parent[t] = -1;
	max_dist[t] = 0;
	for (auto& v: topo)
	{
		for (Vertex w: D.Successors(v))
		{
			if (max_dist[w]+1 > max_dist[v])
			{
				max_dist[v] = max_dist[w]+1;
				parent[v] = w;
			}
		}
	}
	
	// Reconstruct the longest path.
	GraphPath L;
	for (Vertex v = s; v != -1; v = parent[v])
		L.push_back(v);
	
	return L;
}

vector<double> compute_earliest_arrival_time(const Digraph& D, Vertex s, double t0, const function<double(Vertex, Vertex, double)>& tt)
{
	priority_queue<pair<double, Vertex>, vector<pair<double, Vertex>>, greater<>> q;
	vector<bool> visited(D.VertexCount(), false);
	vector<double> EAT(D.VertexCount(), INFTY); // EAT[j] = Earliest arrival time to vertex j
	q.push({t0, s});
	while (!q.empty())
	{
		double t; Vertex v;
		tie(t, v) = q.top();
		q.pop();
		if (visited[v]) continue;
		visited[v] = true;
		EAT[v] = t;
		for (auto& w: D.Successors(v))
		{
			if (!visited[w])
			{
				double travel_time = tt(v, w, t);
				if (travel_time == INFTY) continue;
				q.push({t + travel_time, w});
			}
		}
	}
	return EAT;
}

vector<double> compute_latest_departure_time(const Digraph& D, Vertex s, double t0, const function<double(Vertex, Vertex, double)>& dep)
{
	priority_queue<pair<double, Vertex>> q;
	vector<bool> visited(D.VertexCount(), false);
	vector<double> LDT(D.VertexCount(), -INFTY);
	q.push({t0, s});
	while (!q.empty())
	{
		double t; Vertex v;
		tie(t, v) = q.top();
		q.pop();
		if (visited[v]) continue;
		visited[v] = true;
		if (t == INFTY) continue;
		LDT[v] = t;
		for (auto& w: D.Predecessors(v))
		{
			if (!visited[w])
			{
				double d = dep(w, v, t);
				if (d != INFTY) q.push({d, w});
			}
		}
	}
	return LDT;
}
} // namespace goc