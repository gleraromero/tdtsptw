//
// Created by Gonzalo Lera Romero on 24/04/2020.
//

#include "ngl_info.h"

using namespace std;
using namespace goc;

namespace tdtsptw
{
void create_default_nginfo(const VRPInstance& vrp, int delta, NGLInfo* forward, NGLInfo* backward)
{
	// Compute Neighbourhoods.
	forward->N = vector<VertexSet>(vrp.D.VertexCount());
	vector<vector<Vertex>> N_list(vrp.D.VertexCount());
	for (Vertex i: vrp.D.Vertices())
	{
		vector<pair<double, Vertex>> N_by_dist;
		for (Vertex j: vrp.D.Vertices())
		{
			if (i == j || j == vrp.o || j == vrp.d) continue;
			if (epsilon_bigger(vrp.EAT[j][i], vrp.tw[i].right)) continue;
			if (epsilon_bigger(vrp.EAT[i][j], vrp.tw[j].right)) continue;
			N_by_dist.push_back({vrp.MinimumTravelTime({i,j}), j});
		}
		sort(N_by_dist.begin(), N_by_dist.end());
		for (int k = 0; k < min(delta, (int)N_by_dist.size()); ++k)
		{
			forward->N[i].set(N_by_dist[k].second);
			N_list[i].push_back(N_by_dist[k].second);
		}
	}

	// Compute NGL path
	Digraph P(vrp.D.VertexCount());
	for (Vertex v: vrp.D.Vertices())
		for (Vertex w: vrp.D.Vertices())
			if (vrp.prec[v][w])
				P.AddArc({v, w});
	forward->L = longest_path(P, vrp.o, vrp.d);
	int L_size = forward->L.size();

	VertexSet V_L;
	for (Vertex v: forward->L) V_L.set(v);

	// Compute NGL sets V_i = { v \in V : !(v < L_i) and !(L_{i+1} < v) }.
	forward->V = vector<VertexSet>(L_size);
	for (int i = 0; i < L_size-1; ++i)
	{
		forward->V[i].set(forward->L[i+1]);
		for (Vertex v: vrp.D.Vertices())
			if (!V_L.test(v) && !vrp.prec[v][forward->L[i]] && !vrp.prec[forward->L[i+1]][v])
				forward->V[i].set(v);
	}

	// Assign backward fields.
	forward->delta = backward->delta = delta;
	backward->N = forward->N;
	backward->L = reverse(forward->L);
	backward->V = vector<VertexSet>(L_size);
	for (int i = 0; i < L_size-1; ++i)
	{
		backward->V[i] = forward->V[L_size-2-i];
		backward->V[i].reset(backward->L[i]);
		backward->V[i].set(backward->L[i+1]);
	}
}
} // namespace tdtsptw