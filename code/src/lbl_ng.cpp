//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "lbl_ng.h"

using namespace std;
using namespace goc;

namespace tdtsptw
{
// Generates all subsets of V and puts them into S.
void generate_subsets(const vector<Vertex>& V, vector<VertexSet>& S, VertexSet s=0, int i=0)
{
	if (i == V.size()) { S.push_back(s); return; }
	generate_subsets(V, S, unite(s, {V[i]}), i+1);
	generate_subsets(V, S, s, i+1);
}

NGStructure::NGStructure(const VRPInstance& vrp, int delta) : delta(delta)
{
	int n = vrp.D.VertexCount();
	
	// Compute Neighbourhoods.
	N = vector<VertexSet>(vrp.D.VertexCount());
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
			N[i].set(N_by_dist[k].second);
			N_list[i].push_back(N_by_dist[k].second);
		}
	}
	
	// Compute NGL path
	Digraph P(vrp.D.VertexCount());
	for (Vertex v: vrp.D.Vertices())
		for (Vertex w: vrp.D.Vertices())
			if (vrp.prec[v][w])
				P.AddArc({v, w});
	L = longest_path(P, vrp.o, vrp.d);
	
	VertexSet V_L;
	for (Vertex v: L) V_L.set(v);
	
	// Compute NGL sets V_i = { v \in V : !(v < L_i) and !(L_{i+1} < v) }.
	V = vector<VertexSet>(L.size());
	for (int i = 0; i < (int)L.size()-1; ++i)
	{
		V[i].set(L[i+1]);
		for (Vertex v: vrp.D.Vertices())
			if (!V_L.test(v) && !vrp.prec[v][L[i]] && !vrp.prec[L[i+1]][v])
				V[i].set(v);
	}
	
	// Compute NGSets.
	NGSet = vector<vector<VertexSet>>(n);
	for (Vertex v: vrp.D.Vertices()) generate_subsets(N_list[v], NGSet[v], create_bitset<MAX_N>({v}));
	
	// Compute NGSub.
	NGSub = vector<vector<vector<int>>>(n);
	for (Vertex v: vrp.D.Vertices())
	{
		for (int i = 0; i < NGSet[v].size(); ++i)
		{
			NGSub[v].push_back(vector<int>());
			for (int j = 0; j < NGSet[v].size(); ++j)
				if (i != j && is_subset(NGSet[v][j], NGSet[v][i]))
					NGSub[v][i].push_back(j);
		}
	}
	
	// Compute NGArc.
	NGArc = vector<vector<unordered_map<Vertex, int>>>(n);
	for (Vertex v: vrp.D.Vertices())
	{
		NGArc[v] = vector<unordered_map<Vertex, int>>(NGSet[v].size());
		for (int i = 0; i < NGSet[v].size(); ++i)
		{
			for (Vertex w: vrp.D.Successors(v))
			{
				if (NGSet[v][i].test(w)) continue;
				VertexSet NGw = unite(intersection(NGSet[v][i], N[w]), {w});
				for (int j = 0; j < NGSet[w].size(); ++j)
				{
					if (NGSet[w][j] == NGw)
					{
						NGArc[v][i].insert({w, j});
						break;
					}
				}
			}
		}
	}
}

NGLabel::NGLabel(NGLabel* prev, Vertex v, int S, double Tdur, double Thelp, double lambda) : prev(prev), v(v), S(S), Tdur(Tdur), Thelp(Thelp), lambda(lambda)
{}

GraphPath NGLabel::Path() const
{
	if (!prev) return {v};
	auto p = prev->Path();
	p.push_back(v);
	return p;
}

void NGLabel::Print(ostream& os) const
{
	os << "{P:" << Path() << ", S: " << S << ", Tdur: " << Tdur << ", Thelp: " << Thelp << ", lambda: " << lambda << "}";
}

BoundingStructure::BoundingStructure(VRPInstance* vrp, NGStructure* NG, const vector<double>& penalties)
	: vrp(vrp), NG(NG), penalties(penalties)
{
	// Init structure.
	int n = vrp->D.VertexCount();
	S = Matrix<vector<VectorMap<double, vector<NGLabel>>>>(n+1, n, vector<VectorMap<double, vector<NGLabel>>>(NG->L.size()));
	penalties_sum = sum(penalties);
}

void BoundingStructure::AddBound(int k, int r, const NGLabel& l)
{
	insert_sorted(S[k][l.v][r].Insert(floor(l.Thelp), {}), l, [] (const NGLabel& l1, const NGLabel& l2) { return l1.Tdur < l2.Tdur; });
}

double BoundingStructure::CompletionBound(double Tlambda, const goc::PWLFunction& Tdur, const VertexSet& V, int k, Vertex w)
{
	int n = vrp->D.VertexCount();
	double T = vrp->T;
	
	double LB = UB;
	double Thelp = -(max(dom(Tdur))-Tdur.Value(max(dom(Tdur))))-Tlambda;
	double min_dur = min(img(Tdur))-Tlambda;
	// Bounding labels should have visited all L[r] such that L[r] \not\in S, and not more. It also should
	// have visited L[r'] if there is a L[r'] = w.
	int r = 0;
	for (r = 0; r < (int)NG->L.size()-1 && (!V.test(NG->L[r+1]) || NG->L[r+1] == w); ++r) {}
	
	for (auto& Thelp_entry: S[n-k+1][w][r])
	{
		double Thelpm = Thelp_entry.first;
		if (epsilon_bigger(T + Thelpm + Thelp + penalties[w] + penalties_sum, LB)) break;
		for (auto& m: Thelp_entry.second)
		{
			if (epsilon_bigger(m.Tdur + min_dur + penalties[w] + penalties_sum, LB)) break;
			if (intersection(NG->NGSet[w][m.S], V) != create_bitset<MAX_N>({w})) continue;
			LB = min(LB, max(m.Tdur + min_dur, T + m.Thelp + Thelp) + penalties[w] + penalties_sum);
		}
	}
	return LB;
}

vector<Route> run_ng(const VRPInstance& vrp, const NGStructure& NG, const vector<double>& lambda, double UB,
					 Route* best_route, double* best_cost, MLBExecutionLog* log, BoundingStructure* B)
{
	// Return set of negative cost routes in NEG.
	vector<Route> NEG;
	
	int n = vrp.D.VertexCount();
	int R = NG.L.size();
	*best_route = Route({}, 0.0, INFTY);
	*best_cost = INFTY;
	
	stretch_to_size(*log->count_by_length, n+1, 0);
	Stopwatch rolex(true), rolex_domination(false), rolex_extension(false), rolex_queuing(false);
	
	// Initialize queue.
	Matrix<vector<vector<NGLabel>>> q(n, R, vector<vector<NGLabel>>(n));
	q[1][0][vrp.o].push_back(NGLabel(nullptr, vrp.o, 0, -lambda[vrp.o], -vrp.b[vrp.o]-lambda[vrp.o], lambda[vrp.o]));
	for (int k = 1; k < n; ++k)
	{
		for (int r = 0; r < R-1; ++r)
		{
			for (int v = 0; v < n; ++v)
			{
				rolex_queuing.Resume();
				sort(q[k][r][v].begin(), q[k][r][v].end(), [] (NGLabel& l1, NGLabel& l2) { return make_tuple(l1.Tdur, l1.Thelp, l1.S) < make_tuple(l2.Tdur, l2.Thelp, l1.S); });
				rolex_queuing.Pause();
				
				// D[i][S'] = minimum duration of a label M with: k(M)=k, r(M)=r, i(M)=i, S(M)=S'.
				vector<double> D(NG.NGSet[v].size(), INFTY);
				for (NGLabel& l: q[k][r][v])
				{
					bool is_dominated = false;
					rolex_domination.Resume();
					if (epsilon_smaller_equal(D[l.S], l.Thelp)) is_dominated = true;
					if (!is_dominated)
					{
						D[l.S] = l.Thelp;
						for (int sub: NG.NGSub[v][l.S])
						{
							if (epsilon_smaller_equal(D[sub], l.Thelp))
							{
								is_dominated = true;
								break;
							}
						}
					}
					rolex_domination.Pause();
					if (is_dominated)
					{
						log->dominated_count++;
						continue;
					}
					
					// Save it to B.
					if (B) B->AddBound(k, r, l);
					
					// Extension.
					log->processed_count++;
					rolex_extension.Resume();
					(*log->count_by_length)[k]++;
					
					for (pair<Vertex, int> w_Sw: NG.NGArc[v][l.S])
					{
						Vertex w = w_Sw.first;
						int Sw = w_Sw.second;
						
						// Feasibility check.
						if (k < vrp.prec_count[w]) continue;
						if (n-k-1 < vrp.suc_count[w]) continue;
						if (n-k-1 < vrp.suc_count[NG.L[r+1]]) continue;
						if (!NG.V[r].test(w)) continue;
						double d_vw = vrp.MinimumTravelTime({v,w});
						NGLabel lw(&l, w, Sw,
								   max(l.Tdur+d_vw, l.Thelp+vrp.a[w])-lambda[w], // Tdur(lw)
								   max(l.Tdur+d_vw-vrp.b[w], l.Thelp)-lambda[w], // Thelp(lw)
								   l.lambda + lambda[w]
						);
						
						// Process tour.
						if (w == vrp.d)
						{
							(*log->count_by_length)[k+1]++;
							if (lw.Tdur < *best_cost)
							{
								*best_cost = lw.Tdur;
								*best_route = Route(lw.Path(), -lw.Thelp + lw.lambda, lw.Tdur + lw.lambda);
							}
							if (epsilon_smaller(lw.Tdur, 0.0))
							{
								NEG.emplace_back(Route(lw.Path(), -lw.Thelp + lw.lambda, lw.Tdur + lw.lambda));
							}
							continue;
						}
						
						log->extended_count++;
						q[k+1][r + (w == NG.L[r+1])][w].push_back(lw);
					}
					rolex_extension.Pause();
				}
			}
		}
	}
	log->queuing_time = rolex_queuing.Peek();
	log->domination_time = rolex_domination.Peek();
	log->extension_time = rolex_extension.Peek();
	log->time = rolex.Peek();
	if (*best_cost == INFTY)
	{
		clog << "penalties: " << lambda << endl;
		fail("Should always be a best cost");
	}
	return NEG;
}
} // namespace tdtsptw