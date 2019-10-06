//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "lbl_ng.h"

#include "pwl_domination_function.h"

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

NGStructure::NGStructure(const VRPInstance& vrp, const std::vector<VertexSet>& N, const goc::GraphPath& L, int delta)
	: N(N), L(L), delta(delta)
{
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
}

NGLabel::NGLabel(NGLabel* prev, Vertex v, int S, double Ttime, double Tdur, double Thelp, double lambda) : prev(prev), v(v), S(S), Ttime(Ttime), Tdur(Tdur), Thelp(Thelp), lambda(lambda)
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

vector<Route> run_ngl2res(const VRPInstance& vrp, const NGStructure& NG, const vector<double>& lambda, double UB,
					 Route* best_route, double* best_cost, MLBExecutionLog* log)
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
	q[1][0][vrp.o].push_back(NGLabel(nullptr, vrp.o, 0, 0.0, -lambda[vrp.o], -vrp.b[vrp.o]-lambda[vrp.o], lambda[vrp.o]));
	for (int k = 1; k < n; ++k)
	{
		for (int r = 0; r < R-1; ++r)
		{
			for (int v = 0; v < n; ++v)
			{
				rolex_queuing.Resume();
				sort(q[k][r][v].begin(), q[k][r][v].end(), [] (NGLabel& l1, NGLabel& l2) { return make_tuple(l1.Tdur, l1.Thelp, l1.S) < make_tuple(l2.Tdur, l2.Thelp, l2.S); });
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
								   0.0,
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

class TIState : public Printable
{
public:
	class Piece : public Printable
	{
	public:
		double l, r;
		double cost;

		Piece(double l, double r, double cost) : l(l), r(r), cost(cost)
		{ }

		virtual void Print(ostream& os) const
		{
			os << "(" << l << ", " << r << ", " << cost << ")";
		}

		bool Dominate(Piece& p)
		{
			if (epsilon_smaller(p.cost, cost)) return false;
			double ll = l;
			double rr = r + p.cost - cost;
			if (epsilon_bigger(ll, p.l) && epsilon_smaller(rr, p.r)) return false;
			if (epsilon_smaller_equal(ll, p.l) && epsilon_bigger_equal(rr, p.r)) return true;
			else if (epsilon_smaller_equal(ll, p.l)) p.l = max(p.l, rr);
			else if (epsilon_bigger_equal(rr, p.r)) p.r = min(p.r, ll);
			return false;
		}
	};

	void Merge(vector<Piece>& L2)
	{
		auto L1 = L;
		L.clear();

		// Sweep pieces.
		int i = 0, j = 0, k = 0;
		while (i < L1.size() && j < L2.size())
		{
			if (k++ > 10000)
			{
				clog.precision(17);
				clog << "Merge: " << k << endl;
				clog << L1[i] << " " << L2[j] << endl;
//				fail("Merge");
			}
			auto& p1 = L1[i];
			auto& p2 = L2[j];
			if (epsilon_smaller_equal(p1.l, p2.l) && p1.Dominate(p2)) { ++j; continue; }
			if (epsilon_smaller_equal(p2.l, p1.l) && p2.Dominate(p1)) { ++i; continue; }
			if (epsilon_smaller(p1.l, p2.l))
			{
				if (epsilon_smaller_equal(p1.r, p2.l))
				{
					L.push_back(p1);
					++i;
				}
				else
				{
					L.push_back(Piece(p1.l, p2.l, p1.cost));
					p1.l = p2.l+EPS;
				}
			}
			else
			{
				if (epsilon_smaller_equal(p2.r, p1.l))
				{
					L.push_back(p2);
					++j;
				}
				else
				{
					L.push_back(Piece(p2.l, p1.l, p2.cost));
					p2.l = p1.l+EPS;
				}
			}
		}
		// Add remaining pieces.
		L.insert(L.end(), L1.begin()+i, L1.end());
		L.insert(L.end(), L2.begin()+j, L2.end());
	}

	void DominateBy(const vector<Piece>& P)
	{
		auto L1 = P;
		auto L2 = L;
		L.clear();

		// Sweep pieces.
		int i = 0, j = 0, k = 0;
		while (i < L1.size() && j < L2.size())
		{
			if (k++ > 10000)
			{
				clog.precision(17);
				clog << "DominateBy: " << k << endl;
				clog << L1[i] << " " << L2[j] << endl;
				fail("DominateBy");
			}
			auto& p1 = L1[i];
			auto& p2 = L2[j];
			if (epsilon_smaller_equal(p1.l, p2.l) && p1.Dominate(p2)) { ++j; continue; }
			if (epsilon_smaller_equal(p2.l, p1.l) && p2.Dominate(p1)) { ++i; continue; }
			if (epsilon_smaller(p1.l, p2.l))
			{
				if (epsilon_smaller_equal(p1.r, p2.l))
				{
					++i;
				}
				else
				{
					p1.l = p2.l+EPS;
				}
			}
			else
			{
				if (epsilon_smaller_equal(p2.r, p1.l))
				{
					L.push_back(p2);
					++j;
				}
				else
				{
					L.push_back(Piece(p2.l, p1.l, p2.cost));
					p2.l = p1.l+EPS;
				}
			}
		}
		// Add remaining pieces.
		L.insert(L.end(), L2.begin()+j, L2.end());
	}

	virtual void Print(ostream& os) const
	{
		os << L;
	}

	vector<Piece> L;
};

double run_nglti(const VRPInstance& vrp, const NGStructure& NG, const vector<double>& lambda, double UB,
		Route& best_route, double& best_cost, MLBExecutionLog* log)
{
	int n = vrp.D.VertexCount();
	int R = NG.L.size();
	best_route = Route({}, 0.0, INFTY);
	best_cost = INFTY;
	double best_arrival = INFTY;

	stretch_to_size(*log->count_by_length, n+1, 0);
	Stopwatch rolex(true), rolex_domination(false), rolex_extension(false), rolex_queuing(false);

	Matrix<double> tt(n, n);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			tt[i][j] = vrp.MinimumTravelTime({i,j});

	// Initialize queue.
	Matrix<vector<vector<TIState>>> D(n, R, vector<vector<TIState>>(n, vector<TIState>(8))); // TODO: Change 8 by actual value.
	vector<TIState::Piece> initial_state = {TIState::Piece(vrp.a[vrp.o], vrp.b[vrp.o], -lambda[vrp.o])};
	D[1][0][vrp.o][0].Merge(initial_state);
	for (int k = 1; k < n; ++k)
	{
		for (int r = 0; r < R-1; ++r)
		{
			for (int v = 0; v < n; ++v)
			{
				for (int S = 0; S < NG.NGSet[v].size(); ++S)
				{
					auto& Delta = D[k][r][v][S];
					if (Delta.L.empty()) continue;

					// Dominate by subsets.
					for (int sub: NG.NGSub[v][S])
					{
						Delta.DominateBy(D[k][r][v][sub].L);
						if (Delta.L.empty()) continue;
					}

					double Ttimel = Delta.L.front().l;

					// Extend.
					for (pair<Vertex, int> w_Sw: NG.NGArc[v][S])
					{
						int w = w_Sw.first;
						int Sw = w_Sw.second;

						// Feasibility check.
						if (k < vrp.prec_count[w]) continue;
						if (n-k-1 < vrp.suc_count[w]) continue;
						if (w != NG.L[r + 1] && vrp.suc_count[NG.L[r + 1]] > n - k - 2) continue;
						if (!NG.V[r].test(w)) continue;
						double Ttimew = vrp.ArrivalTime({v,w}, Ttimel);
						if (Ttimew == INFTY) continue;
						vector<TIState::Piece> EXT;
						for (auto& p: Delta.L)
						{
							if (epsilon_bigger(p.l + tt[v][w], vrp.b[w])) break;
							double ll = max(p.l + tt[v][w], Ttimew);
							double rr = min(max(p.r + tt[v][w], Ttimew), vrp.b[w]);
							double cost = p.cost + (epsilon_equal(rr, Ttimew) ? Ttimew - p.r : tt[v][w]) - lambda[w];

							if (w == vrp.d)
							{
								(*log->count_by_length)[k+1]++;
								if (cost < best_cost)
								{
									best_cost = min(best_cost, cost);
									best_arrival = ll;
								}
								continue;
							}
							EXT.emplace_back(TIState::Piece(ll, rr, cost));
						}

						if (!EXT.empty())
						{
							D[k + 1][r + (NG.L[r + 1] == w)][w][Sw].Merge(EXT);
						}
					}
				}
			}
		}
	}

	GraphPath P = {vrp.d};
	int S = 0;
	int r = (int)NG.L.size()-2;
	Vertex w = vrp.d;
	double t = best_arrival;
	double c = best_cost;
	for (int k = n-1; k >= 0; --k)
	{
		bool found = false;
		for (Vertex v: vrp.D.Predecessors(w))
		{
			for (int Sv = 0; Sv < NG.NGSet[v].size(); Sv++)
			{
				if (NG.NGArc[v][Sv].at(w) != S) continue;
				for (auto& p: D[k][r][v][Sv].L)
				{
					found |= epsilon_smaller(p.r + tt[v][w], t) && epsilon_smaller_equal(p.cost+t-p.r-lambda[w], c);
					found |= epsilon_smaller_equal(p.l + tt[v][w], t) && epsilon_bigger_equal(p.r + tt[v][w], t) && epsilon_smaller_equal(p.cost+tt[v][w]-lambda[w], c);
					if (found)
					{
						t = min(p.r, t - tt[v][w]);
						c = p.cost;
						S = Sv;
						P.push_back(v);
						w = v;
						if (v == NG.L[r-1]) r--;
						break;
					}
				}
				if (found) break;
			}
			if (found) break;
		}
		if (!found) fail("Could not rebuild path in NGLTI.");
	}
	P = reverse(P);

	// Calculate path departure.
	double P_t0 = best_arrival;
	for (int k = n-1; k > 0; --k)
		P_t0 = min(vrp.b[P[k-1]], P_t0 - tt[P[k-1]][P[k]]);
	double P_duration = best_arrival - P_t0;
	double P_lambda = sum<Vertex>(P, [&](Vertex v) { return lambda[v]; });
	if (epsilon_bigger(P_duration - P_lambda, best_cost)) fail("Wrong cost");
	best_cost = P_duration - P_lambda;
	best_route = Route(P, P_t0, P_duration);

	log->queuing_time = rolex_queuing.Peek();
	log->domination_time = rolex_domination.Peek();
	log->extension_time = rolex_extension.Peek();
	log->time = rolex.Peek();
	if (best_cost == INFTY)
	{
		clog << "penalties: " << lambda << endl;
		fail("Should always be a best cost");
	}
	return best_cost + sum(lambda);
}
} // namespace tdtsptw