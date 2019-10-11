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
	
	AdjustTimeWindows(vrp);
}

NGStructure::NGStructure(const VRPInstance& vrp, const std::vector<VertexSet>& N, const goc::GraphPath& L, int delta)
	: L(L), delta(delta), N(N)
{
	int n = vrp.D.VertexCount();
	
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
	
	vector<vector<Vertex>> N_list(vrp.D.VertexCount());
	for (Vertex i: vrp.D.Vertices())
		for (Vertex j: vrp.D.Vertices())
			if (N[i].test(j))
				N_list[i].push_back(j);
		
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
	
	AdjustTimeWindows(vrp);
}

void NGStructure::AdjustTimeWindows(const VRPInstance& vrp)
{
	int n = vrp.D.VertexCount();
	tw = Matrix<Interval>(n+1, n);
	for (int k = 1; k <= n; ++k)
		for (Vertex v: vrp.D.Vertices())
			tw[k][v] = vrp.tw[v];
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
						if (w != NG.L[r + 1] && vrp.suc_count[NG.L[r + 1]] > n - k - 2) continue;
						if (!NG.V[r].test(w)) continue;
						double d_vw = vrp.MinimumTravelTime({v,w});
						double b_w = NG.tw[k+1][w].right;
						NGLabel lw(&l, w, Sw,
								   0.0,
								   max(l.Tdur+d_vw, l.Thelp+vrp.a[w])-lambda[w], // Tdur(lw)
								   max(l.Tdur+d_vw-b_w, l.Thelp)-lambda[w], // Thelp(lw)
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
		Piece* prev;
		Vertex v;

		Piece(double l, double r, double cost, Piece* prev, Vertex v) : l(l), r(r), cost(cost), prev(prev), v(v)
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
		
		GraphPath Path() const
		{
			if (!prev) return {v};
			auto P = prev->Path();
			P.push_back(v);
			return P;
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
				clog << "Merge: " << k << " - " << i << " " << j << endl;
				clog << L1[i] << " " << L2[j] << endl;
				fail("Merge");
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
					L.push_back(Piece(p1.l, p2.l, p1.cost, p1.prev, p1.v));
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
					L.push_back(Piece(p2.l, p1.l, p2.cost, p2.prev, p2.v));
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
					L.push_back(Piece(p2.l, p1.l, p2.cost, p2.prev, p2.v));
					p2.l = p1.l+EPS;
				}
			}
		}
		// Add remaining pieces.
		L.insert(L.end(), L2.begin()+j, L2.end());
	}
	
	void Normalize()
	{
		auto L2 = L;
		L.clear();
		for (auto& p: L2)
		{
			if (L.empty() || p.cost != L.back().cost || epsilon_smaller(L.back().r, p.l))
			{
				L.push_back(p);
			}
			else
			{
				L.back().r = p.r;
			}
		}
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
	TIState::Piece best_piece(0,0,0,nullptr, vrp.o);

	stretch_to_size(*log->count_by_length, n+1, 0);
	Stopwatch rolex(true), rolex_domination(false), rolex_extension(false), rolex_queuing(false);

	// Initialize queue.
	Matrix<vector<vector<TIState>>> D(n, R, vector<vector<TIState>>(n, vector<TIState>(8))); // TODO: Change 8 by actual value.
	vector<TIState::Piece> initial_state = {TIState::Piece(vrp.a[vrp.o], vrp.b[vrp.o], -lambda[vrp.o], nullptr, vrp.o)};
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
					log->processed_count++;
					
					// Dominate by subsets.
					rolex_domination.Resume();
					for (int sub: NG.NGSub[v][S])
					{
						if (Delta.L.empty()) break;
						Delta.DominateBy(D[k][r][v][sub].L);
					}
					Delta.Normalize();
					rolex_domination.Pause();
					if (Delta.L.empty()) continue;
					log->extended_count++;
					log->enumerated_count += Delta.L.size();
					(*log->count_by_length)[k] += Delta.L.size();
					
					// Extend.
					rolex_extension.Resume();
					double Ttimel = Delta.L.front().l;
					for (pair<Vertex, int> w_Sw: NG.NGArc[v][S])
					{
						int w = w_Sw.first;
						int Sw = w_Sw.second;

						// Feasibility check.
						if (k < vrp.prec_count[w]) continue;
						if (n-k-1 < vrp.suc_count[w]) continue;
						if (w != NG.L[r + 1] && vrp.suc_count[NG.L[r + 1]] > n - k - 2) continue;
						if (!NG.V[r].test(w)) continue;
						double a_w = NG.tw[k+1][w].left;
						double b_w = NG.tw[k+1][w].right;
						double Ttimew = max(a_w, vrp.ArrivalTime({v, w}, Ttimel));
						if (epsilon_bigger(Ttimew, b_w)) continue;
						vector<TIState::Piece> EXT;
						int j = 0;
						for (auto& p: Delta.L)
						{
							if (epsilon_bigger(p.l, vrp.tau[v][w].Domain().right)) continue;
							double tt = INFTY;
							double ll = INFTY;
							while (j < vrp.tau[v][w].PieceCount())
							{
								if (epsilon_smaller(vrp.tau[v][w][j].domain.right, p.l)) { ++j; continue; }
								if (epsilon_bigger(vrp.tau[v][w][j].domain.left, p.r)) break;
								if (vrp.tau[v][w][j].domain.Includes(p.l)) ll = p.l + vrp.tau[v][w][j](p.l);
								if (epsilon_smaller_equal(p.l, vrp.tau[v][w][j].domain.left) && epsilon_bigger_equal(p.r, vrp.tau[v][w][j].domain.right))
								{
									tt = min(tt, vrp.tau[v][w][j].image.left);
								}
								else
								{
									double tl = max(vrp.tau[v][w][j].domain.left, p.l);
									double tr = min(vrp.tau[v][w][j].domain.right, p.r);
									tt = min(tt, vrp.tau[v][w][j](tl));
									tt = min(tt, vrp.tau[v][w][j](tr));
								}
								if (epsilon_bigger_equal(vrp.tau[v][w][j].domain.right, p.r)) break;
								++j;
							}
							if (epsilon_bigger(ll, b_w)) break;
							double rr = min(b_w, p.r + tt);
							double cost = p.cost + tt - lambda[w];
							if (!EXT.empty() && epsilon_smaller_equal(EXT.back().cost, cost) && epsilon_bigger_equal(EXT.back().r, rr)) continue;
							if (w == vrp.d)
							{
								if (cost < best_cost)
								{
									best_cost = min(best_cost, cost);
									best_piece = TIState::Piece(ll, rr, cost, &p, w);
								}
								continue;
							}
							if (!EXT.empty() && epsilon_equal(EXT.back().l, ll) && epsilon_bigger_equal(EXT.back().cost, cost)) EXT.pop_back();
							EXT.emplace_back(TIState::Piece(ll, rr, cost, &p, w));
						}

						if (!EXT.empty())
						{
							rolex_extension.Pause();
							rolex_domination.Resume();
							D[k + 1][r + (NG.L[r + 1] == w)][w][Sw].Merge(EXT);
							rolex_domination.Pause();
							rolex_extension.Resume();
						}
					}
					rolex_extension.Pause();
				}
			}
		}
	}

	// Rebuild path.
	GraphPath P;
	TIState::Piece* p = &best_piece;
	for (int k = 0; k < n; ++k)
	{
		P.push_back(p->v);
		p = p->prev;
	}
	P = reverse(P);

	// Calculate path departure.
	double P_lambda = sum<Vertex>(P, [&](Vertex v) { return lambda[v]; });
	double P_duration = best_cost + P_lambda;
	if (epsilon_bigger(P_duration - P_lambda, best_cost)) fail("Wrong cost");
	best_cost = P_duration - P_lambda;
	best_route = Route(P, -1, P_duration);

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

struct Merger
{
	vector<double> lambda;
	GraphPath L;
	Matrix<vector<pair<VertexSet, TIState>>> F;
	
	Merger(int n, const vector<double>& lambda, const GraphPath& L) : lambda(lambda), L(L)
	{
		F = Matrix<vector<pair<VertexSet, TIState>>>(n, n);
	}
	
	void AddForward(const VertexSet& S, TIState& s, Vertex v, int r)
	{
		F[r][v].push_back({S, s});
	}
	
	void MergeBackward(const VertexSet& S, TIState& s, Vertex v, int r, double& min_cost, GraphPath& min_path)
	{
		if (s.L.empty()) return;
		int F_r = L.size() - 1 - r;
		if (L[F_r] != v) F_r--;
		vector<LinearFunction> Delta2;
		for (int i = s.L.size()-1; i >= 0; --i) Delta2.push_back(LinearFunction({-s.L[i].r, s.L[i].cost+lambda[v]}, {-s.L[i].l, s.L[i].cost+lambda[v]}));
		for (auto& S_s: F[F_r][v])
		{
			auto& Sf = S_s.first;
			auto& Delta = S_s.second;
			if ((Sf & S) != create_bitset<MAX_N>({v})) continue;
			int i = 0, j = 0;
			while (i < Delta.L.size() && j < Delta2.size())
			{
				auto& p_i = Delta.L[i];
				auto& p_j = Delta2[j];
				
				if (epsilon_bigger(p_i.l, max(dom(p_j))))
				{
					// Case 1: p_i can not be extended with p_j, skip p_j.
					++j;
				}
				else if (epsilon_bigger(min(dom(p_j)), p_i.r))
				{
					// Case 2: p_j comes after p_i then p_i is bounded by waiting.
					double ti = p_i.cost, tj = p_j(min(dom(p_j)));
					double wait = min(dom(p_j)) - p_i.r;
					if (ti+wait+tj < min_cost)
					{
						min_cost = ti + wait + tj;
						min_path = Delta.L[i].Path();
						auto Pb = s.L[s.L.size()-j-1].Path();
						Pb.pop_back();
						Pb = reverse(Pb);
						min_path.insert(min_path.end(), Pb.begin(), Pb.end());
					}
					++i;
				}
				else
				{
					// Case 3: p_i overlaps p_j, we must get the best bound.
					double t0 = max(p_i.l, min(dom(p_j)));
					double t1 = min(p_i.r, max(dom(p_j)));
					if (min_cost > min(p_i.cost+p_j(t0), p_i.cost+p_j(t1)))
					{
						min_cost = min(p_i.cost+p_j(t0), p_i.cost+p_j(t1));
						min_path = Delta.L[i].Path();
						auto Pb = s.L[s.L.size()-j-1].Path();
						Pb.pop_back();
						Pb = reverse(Pb);
						min_path.insert(min_path.end(), Pb.begin(), Pb.end());
					}
					if (epsilon_smaller(p_i.r, max(dom(p_j)))) ++i;
					else ++j;
				}
			}
		}
	}
};

double bidirectional_run_nglti(const VRPInstance& fvrp, const VRPInstance& bvrp, const NGStructure& fNG, const NGStructure& bNG, const vector<double>& lambda, double UB,
				 Route& best_route, double& best_cost, BLBExecutionLog* blb_log)
{
	best_route = Route({}, 0.0, INFTY);
	best_cost = INFTY;
	Stopwatch rolex_blb(true);
	int n = fvrp.D.VertexCount();
	int R = fNG.L.size();
	auto& N = fNG.N;
	Matrix<vector<vector<TIState>>> DS[2] = {
		Matrix<vector<vector<TIState>>>(n, R, vector<vector<TIState>>(n, vector<TIState>(8))),
		Matrix<vector<vector<TIState>>>(n, R, vector<vector<TIState>>(n, vector<TIState>(8)))
	};
	int nn[] = {n / 2+1, n - n/2};
	const NGStructure* NGS[] = {&fNG, &bNG};
	const VRPInstance* vrps[] = {&fvrp, &bvrp};
	for (int d: {0, 1})
	{
		MLBExecutionLog log(true);
		auto& vrp = *vrps[d];
		auto& NG = *NGS[d];
		auto& L = NG.L;
		auto& V = NG.V;
		auto& D = DS[d];
		
		stretch_to_size(*log.count_by_length, n + 1, 0);
		Stopwatch rolex(true), rolex_domination(false), rolex_extension(false), rolex_queuing(false);
		
		// Initialize queue.
		vector<TIState::Piece> initial_state = {
			TIState::Piece(vrp.a[vrp.o], vrp.b[vrp.o], -lambda[vrp.o], nullptr, vrp.o)};
		D[1][0][vrp.o][0].Merge(initial_state);
		for (int k = 1; k <= nn[d]; ++k)
		{
			for (int r = 0; r < R - 1; ++r)
			{
				for (int v = 0; v < n; ++v)
				{
					for (int S = 0; S < NG.NGSet[v].size(); ++S)
					{
						auto& Delta = D[k][r][v][S];
						if (Delta.L.empty()) continue;
						log.processed_count++;
						
						// Dominate by subsets.
						rolex_domination.Resume();
						for (int sub: NG.NGSub[v][S])
						{
							if (Delta.L.empty()) break;
							Delta.DominateBy(D[k][r][v][sub].L);
						}
						Delta.Normalize();
						rolex_domination.Pause();
						if (Delta.L.empty()) continue;
						log.extended_count++;
						log.enumerated_count += Delta.L.size();
						(*log.count_by_length)[k] += Delta.L.size();
						
						// Extend.
						rolex_extension.Resume();
						double Ttimel = Delta.L.front().l;
						for (pair<Vertex, int> w_Sw: NG.NGArc[v][S])
						{
							int w = w_Sw.first;
							int Sw = w_Sw.second;
							
							// Feasibility check.
							if (k < vrp.prec_count[w]) continue;
							if (n - k - 1 < vrp.suc_count[w]) continue;
							if (w != NG.L[r + 1] && vrp.suc_count[NG.L[r + 1]] > n - k - 2) continue;
							if (!NG.V[r].test(w)) continue;
							double Ttimew = vrp.ArrivalTime({v, w}, Ttimel);
							double b_w = NG.tw[k + 1][w].right;
							if (epsilon_bigger(Ttimew, b_w)) continue;
							vector<TIState::Piece> EXT;
							int j = 0;
							for (auto& p: Delta.L)
							{
								double tt = INFTY;
								while (j < vrp.tau[v][w].PieceCount())
								{
									if (epsilon_bigger(vrp.tau[v][w][j].domain.left, p.r)) break;
									if (vrp.tau[v][w][j].domain.Intersects({p.l, p.r}))
										tt = min(tt, vrp.tau[v][w][j].image.left);
									if (epsilon_bigger_equal(vrp.tau[v][w][j].domain.right, p.r)) break;
									++j;
								}
								if (epsilon_bigger(p.l + tt, b_w)) break;
								double ll = max(p.l + tt, Ttimew);
								double rr = min(max(p.r + tt, Ttimew), b_w);
								if (!EXT.empty())
								{
									ll = max(ll, EXT.back().r);
									rr = max(rr, EXT.back().r);
								}
								tt = max(tt, ll - p.r);
								double cost =
									p.cost + (epsilon_equal(rr, Ttimew) ? Ttimew - min(p.r, b_w - tt) : tt) - lambda[w];
								if (!EXT.empty() && epsilon_smaller_equal(EXT.back().cost, cost) &&
									epsilon_bigger_equal(EXT.back().r, rr))
									continue;
								if (!EXT.empty() && epsilon_equal(EXT.back().l, ll) &&
									epsilon_bigger_equal(EXT.back().cost, cost))
									EXT.pop_back();
								EXT.emplace_back(TIState::Piece(ll, rr, cost, &p, w));
							}
							
							if (!EXT.empty())
							{
								rolex_extension.Pause();
								rolex_domination.Resume();
								D[k + 1][r + (NG.L[r + 1] == w)][w][Sw].Merge(EXT);
								rolex_domination.Pause();
								rolex_extension.Resume();
							}
						}
						rolex_extension.Pause();
					}
				}
			}
		}
		
		log.queuing_time = rolex_queuing.Peek();
		log.domination_time = rolex_domination.Peek();
		log.extension_time = rolex_extension.Peek();
		log.time = rolex.Peek();
		if (d == 0) blb_log->forward_log = log;
		else blb_log->backward_log = log;
	}
	
	// Merge.
	Stopwatch rolex_merge(true);
	Merger M(n, lambda, fNG.L);
	double min_cost = INFTY;
	for (int r = 0; r < R-1; ++r)
		for (Vertex v: fvrp.D.Vertices())
			for (int S = 0; S < DS[0][nn[0]][r][v].size(); ++S)
				M.AddForward(fNG.NGSet[v][S], DS[0][nn[0]][r][v][S], v, r);
	
	GraphPath min_path;
	for (int r = 0; r < R-1; ++r)
		for (Vertex v: fvrp.D.Vertices())
			for (int S = 0; S < DS[1][nn[1]][r][v].size(); ++S)
				M.MergeBackward(bNG.NGSet[v][S], DS[1][nn[1]][r][v][S], v, r, min_cost, min_path);
	
	VertexSet NGSet;
	for (auto& v: min_path)
	{
		NGSet = (NGSet & fNG.N[v]);
		NGSet.set(v);
	}
	blb_log->merge_time = rolex_merge.Peek();
	blb_log->time = rolex_blb.Peek();
	best_cost = min_cost;
	best_route = Route(min_path, 0.0, min_cost + sum<Vertex>(min_path, [&](Vertex v) { return lambda[v]; }));
	return min_cost + sum(lambda);
}
} // namespace tdtsptw