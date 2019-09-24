//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "lbl_exact.h"

#include "pwl_domination_function.h"

using namespace std;
using namespace goc;

namespace tdtsptw
{
Route run_exact(const VRPInstance& vrp, const NGStructure& NG, BoundingStructure& B, const vector<double>& lambda,
				const Route& UB, double lb, MLBExecutionLog* log)
{
	Stopwatch rolex(true), rolex_domination(false), rolex_queuing(false), rolex_extension(false), rolex_bounding(false);
	Route best = UB;
	int n = vrp.D.VertexCount();
	stretch_to_size(*log->count_by_length, n+1, 0);
	double Lambda = sum(lambda);
	TimeUnit T = vrp.T;
	
	// q is the queue where q[LB][k][v] is the bucket of labels with duration bound LB, having k vertices and
	// ending at vertex v.
	Matrix<vector<vector<Label*>>> q(UB.duration+1, n, vector<vector<Label*>>(n, vector<Label*>()));
	q[0][1][vrp.o].push_back(new Label(nullptr, vrp.o, create_bitset<MAX_N>({vrp.o}), PWLFunction::ConstantFunction(0.0, vrp.tw[vrp.o]), vrp.a[vrp.o], lambda[vrp.o]));
	B.UB = best.duration;
	// Dominance structure: D[v][S] contains labels with v(l)=v, S(l)=S sorted by cost.
	vector<unordered_map<VertexSet, vector<Label*>>> D(n);
	for (int LB = 0; LB < best.duration; ++LB)
	{
		for (int k = 1; k < n; ++k)
		{
			for (Vertex v: vrp.D.Vertices())
			{
				// Check if LB has reached best.
				if (epsilon_bigger_equal(LB, best.duration)) break;
				
				rolex_queuing.Resume();
				sort(q[LB][k][v].begin(), q[LB][k][v].end(), [](Label* l1, Label* l2) { return make_tuple(l1->Ttime, min(img(l1->Tdur))) < make_tuple(l2->Ttime, min(img(l2->Tdur))); });
				rolex_queuing.Pause();
				
				for (Label* l: q[LB][k][v])
				{
					// Check if LB has reached best.
					if (epsilon_bigger_equal(LB, best.duration)) break;
					
					// Domination.
					rolex_domination.Resume();
					bool is_dominated = false;
					PWLDominationFunction l_D = l->Tdur;
					for (Label* m: D[l->v][l->S])
					{
						if (epsilon_bigger(min(img(m->Tdur)), max(img(l->Tdur)))) break;
						if (!l_D.DominatePieces(m->Tdur)) continue;
						is_dominated = true;
						break;
					}
					rolex_domination.Pause();
					if (is_dominated)
					{
						log->dominated_count++;
						rolex_domination.Pause();
						delete l;
						continue;
					}
					
					// Get non dominated pieces.
					l->Tdur = (PWLFunction)l_D;
					l->Ttime = min(dom(l->Tdur));
					insert_sorted(D[l->v][l->S], l, [] (Label* l1, Label* l2) { return min(img(l1->Tdur)) < min(img(l2->Tdur));});
					rolex_domination.Pause();
					
					// Extension.
					(*log->count_by_length)[k]++;
					log->processed_count++;
					for (Vertex w: vrp.D.Successors(v))
					{
						// Feasibility check.
						if (l->S.test(w)) continue;
						rolex_extension.Resume();
						VertexSet Sw = l->S;
						Sw.set(w);
						double Ttimew = vrp.ArrivalTime({v,w}, l->Ttime);
						double LDTw = vrp.b[w];
						for (Vertex u: vrp.D.Vertices())
							if (!Sw.test(u))
								LDTw = min(LDTw, vrp.LDT[w][u]);
						if (epsilon_smaller(LDTw, Ttimew)) continue;
						
						// Extension.
						PWLFunction Tdurw = epsilon_smaller(max(dom(l->Tdur)), min(img(vrp.dep[v][w])))
										  ? PWLFunction::ConstantFunction(l->Tdur(max(dom(l->Tdur))) + vrp.a[w] - max(dom(l->Tdur)), {vrp.a[w], vrp.a[w]})
										  : (l->Tdur + vrp.tau[v][w]).Compose(vrp.dep[v][w]);
						
						// Optimization: We do not need points after LDTw.
						Tdurw = Tdurw.RestrictDomain({min(dom(Tdurw)), LDTw});
						if (Tdurw.Empty()) continue;
						auto lw = new Label(l, w, Sw, Tdurw, Ttimew, l->lambda + lambda[w]);
						
						// Process tour.
						if (w == vrp.d)
						{
							if (min(img(lw->Tdur)) < best.duration) best = Route(lw->Path(), max(dom(lw->Tdur))-lw->Tdur(max(dom(lw->Tdur))), min(img(lw->Tdur)));
							break;
						}
						log->enumerated_count++;
						rolex_extension.Pause();
						
						// Domination.
						rolex_domination.Resume();
						bool is_dominated = false;
						PWLDominationFunction lw_D = lw->Tdur;
						for (Label* m: D[lw->v][lw->S])
						{
							if (epsilon_bigger(min(img(m->Tdur)), max(img(lw->Tdur)))) break;
							if (!lw_D.DominatePieces(m->Tdur)) continue;
							is_dominated = true;
							break;
						}
						rolex_domination.Pause();
						if (is_dominated)
						{
							log->dominated_count++;
							rolex_domination.Pause();
							continue;
						}
						
						// Get non dominated pieces.
						lw->Tdur = (PWLFunction)lw_D;
						lw->Ttime = min(dom(lw->Tdur));
						rolex_domination.Pause();
						
						// Bounding.
						rolex_bounding.Resume();
						double LBw = B.CompletionBound(k+1, -1, *lw);
						if (epsilon_bigger_equal(LBw, best.duration))
						{
							log->bounded_count++;
							rolex_bounding.Pause();
							continue;
						}
						LBw = max(LB, (int)floor(LBw));
						q[LBw][k+1][w].push_back(lw);
						rolex_bounding.Pause();
						log->extended_count++;
					}
					rolex_extension.Pause();
				}
				
				q[LB][k][v].clear();
			}
		}
	}
	log->domination_time = rolex_domination.Peek();
	log->extension_time = rolex_extension.Peek();
	log->queuing_time = rolex_queuing.Peek();
	log->bounding_time = rolex_bounding.Peek();
	log->time = rolex.Peek();
	
	return best;
}

// ------------------ Labeling by pieces ----------------------.
// ------------------ Labeling by pieces ----------------------.
// ------------------ Labeling by pieces ----------------------.
// ------------------ Labeling by pieces ----------------------.
// ------------------ Labeling by pieces ----------------------.


State::Piece::Piece(double lb, const LinearFunction& f)
	: lb(lb), f(f)
{ }

bool State::Piece::Dominate(Piece& p2)
{
	auto& p = f;
	auto& q = p2.f;
	// Remove the case when p is after q.
	if (epsilon_bigger(min(dom(p)), max(dom(q)))) return false;

	double t1 = -INFTY, t2 = INFTY;
	if (epsilon_equal(p.slope, q.slope))
	{
		if (epsilon_bigger(p.intercept, q.intercept))
			return false;
	}
	else if (epsilon_smaller(p.slope, q.slope))
	{
		t1 = p.Intersection(q);
	}
	else
	{
		t2 = p.Intersection(q);
	}

	// Dominate by waiting time.
	if (epsilon_bigger(max(dom(q)), max(dom(p))) && epsilon_bigger(t2, max(dom(p))))
	{
		LinearFunction ID({max(dom(p)), p(max(dom(p)))}, {max(dom(q)), max(dom(q))-max(dom(p))+p(max(dom(p)))});
		t2 = ID.Intersection(q);
	}
	t1 = max(t1, min(dom(p)));
	t2 = min(t2, max(dom(q)));

	// q is dominated in [t1, t2].
	// If all q is dominated, return true.
	if (epsilon_smaller_equal(t1, min(dom(q))) && epsilon_bigger_equal(t2, max(dom(q))))
	{
		return true;
	}
	// A prefix is dominated.
	else if (epsilon_smaller_equal(t1, min(dom(q))))
	{
		q = q.RestrictDomain({t2+EPS, max(dom(q))});
	}
	// A suffix is dominated
	else if (epsilon_bigger_equal(t2, max(dom(q))))
	{
		q = q.RestrictDomain({min(dom(q)), t1-EPS});
	}
	return false;
}

void State::Piece::Print(ostream& os) const
{
	os << f;
}

void State::Merge(vector<Piece>& P)
{
	auto P1 = F;
	auto& P2 = P;
	F.clear();

	// Sweep pieces.
	int i = 0, j = 0;
	for (; i < P1.size() && j < P2.size();)
	{
		if (P1[i].Dominate(P2[j])) { ++j; continue; }
		if (P2[j].Dominate(P1[i])) { ++i; continue; }

		if (epsilon_smaller(min(dom(P1[i].f)), min(dom(P2[j].f))))
		{
			if (epsilon_smaller_equal(max(dom(P1[i].f)), min(dom(P2[j].f))))
			{
				F.push_back(P1[i]);
				++i;
			}
			else
			{
				F.push_back(Piece(P1[i].lb, P1[i].f.RestrictDomain({min(dom(P1[i].f)), min(dom(P2[j].f))})));
				P1[i].f.domain.left = P2[j].f.domain.left;
			}
		}
		else
		{
			if (epsilon_smaller_equal(max(dom(P2[j].f)), min(dom(P1[i].f))))
			{
				F.push_back(P2[j]);
				++j;
			}
			else
			{
				F.push_back(Piece(P2[j].lb, P2[j].f.RestrictDomain({min(dom(P2[j].f)), min(dom(P1[i].f))})));
				P2[j].f.domain.left = P1[i].f.domain.left;
			}
		}
	}
	// Add remaining pieces.
	F.insert(F.end(), P1.begin()+i, P1.end());
	F.insert(F.end(), P2.begin()+j, P2.end());

	// Checks.
	// Assert is increasing in domain.
	for (int i = 0; i < (int)F.size()-1; ++i)
	{
		assert(epsilon_smaller_equal(max(dom(F[i].f)), min(dom(F[i+1].f))));
	}
}
	
void State::DominateBy(const State& s2)
{
	auto P1 = s2.F;
	auto P2 = F;
	F.clear();
	
	// Sweep pieces.
	int i = 0, j = 0;
	for (; i < P1.size() && j < P2.size();)
	{
		if (P1[i].Dominate(P2[j])) { ++j; continue; }
		if (P2[j].Dominate(P1[i])) { ++i; continue; }
		
		if (epsilon_smaller(min(dom(P1[i].f)), min(dom(P2[j].f))))
		{
			if (epsilon_smaller_equal(max(dom(P1[i].f)), min(dom(P2[j].f))))
			{
				++i;
			}
			else
			{
				P1[i].f.domain.left = P2[j].f.domain.left;
			}
		}
		else
		{
			if (epsilon_smaller_equal(max(dom(P2[j].f)), min(dom(P1[i].f))))
			{
				F.push_back(P2[j]);
				++j;
			}
			else
			{
				F.push_back(Piece(P2[j].lb, P2[j].f.RestrictDomain({min(dom(P2[j].f)), min(dom(P1[i].f))})));
				P2[j].f.domain.left = P1[i].f.domain.left;
			}
		}
	}
	// Add remaining pieces.
	F.insert(F.end(), P2.begin()+j, P2.end());
	
	// Checks.
	// Assert is increasing in domain.
	for (int i = 0; i < (int)F.size()-1; ++i)
	{
		assert(epsilon_smaller_equal(max(dom(F[i].f)), min(dom(F[i+1].f))));
	}
}

Bounding::Bounding(const VRPInstance& vrp, const NGStructure& NG, const std::vector<double>& lambda)
	: vrp(vrp), NG(NG), lambda(lambda)
{
	int n = vrp.D.VertexCount();
	B = vector<vector<vector<vector<pair<VertexSet, vector<LinearFunction>>>>>>(n+1, vector<vector<vector<pair<VertexSet, vector<LinearFunction>>>>>(NG.L.size()+1, vector<vector<pair<VertexSet, vector<LinearFunction>>>>(n)));
	for (auto& v: NG.L) LSet.set(v);
	Lambda = sum(lambda);
}

void Bounding::AddBound(int k, int r, goc::Vertex v, const VertexSet& S, const State& Delta)
{
	// Reverse pieces in Delta.
	vector<LinearFunction> DeltaR;
	for (int j = (int)Delta.F.size()-1; j >= 0; --j)
	{
		auto& p = Delta.F[j].f;
		DeltaR.push_back({{-max(dom(p)), p(max(dom(p)))+Lambda+lambda[v]}, {-min(dom(p)), p(min(dom(p)))+Lambda+lambda[v]}});
	}
	
	int n = vrp.D.VertexCount();
	int R = NG.L.size();
//	clog << n-k+1 << "\t" << R-r-(!LSet.test(v)) << "\t" << v << "\t" << S << "\t" << DeltaR << endl;
	B[n-k+1][R-r-(!LSet.test(v))][v].push_back({S, DeltaR});
}

void Bounding::Bound(goc::Vertex v, VertexSet S, State& Delta)
{
//	clog << S.count() << "\t" << v << "\t" << S << "\t" << Delta.F << endl;
	int k = S.count();
	int r = (LSet & S).count();
	
	bool all_bounded = all_of(Delta.F.begin(), Delta.F.end(), [] (State::Piece& p) { return p.lb > -1; });
	if (all_bounded) return;
	
	for (auto& p: Delta.F) if (p.lb == -1) p.lb = INFTY;
	
	for (auto& S2_Delta2: B[k][r][v])
	{
		auto& S2 = S2_Delta2.first;
		auto& Delta2 = S2_Delta2.second;
		if ((S2 & S) != create_bitset<MAX_N>({v})) continue;
		
		// Bound.
		int i = 0, j = 0;
		while (i < Delta.F.size() && j < Delta2.size())
		{
			auto& p_i = Delta.F[i].f;
			auto& p_j = Delta2[j];
			
			if (epsilon_bigger(min(dom(p_i)), max(dom(p_j))))
			{
				// Case 1: p_i can not be extended with p_j, skip p_j.
				++j;
			}
			else if (epsilon_bigger(min(dom(p_j)), max(dom(p_i))))
			{
				// Case 2: p_j comes after p_i.
				// Case 2a: if we do not have p_{i+1} which comes before p_j, then p_i is bounded by waiting.
				if (i == Delta.F.size()-1 || epsilon_bigger(min(dom(Delta.F[i+1].f)), min(dom(p_j))))
				{
					double ti = p_i(max(dom(p_i))), tj = p_j(min(dom(p_j)));
					double wait = min(dom(p_j)) - max(dom(p_i));
					Delta.F[i].lb = min(Delta.F[i].lb, ti+wait+tj);
				}
				++i;
			}
			else
			{
				// Case 3: p_i overlaps p_j, we must get the best bound.
				double t0 = max(min(dom(p_i)), min(dom(p_j)));
				double t1 = min(max(dom(p_i)), max(dom(p_j)));
				Delta.F[i].lb = min(Delta.F[i].lb, min(p_i(t0)+p_j(t0), p_i(t1)+p_j(t1)));
				if (epsilon_smaller(max(dom(p_i)), max(dom(p_j)))) ++i;
				else ++j;
			}
		}
	}
}

double run_ngl(const VRPInstance& vrp, const NGStructure& NG, const vector<double>& lambda, MLBExecutionLog* log, Bounding* B)
{
	int n = vrp.D.VertexCount();
	auto& N = NG.N;
	auto& L = NG.L;
	auto& V = NG.V;
	
	// Init structures.
	auto D = vector<vector<vector<unordered_map<VertexSet, State>>>>(n, vector<vector<unordered_map<VertexSet, State>>>(L.size(), vector<unordered_map<VertexSet, State>>(n)));
	double LB = INFTY;
	
	Stopwatch rolex(true), rolex_domination(false), rolex_extension(false), rolex_bounding(false);
	State::Piece p0(-1, {{vrp.a[vrp.o], -lambda[vrp.o]}, {vrp.b[vrp.o], -lambda[vrp.o]}});
	vector<State::Piece> P0 = {p0};
	D[1][0][vrp.o].insert({create_bitset<MAX_N>({vrp.o}), State()}).first->second.Merge(P0);
	
	for (int k = 1; k < n; ++k)
	{
		for (int r = 0; r < (int)L.size(); ++r)
		{
			for (int v = 0; v < n; ++v)
			{
				for (auto& S_Delta: D[k][r][v])
				{
					rolex_domination.Resume();
					// Get non dominated pieces (by the same S).
					auto& S = S_Delta.first;
					auto& Delta = S_Delta.second;
					
					// Dominate pieces by subsets of S.
					for (auto& S2_Delta2: D[k][r][v])
					{
						auto& S2 = S2_Delta2.first;
						auto& Delta2 = S2_Delta2.second;
						if (S2 == S || !is_subset(S2, S)) continue;
						Delta.DominateBy(Delta2);
					}
					rolex_domination.Pause();
					
					if (Delta.F.empty()) continue;
					
					// Add bound to the bounding structure.
					B->AddBound(k, r, v, S, Delta);
					
					log->processed_count++;
					log->extended_count += Delta.F.size();
					log->enumerated_count += Delta.F.size();
					
					// Extension of pieces.
					rolex_extension.Resume();
					for (Vertex w: vrp.D.Successors(v))
					{
						if (S.test(w)) continue;
						if (!V[r].test(w)) continue;
						if (vrp.prec_count[w] > k) continue;
						if (vrp.suc_count[w] > n - k - 1) continue;
						if (w != L[r + 1] && vrp.suc_count[L[r + 1]] > n - k - 2) continue;
						
						double LDTw_at_v = vrp.LDT[v][w];
						
						vector<State::Piece> EXT_P; // Extended pieces.
						int j = 0;
						for (auto& p: Delta.F)
						{
							auto& p_i = p.f;
							if (epsilon_bigger(min(dom(p_i)), LDTw_at_v)) break;
							if (epsilon_bigger(min(dom(vrp.tau[v][w][j])), max(dom(p_i)))) break;
							for (; j < vrp.tau[v][w].PieceCount(); ++j)
							{
								auto& tau_j = vrp.tau[v][w][j];
								
								// Check that tau_j \cap p_i \neq \emptyset.
								if (epsilon_smaller(max(dom(tau_j)), min(dom(p_i)))) continue;
								if (epsilon_bigger(min(dom(tau_j)), max(dom(p_i)))) break;
								
								// Find intersection of pieces.
								double t1 = max(min(dom(tau_j)), min(dom(p_i)));
								double t2 = min(LDTw_at_v, min(max(dom(tau_j)), max(dom(p_i))));
								
								// Extend.
								LinearFunction pp({t1 + tau_j(t1), p_i(t1) + tau_j(t1) - lambda[w]},
												  {t2 + tau_j(t2), p_i(t2) + tau_j(t2) - lambda[w]});
								if (k == n - 1) // If complete tour, add it to the solution.
								{
									LB = min(LB, min(img(pp)));
								}
								else // Otherwise, add it to the queue.
								{
									if (!EXT_P.empty() && epsilon_equal(min(dom(EXT_P.back().f)), min(dom(pp))))
										EXT_P.pop_back(); // Remove waiting time piece.
									EXT_P.push_back(State::Piece(-1, pp));
								}
								
								// Stop if tau_vw exceeds the latest departure time.
								if (epsilon_bigger_equal(max(dom(tau_j)), max(dom(p_i)))) break;
							}
						}
						if (!EXT_P.empty())
						{
							rolex_extension.Pause();
							rolex_domination.Resume();
							auto S_w = S & N[w];
							S_w.set(w);
							D[k + 1][r + (w == L[r + 1])][w].insert({S_w, State()}).first->second.Merge(EXT_P);
							rolex_domination.Pause();
							rolex_extension.Resume();
						}
					}
					rolex_extension.Pause();
				}
			}
		}
	}
	
	return LB+sum(lambda);
}

Route run_exact_piecewise(const VRPInstance& vrp, const GraphPath& L, const vector<double>& lambda,
	double LB, double UB, MLBExecutionLog* log, Bounding* B)
{
	int n = vrp.D.VertexCount();
	int BASE = floor(LB), TOP = ceil(UB);
	
	// Init structures.
	auto D = vector<vector<unordered_map<VertexSet, State>>>(n, vector<unordered_map<VertexSet, State>>(n));
	vector<LinearFunction> R;

	Stopwatch rolex(true), rolex_domination(false), rolex_extension(false), rolex_bounding(false);
	State::Piece p0(LB, {{vrp.a[vrp.o], -lambda[vrp.o]}, {vrp.b[vrp.o], -lambda[vrp.o]}});
	vector<State::Piece> P0 = {p0};
	D[1][vrp.o].insert({create_bitset<MAX_N>({vrp.o}), State()}).first->second.Merge(P0);
	for (int lb = 0; lb <= TOP-BASE; ++lb)
	{
		for (int k = 1; k < n; ++k)
		{
			for (int v = 0; v < n; ++v)
			{
				// All labels in D[k][v] are non dominated for lower bound lb.
				for (auto& S_Delta: D[k][v])
				{
					rolex_extension.Resume();
					
					// Get non dominated pieces.
					auto& S = S_Delta.first;
					auto& Delta = S_Delta.second;
					log->processed_count++;
					
					// Apply bounds.
					vector<LinearFunction> EXT; // Pieces with lb = lb that must be extended.
					rolex_extension.Pause();
					rolex_bounding.Resume();
					for (auto& p: Delta.F) if (p.lb == -1) log->enumerated_count++;
					if (B) B->Bound(v, S, Delta);
					for (auto& p: Delta.F)
						if (floor(p.lb) == lb+BASE || !B)
							EXT.push_back(p.f);
					rolex_bounding.Pause();
					rolex_extension.Resume();
					log->extended_count += EXT.size();

					// Extension of pieces.
					for (Vertex w: vrp.D.Successors(v))
					{
						if (S.test(w)) continue;
						double LDT_w = INFTY;
						for (Vertex u: vrp.D.Vertices()) if (u != w && !S.test(u)) LDT_w = min(LDT_w, vrp.LDT[w][u]);
						double LDTw_at_v = vrp.DepartureTime({v,w}, LDT_w);
						if (LDTw_at_v == INFTY) continue;
						auto S_w = S;
						S_w.set(w);

						vector<State::Piece> EXT_P; // Extended pieces.
						int j = 0;
						for (auto& p_i: EXT)
						{
							if (epsilon_bigger(min(dom(p_i)), LDTw_at_v)) break;
							if (epsilon_bigger(min(dom(vrp.tau[v][w][j])), max(dom(p_i)))) break;
							for (; j < vrp.tau[v][w].PieceCount(); ++j)
							{
								auto& tau_j = vrp.tau[v][w][j];
								
								// Check that tau_j \cap p_i \neq \emptyset.
								if (epsilon_smaller(max(dom(tau_j)), min(dom(p_i)))) continue;
								if (epsilon_bigger(min(dom(tau_j)), max(dom(p_i)))) break;
								
								// Find intersection of pieces.
								double t1 = max(min(dom(tau_j)), min(dom(p_i)));
								double t2 = min(LDTw_at_v, min(max(dom(tau_j)), max(dom(p_i))));
								
								// Extend.
								LinearFunction pp({t1+tau_j(t1), p_i(t1)+tau_j(t1)-lambda[w]}, {t2+tau_j(t2), p_i(t2)+tau_j(t2)-lambda[w]});
								if (k == n-1) // If complete tour, add it to the solution.
								{
									R.push_back(pp);
								}
								else // Otherwise, add it to the queue.
								{
									if (!EXT_P.empty() && epsilon_equal(min(dom(EXT_P.back().f)) , min(dom(pp)))) EXT_P.pop_back(); // Remove waiting time piece.
									EXT_P.push_back(State::Piece(-1, pp));
								}
								
								// Stop if tau_vw exceeds the latest departure time.
								if (epsilon_bigger_equal(max(dom(tau_j)), max(dom(p_i)))) break;
							}
						}
						if (!EXT_P.empty())
						{
							rolex_extension.Pause();
							rolex_domination.Resume();
							D[k + 1][w].insert({S_w, State()}).first->second.Merge(EXT_P);
							rolex_domination.Pause();
							rolex_extension.Resume();
						}
					}
					rolex_extension.Pause();
				}
			}
		}
		if (!R.empty()) break;
	}
	log->bounded_count = log->enumerated_count - log->extended_count;
	log->extension_time = rolex_extension.Peek();
	log->bounding_time = rolex_bounding.Peek();
	log->domination_time = rolex_domination.Peek();
	log->time = rolex.Peek();
	
	// Rebuild solution.
	double best_t = INFTY;
	for (auto& p: R) best_t = min(best_t, p.image.left);
	best_t += sum(lambda);
	return Route({}, 0.0, best_t);
}
} // namespace tdtsptw