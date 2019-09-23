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
vector<LinearFunction> MinFunc::Merge(const vector<LinearFunction>& F)
{
	if (F.empty()) return {};
	PWLFunction P;
	for (auto& p: F) P = Min(P, PWLFunction({p}));
	PWLDominationFunction D(P);
	D.DominatePieces(PWLFunction(pieces_));
	if (D.Empty()) return {};
	
	// Now we have in D all non dominated pieces from F, we should merge them with pieces_ and return the
	// new pieces.
	PWLFunction FD(D);
	pieces_ = Min(pieces_, FD);
	return FD.Pieces();
}

const LinearFunction& MinFunc::operator[](int index)
{
	return pieces_[index];
}

struct State
{
public:
	class Piece : public Printable
	{
	public:
		int lb;
		LinearFunction f;

		Piece(int lb, const LinearFunction& f)
			: lb(lb), f(f)
		{ }

		// Reduces domain of p2 if any prefix or suffix is dominated (this(x) <= p2(x)).
		// If p2 is fully dominated the function returns true.
		bool Dominate(Piece& p2)
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

		virtual void Print(ostream& os) const
		{
			os << f;
		}
	};

	// Merges pieces in F with pieces with P, preserving pieces in F when matches occur.
	// Precondition: pieces in P are disjunt and sorted by domain.
	void Merge(vector<Piece>& P)
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

	vector<Piece> F;
};

Route run_exact_piecewise(const VRPInstance& vrp, const GraphPath& L, const vector<double>& lambda,
	double LB, double UB, MLBExecutionLog* log)
{
	int n = vrp.D.VertexCount();
	int BASE = floor(LB), TOP = ceil(UB);
	TOP = floor(LB);
	
	// Init structures.
	auto D = vector<vector<vector<unordered_map<VertexSet, State>>>>(n, vector<vector<unordered_map<VertexSet, State>>>(L.size(), vector<unordered_map<VertexSet, State>>(n)));
	vector<LinearFunction> R;

	State::Piece p0(-1, {{vrp.a[vrp.o], -lambda[vrp.o]}, {vrp.b[vrp.o], -lambda[vrp.o]}});
	vector<State::Piece> P0 = {p0};
	D[1][0][0].insert({create_bitset<MAX_N>({vrp.o}), State()}).first->second.Merge(P0);
	for (int lb = 0; lb <= TOP-BASE; ++lb)
	{
		for (int k = 1; k < n; ++k)
		{
			for (int r = 0; r < (int)L.size()-1; ++r)
			{
				for (int v = 0; v < n; ++v)
				{
					// All labels in D[k][r][v] are non dominated for lower bound lb.
					for (auto& S_f: D[k][r][v])
					{
						log->enumerated_count++;
						// Get non dominated pieces.
						auto& S = S_f.first;
						auto& Delta = S_f.second;
						// Apply bounds.
						vector<LinearFunction> EXT; // Pieces with lb = lb that must be extended.
						for (auto& p: Delta.F)
						{
							if (p.lb == -1)
							{
								p.lb = lb; // TODO: Calculate bound here.
								if (p.lb == lb) EXT.push_back(p.f);
							}
						}
						log->extended_count += Delta.F.size();

						// Extension of pieces.
						for (Vertex w: vrp.D.Successors(v))
						{
							if (S.test(w)) continue;
							double LDT_w = INFTY;
							for (Vertex u: vrp.D.Vertices()) if (u != w && !S.test(u)) LDT_w = min(LDT_w, vrp.LDT[w][u]);
							int j = 0; // index of tau_vw.
							double LDTw_at_v = vrp.DepartureTime({v,w}, LDT_w);
							if (LDTw_at_v == INFTY) continue;
							auto S_w = S;
							S_w.set(w);

							vector<State::Piece> EXT_P; // Extended pieces.
							for (auto& p_i: EXT)
							{
								int lb_i = lb;
								if (epsilon_bigger(min(dom(p_i)), LDTw_at_v)) break;
								for (j = 0; j < vrp.tau[v][w].PieceCount(); ++j)
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
									if (epsilon_bigger_equal(max(dom(tau_j)), LDTw_at_v)) j = vrp.tau[v][w].PieceCount();
								}
							}
							if (!EXT_P.empty())
							{
								D[k + 1][r + (w == L[r + 1])][w].insert({S_w, State()}).first->second.Merge(EXT_P);
							}
						}
					}
				}
			}
		}
		if (!R.empty()) break;
	}

	clog << "#States: " << log->enumerated_count << endl;
	clog << "#Pieces: " << log->extended_count << endl;

	// Rebuild solution.
	double best_t = INFTY;
	for (auto& p: R) best_t = min(best_t, p.image.left);
	return Route({}, 0.0, best_t);
}
} // namespace tdtsptw