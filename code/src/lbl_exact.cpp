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

struct PWLState
{
	void Add(const LinearFunction& f)
	{
		LinearFunction p = f;
		vector<LinearFunction> ND;
		
		// If p is a point.
		if (epsilon_equal(min(dom(p)), max(dom(p))))
		{
			double t = min(dom(p)), v = min(img(p));
			bool dominated = false;
			for (auto& q: F)
			{
				if (epsilon_bigger(min(dom(q)), t)) break;
				// Check if q includes p domain, then check if it is smaller equal.
				if (dom(q).Includes(t)) dominated = epsilon_smaller_equal(q(t), v);
				// Check if waiting after q dominates t.
				else dominated = epsilon_smaller_equal(q(max(dom(q)))+t-max(dom(q)), v);
				if (dominated) break;
			}
			if (!dominated) ND.push_back(p);
		}
		// Otherwise, p is not a point.
		else
		{
			for (auto& q: F)
			{
				// If p is a point, then it is dominated.
				if (epsilon_equal(min(dom(p)), max(dom(p)))) break;
				if (epsilon_bigger(min(dom(q)), max(dom(p)))) break;
				
				for (int i: {0, 1})
				{
					// If p is a point, then it is dominated.
					if (epsilon_equal(min(dom(p)), max(dom(p)))) break;
					
					LinearFunction q_i = i == 0 ? q : LinearFunction({max(dom(q)), q(max(dom(q)))}, {max(dom(p)), q(max(dom(q)))+max(dom(p))-max(dom(q))});
					if (epsilon_smaller(max(dom(q_i)), min(dom(p)))) continue;
					
					// We have intersection between p and q.
					double l = max(min(dom(p)), min(dom(q_i)));
					double r = min(max(dom(p)), max(dom(q_i)));
					double m = p.Intersection(q_i);
					
					double pl = p(l), pr = p(r);
					double ql = q_i(l), qr = q_i(r);
					
					// Leave in [l, r] the dominated portion.
					if (epsilon_smaller_equal(ql, pl) && epsilon_smaller_equal(qr, pr))
					{
						// [l, r] is dominated.
					}
					else if (epsilon_bigger_equal(ql, pl) && epsilon_bigger_equal(qr, pr))
					{
						// no domination.
						l = r + 1;
					}
					else if (epsilon_bigger(m, l) && epsilon_smaller(m, r))
					{
						// [l, m] is dominated, but [m, r] is dominated for q.
						if (epsilon_bigger_equal(pl, ql))
						{
							if (i == 0) q = q.RestrictDomain({min(dom(q)), m});
							r = m;
						}
						// [m, r] is dominated, but [l, m] is dominated for q.
						else
						{
							if (i == 0) q = q.RestrictDomain({m, max(dom(q))});
							l = m;
						}
					}
					
					// If some domination occured...
					if (l != r + 1)
					{
						// A prefix is dominated.
						if (epsilon_equal(l, min(dom(p)))) p.domain.left = r;
						// A suffix is dominated.
						else if (epsilon_equal(r, max(dom(p)))) p.domain.right = l;
						// Middle is dominated, therefore the left side will not be dominated.
						else
						{
							ND.push_back(LinearFunction({min(dom(p)), p(min(dom(p)))}, {l, p(l)}));
							p.domain.left = r;
						}
					}
				}
			}
			if (epsilon_different(min(dom(p)), max(dom(p)))) ND.push_back(p);
		}
		
//		// Dominate F with waiting time in p.
//		if (!ND.empty() && max(dom(ND.back())) == max(dom(p)) && !F.empty())
//		{
//			double t = max(dom(p)), v = p(max(dom(p)));
//			LinearFunction p_w({t, v}, {max(dom(F.back())), v+max(dom(F.back()))-t});
//			for (int i = (int)F.size() - 1; i >= 0; --i)
//			{
//				auto& q = F[i];
//				if (epsilon_smaller(max(dom(q)), t)) break;
//				if (epsilon_bigger_equal(min(dom(q)), t))
//				{
//					if (epsilon_bigger_equal(q(min(dom(q))), p_w(min(dom(q)))))
//					{
//
//					}
//				}
//				if (epsilon_bigger_equal(F[i](max(dom(F[i])))))
//			}
//		}
		
		// Now we have in ND the non dominated pieces of p, we have to insert them into F.
		for (auto& r: ND)
		{
			F.push_back(r);
			LB.push_back(INT_MAX);
			
			int i = (int)F.size()-1;
			while (i > 0 && epsilon_smaller(min(dom(F[i])), min(dom(F[i-1]))))
			{
				swap(F[i], F[i-1]);
				swap(LB[i], LB[i-1]);
				--i;
			}
		}
	}
	
	LinearFunction& operator[](int i)
	{
		return F[i];
	}
	
	int& lb(int i)
	{
		return LB[i];
	}
	
	void Show()
	{
		clog << F << endl;
	}
	
	vector<LinearFunction> F;
	vector<int> LB;
};

Route run_exact_piecewise(const VRPInstance& vrp, const GraphPath& L, const vector<double>& lambda,
	double LB, double UB, MLBExecutionLog* log)
{
	PWLState S;
	S.Add(LinearFunction({0, 10}, {0, 100}));
	S.Add(LinearFunction({4, 20}, {6, 20}));
	S.Show();
	exit(0);
	
	int n = vrp.D.VertexCount();
	int BASE = floor(LB), TOP = ceil(UB);
	
	// Init structures.
	auto Q = vector<vector<vector<vector<unordered_map<VertexSet, vector<LinearFunction>>>>>>(TOP-BASE+1, vector<vector<vector<unordered_map<VertexSet, vector<LinearFunction>>>>>(n, vector<vector<unordered_map<VertexSet, vector<LinearFunction>>>>(L.size(), vector<unordered_map<VertexSet, vector<LinearFunction>>>(n))));
	auto D = vector<vector<vector<unordered_map<VertexSet, MinFunc>>>>(n, vector<vector<unordered_map<VertexSet, MinFunc>>>(L.size(), vector<unordered_map<VertexSet, MinFunc>>(n)));
	vector<LinearFunction> R;
	
	Q[0][1][0][0].insert({create_bitset<MAX_N>({vrp.o}), {}}).first->second.push_back({{vrp.a[vrp.o], -lambda[vrp.o]}, {vrp.b[vrp.o], -lambda[vrp.o]}});
	for (int lb = 0; lb <= ceil(UB)-BASE; ++lb)
	{
		for (int k = 1; k < n; ++k)
		{
			clog << k << endl;
			for (int r = 0; r < (int)L.size()-1; ++r)
			{
				for (int v = 0; v < n; ++v)
				{
					for (auto& S_f: Q[lb][k][r][v])
					{
						// Domination.
						auto& S = S_f.first;
						auto& f = S_f.second;
						auto& D_S = D[k][r][v].insert({S, MinFunc()}).first->second;
						auto P = D_S.Merge(f);
						if (P.size() > 30)
						{ clog << P << endl; exit(0);}
						clog << P.size() << endl;
						if (P.empty()) continue; // All pieces were dominated by previous pieces.
						
						// Extension.
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
							
							for (auto& p_i: P)
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
										Q[lb_i][k + 1][r + (w == L[r + 1])][w].insert({S_w, {}}).first->second.push_back(pp);
									}
									
									// Stop if tau_vw exceeds the latest departure time.
									if (epsilon_bigger_equal(max(dom(tau_j)), LDTw_at_v)) j = vrp.tau[v][w].PieceCount();
								}
							}
						}
					}
				}
			}
		}
		if (!R.empty()) break;
	}
	
	// Rebuild solution.
	double best_t = INFTY;
	for (auto& p: R) best_t = min(best_t, p.image.left);
	return Route({}, 0.0, best_t);
}
} // namespace tdtsptw