//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "dssr.h"
#include "lbl_exact.h"
#include "pwl_domination_function.h"

using namespace std;
using namespace goc;

namespace tdtsptw
{
namespace
{
VRPInstance reverse_instance(const VRPInstance& vrp)
{
	int n = vrp.D.VertexCount();
	
	VRPInstance rev;
	rev.D = vrp.D.Reverse();
	rev.o = vrp.d, rev.d = vrp.o;
	rev.T = vrp.T;
	for (Vertex v: vrp.D.Vertices()) rev.tw.push_back({rev.T-vrp.tw[v].right, rev.T-vrp.tw[v].left});
	for (Vertex v: vrp.D.Vertices()) rev.a.push_back(rev.tw[v].left);
	for (Vertex v: vrp.D.Vertices()) rev.b.push_back(rev.tw[v].right);
	rev.prec = Matrix<bool>(n,n, false);
	rev.prec_count = vector<int>(n, 0);
	rev.suc_count = vector<int>(n, 0);
	for (Vertex v: vrp.D.Vertices())
	{
		for (Vertex w: vrp.D.Vertices())
		{
			if (vrp.prec[w][v])
			{
				rev.prec[v][w] = true;
				rev.prec_count[w]++;
				rev.suc_count[v]++;
			}
		}
	}
	rev.LDT = rev.EAT = Matrix<double>(n, n);
	for (Vertex v: vrp.D.Vertices())
	{
		for (Vertex w: vrp.D.Vertices())
		{
			rev.EAT[v][w] = vrp.T - vrp.LDT[w][v];
			rev.LDT[v][w] = vrp.T - vrp.EAT[w][v];
		}
	}
	rev.arr = rev.tau = rev.dep = rev.pretau = Matrix<PWLFunction>(n, n);
	for (Vertex u: vrp.D.Vertices())
	{
		for (Vertex v: vrp.D.Successors(u))
		{
			// Compute reverse travel functions.
			rev.arr[v][u] = vrp.T - vrp.dep[u][v].Compose(vrp.T - PWLFunction::IdentityFunction({0.0, vrp.T}));
			rev.arr[v][u] = Min(PWLFunction::ConstantFunction(min(img(rev.arr[v][u])), {min(rev.tw[v]), min(dom(rev.arr[v][u]))}), rev.arr[v][u]);
			rev.tau[v][u] = rev.arr[v][u] - PWLFunction::IdentityFunction({0.0, vrp.T});
			rev.dep[v][u] = rev.arr[v][u].Inverse();
			rev.pretau[v][u] = PWLFunction::IdentityFunction(dom(rev.dep[v][u])) - rev.dep[v][u];
		}
	}
	// Add travel functions for (i, i) (for boundary reasons).
	for (Vertex u: rev.D.Vertices())
	{
		rev.tau[u][u] = rev.pretau[u][u] = PWLFunction::ConstantFunction(0.0, rev.tw[u]);
		rev.dep[u][u] = rev.arr[u][u] = PWLFunction::IdentityFunction(rev.tw[u]);
	}
	return rev;
}
}

Label::Label(Label* prev, Vertex v, VertexSet S, const PWLFunction& Tdur, TimeUnit Ttime, double lambda)
	: prev(prev), v(v), S(S), Tdur(Tdur), Ttime(Ttime), lambda(lambda)
{ }

GraphPath Label::Path() const
{
	if (!prev) return {v};
	auto p = prev->Path();
	p.push_back(v);
	return p;
}

void Label::Print(ostream& os) const
{
	os << "{P:" << Path() << ", S: " << S << ", Tdur: " << Tdur << ", lambda: " << lambda << "}";
}

TILabel Label::ToTI() const
{
	return TILabel(nullptr, v, S, min(img(Tdur))-lambda, -(max(dom(Tdur))-Tdur(max(dom(Tdur))))-lambda, lambda);
}

BoundingStructure::BoundingStructure(VRPInstance* vrp, NGStructure* NG, const vector<double>& penalties, double UB)
	: vrp(vrp), NG(NG), penalties(penalties), UB(UB)
{
	int n = vrp->D.VertexCount();
	S = Matrix<vector<VectorMap<double, vector<Label>>>>(n+1, n, vector<VectorMap<double, vector<Label>>>(NG->L.size()));
	penalties_sum = sum(penalties);
}

void BoundingStructure::AddBound(int k, int r, const Label& l)
{
	insert_sorted(S[k][l.v][r].Insert(floor(l.Thelp()), {}), l,
		[] (const Label& l1, const Label& l2) { return l1.Tdurnum() < l2.Tdurnum(); });
}

double BoundingStructure::CompletionBound(int k, int r, const Label& l) const
{
//	if (!enabled) return l.Tdur;
	int n = vrp->D.VertexCount();
	double T = vrp->T;

	double LB = UB;
	Vertex w = l.v;
	r = (int)NG->L.size()-r-1;
	r -= NG->L[r] != w;
	double Thelpl = l.Thelp();
	double Tdurl = l.Tdurnum();
	for (auto& Thelp_entry: S[n-k+1][w][r])
	{
		double Thelpm = Thelp_entry.first;
		if (epsilon_bigger(T + Thelpm + Thelpl + penalties[w] + penalties_sum, LB)) break;
		for (auto& m: Thelp_entry.second)
		{
			if (epsilon_bigger(m.Tdurnum() + Tdurl + penalties[w] + penalties_sum, LB)) break;
			if (epsilon_smaller(T-m.Ttime, l.Ttime)) continue;
			if (intersection(m.S, l.S) != create_bitset<MAX_N>({w})) continue;
			// Uncomment for TI-Bounding
			if (epsilon_smaller(max(dom(l.Tdur)), T-max(dom(m.Tdur))))
			{
				LB = min(LB, T + m.Thelp() + Thelpl + penalties[w] + penalties_sum);
			}
			else
			{
				PWLFunction lm_duration = l.Tdur + m.Tdur.Compose(T - PWLFunction::IdentityFunction({0.0, T}));
				double dur = lm_duration.Empty() ? INFTY : min(img(lm_duration)) - m.lambda - l.lambda;
				LB = min(LB, dur + penalties[w] + penalties_sum);
			}
//			LB = min(LB, max(m.Tdurnum() + Tdurl + penalties[w] + penalties_sum, T + Thelpm + Thelpl + penalties[w] + penalties_sum));
			continue;

			// TD-Bounding.
			// Merge l and m duration functions lm_d(t) = l_d(t) + m_d(T-t).
			PWLFunction lm_duration = l.Tdur + m.Tdur.Compose(T - PWLFunction::IdentityFunction({0.0, T}));
			double lb = max(T+m.Thelp()+Thelpl, lm_duration.Empty() ? -INFTY : min(img(lm_duration)) - l.lambda - m.lambda) + penalties[w] + penalties_sum;
			LB = min(LB, lb);
		}
	}
	
	return LB;
}

// Expands all non dominated NG routes to get bounds, we do not consider labels which completion bound is
// greater than UB.
// @param vrp: VRP Instance
// @param NG: structure with NG information.
// @param lambda: penalties for vertices.
// @param B: bounding structure.
// @param [out] best_route: route with best cost.
// @param [out] best_cost: cost of best_route.
// @param [out] log: output log to save the execution information.
BoundingStructure run_labeling(const VRPInstance& vrp, const NGStructure& NG, const vector<double>& lambda, Route* best_route, double* best_cost, MLBExecutionLog* log)
{
	BoundingStructure BS((VRPInstance*)&vrp, (NGStructure*)&NG, lambda, 0.0);
	
	int n = vrp.D.VertexCount();
	int R = NG.L.size();
	*best_route = Route({}, 0.0, INFTY);
	*best_cost = INFTY;
	
	stretch_to_size(*log->count_by_length, n+1, 0);
	Stopwatch rolex(true), rolex_domination(false), rolex_extension(false), rolex_queuing(false);
	
	// Initialize queue.
	Matrix<vector<vector<Label>>> q(n, R, vector<vector<Label>>(n));
	q[1][0][vrp.o].push_back(Label(nullptr, vrp.o, create_bitset<MAX_N>({vrp.o}), PWLFunction::ConstantFunction(0.0, vrp.tw[vrp.o]), vrp.a[vrp.o], lambda[vrp.o]));
	for (int k = 1; k < n; ++k)
	{
		for (int r = 0; r < R-1; ++r)
		{
			for (int v = 0; v < n; ++v)
			{
				rolex_queuing.Resume();
				sort(q[k][r][v].begin(), q[k][r][v].end(), [] (Label& l1, Label& l2) { return make_tuple(min(img(l1.Tdur))-l1.lambda, min(img(l1.Tdur))) < make_tuple(min(img(l2.Tdur))-l2.lambda, min(img(l2.Tdur))); });
				rolex_queuing.Pause();
				
				// D[i][S'] = minimum duration of a label M with: k(M)=k, r(M)=r, i(M)=i, S(M)=S'.
				vector<Label*> D; // D[S] are the labels sorted by min(img(Tdur))-lambda.
				for (Label& l: q[k][r][v])
				{
					bool is_dominated = false;
					rolex_domination.Resume();
					PWLDominationFunction Tdurl = l.Tdur;
					for (auto& m: D)
					{
						if (epsilon_bigger(min(img(m->Tdur))-m->lambda, max(img(l.Tdur))-l.lambda)) break;
						if (!is_subset(m->S, l.S)) continue;
						if (Tdurl.DominatePieces(m->Tdur, l.lambda-m->lambda))
						{
							is_dominated = true;
							break;
						}
					}
					rolex_domination.Pause();
					if (is_dominated)
					{
						log->dominated_count++;
						continue;
					}
					l.Tdur = (PWLFunction)Tdurl;
					insert_sorted(D, &l, [] (Label* l1, Label* l2) { return min(img(l1->Tdur))-l1->lambda < min(img(l2->Tdur))-l2->lambda; });
					
					// Save it to B.
//					clog << k << " " << r << " " << l.Path() << " " << l.ToTI() << endl;
					BS.AddBound(k, r, l);
					
					// Extension.
					log->processed_count++;
					rolex_extension.Resume();
					(*log->count_by_length)[k]++;
					Label* ll = new Label(l);
					for (Vertex w: vrp.D.Successors(v))
					{
						// Feasibility check.
						if (l.S.test(w)) continue;
						if (k < vrp.prec_count[w]) continue;
						if (n-k-1 < vrp.suc_count[w]) continue;
						if (n-k-1 < vrp.suc_count[NG.L[r+1]]) continue;
						if (!NG.V[r].test(w)) continue;
						if (epsilon_smaller(vrp.LDT[v][w], min(dom(l.Tdur)))) continue;
						PWLFunction Tdurw = epsilon_smaller(max(dom(l.Tdur)), min(img(vrp.dep[v][w])))
											? PWLFunction::ConstantFunction(l.Tdur(max(dom(l.Tdur))) + vrp.a[w] - max(dom(l.Tdur)), {vrp.a[w], vrp.a[w]})
											: (l.Tdur + vrp.tau[v][w]).Compose(vrp.dep[v][w]);
						if (Tdurw.Empty()) continue;
						VertexSet Sw = intersection(l.S, NG.N[w]);
						Sw.set(w);
						Label lw(ll, w, Sw, Tdurw, min(dom(l.Tdur)), l.lambda + lambda[w]);
						
						// Process tour.
						if (w == vrp.d)
						{
							(*log->count_by_length)[k+1]++;
							double lw_cost = min(img(lw.Tdur))-lw.lambda;
							if (lw_cost < *best_cost)
							{
								*best_cost = lw_cost;
								*best_route = Route(lw.Path(), lw.Tdur.PreValue(min(img(lw.Tdur))), min(img(lw.Tdur)));
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
	return BS;
}

TILabel::TILabel(TILabel* prev, goc::Vertex v, VertexSet S, double Tdur, TimeUnit Thelp, double lambda)
	:prev(prev), v(v), S(S), Tdur(Tdur), Thelp(Thelp), lambda(lambda)
{ }

GraphPath TILabel::Path() const
{
	if (!prev) return {v};
	auto p = prev->Path();
	p.push_back(v);
	return p;
}

void TILabel::Print(ostream& os) const
{
	os << "{P:" << Path() << ", S: " << S << ", Tdur: " << Tdur << ", Thelp: " << Thelp << ", lambda: " << lambda << "}";
}

void run_ti_labeling(const VRPInstance& vrp, const NGStructure& NG, const vector<double>& lambda, Route* best_route, double* best_cost, MLBExecutionLog* log)
{
	int n = vrp.D.VertexCount();
	int R = NG.L.size();
	*best_route = Route({}, 0.0, INFTY);
	*best_cost = INFTY;
	
	stretch_to_size(*log->count_by_length, n+1, 0);
	Stopwatch rolex(true), rolex_domination(false), rolex_extension(false), rolex_queuing(false);
	
	// Initialize queue.
	Matrix<vector<vector<TILabel>>> q(n, R, vector<vector<TILabel>>(n));
	q[1][0][vrp.o].push_back(TILabel(nullptr, vrp.o, create_bitset<MAX_N>({vrp.o}), 0.0-lambda[vrp.o], -vrp.b[vrp.o]-lambda[vrp.o], lambda[vrp.o]));
	for (int k = 1; k < n; ++k)
	{
		for (int r = 0; r < R-1; ++r)
		{
			for (int v = 0; v < n; ++v)
			{
				rolex_queuing.Resume();
				sort(q[k][r][v].begin(), q[k][r][v].end(), [] (TILabel& l1, TILabel& l2) { return make_tuple(l1.Tdur, l1.Thelp) < make_tuple(l2.Tdur, l2.Thelp); });
				rolex_queuing.Pause();
				
				// D[i][S'] = minimum duration of a label M with: k(M)=k, r(M)=r, i(M)=i, S(M)=S'.
				vector<TILabel*> D; // D[S] are the labels sorted by Thelp.
				for (TILabel& l: q[k][r][v])
				{
					bool is_dominated = false;
					rolex_domination.Resume();
					for (auto& m: D)
					{
						if (epsilon_bigger(m->Thelp, l.Thelp)) break;
						if (!is_subset(m->S, l.S)) continue;
						is_dominated = true;
						break;
					}
					rolex_domination.Pause();
					if (is_dominated)
					{
						log->dominated_count++;
						continue;
					}
					insert_sorted(D, &l, [] (TILabel* l1, TILabel* l2) { return l1->Thelp < l2->Thelp; });
					
					// Extension.
					log->processed_count++;
					rolex_extension.Resume();
					(*log->count_by_length)[k]++;
					
					for (Vertex w: vrp.D.Successors(v))
					{
						// Feasibility check.
						if (l.S.test(w)) continue;
						if (k < vrp.prec_count[w]) continue;
						if (n-k-1 < vrp.suc_count[w]) continue;
						if (n-k-1 < vrp.suc_count[NG.L[r+1]]) continue;
						if (!NG.V[r].test(w)) continue;
						double d_vw = vrp.MinimumTravelTime({v,w});
						VertexSet Sw = intersection(l.S, NG.N[w]);
						Sw.set(w);
						TILabel lw(&l, w, Sw,
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
}

BoundingStructure run_dssr(const VRPInstance& vrp, const NGStructure& NG, int max_iter, const vector<double>& lambda,
	Route& UB, double& LB, CGExecutionLog* dssr_log, MLBExecutionLog* exact_log)
{
	Stopwatch rolex(true);
	int n = vrp.D.VertexCount();
	double Lambda = sum(lambda);
	
	// Init Log.
	dssr_log->iteration_count = 0;
	dssr_log->iterations = vector<nlohmann::json>();
	
	// Reverse NG structure and instance.
	((NGStructure*)&NG)->delta = 14;
	VRPInstance back_vrp = reverse_instance(vrp);
	NGStructure back_NG(back_vrp, NG.N, reverse(NG.L), 14);
	
	VRPInstance* vrps[2] = { (VRPInstance*)&vrp, (VRPInstance*)&back_vrp };
	NGStructure* ngs[2] = { (NGStructure*)&NG, &back_NG };
	
	bool changed = true;
	int d = 0;
	while (changed && dssr_log->iteration_count < max_iter)
	{
		changed = false;
		d = 1 - d;

		// Solve labeling.
		double best_cost;
		Route best_route;
		MLBExecutionLog iteration_log(true);
		run_ti_labeling(*vrps[d], *ngs[d], lambda, &best_route, &best_cost, &iteration_log);
		dssr_log->iteration_count++;
		dssr_log->iterations->push_back(iteration_log);

		// Update LB.
		LB = max(LB, best_cost+Lambda);

		// Find cycles.
		auto& P = best_route.path;
		vector<int> last(n, -1);
		vector<pair<int, int>> C(n, {-1, -1});
		bool found_cycles = false;
		for (int i = 1; i < (int)P.size()-1; ++i)
		{
			// If we have already visited vertex P[i], then we have a cycle.
			if (last[P[i]] > -1)
			{
				C[P[i]] = {last[P[i]], i};
				found_cycles = true;
			}
			// Set the last time we visited P[i].
			last[P[i]] = i;
		}

		// If no cycles were found, we have the optimum.
		if (!found_cycles)
		{
			break;
		}

		// Remove cycles from NG.
		for (int i = 1; i < n-1; ++i)
		{
			if (C[i].first == -1 || NG.N[P[i]].count() >= NG.delta) continue;
			for (int j = C[i].first+1; j < C[i].second; ++j)
			{
				ngs[0]->N[P[j]].set(i);
				ngs[1]->N[P[j]].set(i);
				changed = true;
			}
		}

		clog << "DSSR Iteration " << dssr_log->iteration_count << ", LB: " << LB << endl;
	}
	clog << "Finished DSSR with LB: " << LB << " and UB: " << UB.duration << endl;
	Route best_route;
	double best_cost;
	MLBExecutionLog bound_log(true);
	BoundingStructure B(&back_vrp, &back_NG, lambda, UB.duration);
	B = run_labeling(back_vrp, back_NG, lambda, &best_route, &best_cost,&bound_log);
	dssr_log->iteration_count++;
	dssr_log->iterations->push_back(bound_log);
	LB = best_cost + Lambda;
	B.UB = UB.duration;
	clog << "Improvement with NGL finished LB: " << LB << " and UB: " << UB.duration << endl;
	dssr_log->status = CGStatus::Optimum;
	dssr_log->incumbent_value = LB;
	dssr_log->time = rolex.Peek();
	
	UB = run_exact(vrp, back_NG, B, lambda, UB, LB, exact_log);
	clog << "Exact finished in " << exact_log->time << "s, with UB: " << UB.duration << endl;
	return B;
}
} // namespace tdtsptw