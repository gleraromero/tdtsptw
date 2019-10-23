//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "lbl_exact.h"

using namespace std;
using namespace goc;

namespace tdtsptw
{
State::Piece::Piece(double lb, const LinearFunction& f, Piece* prev, Vertex v)
	: lb(lb), f(f), prev(prev), v(v)
{ }

bool State::Piece::Dominate(Piece& p2)
{
	auto& p = f;
	auto& q = p2.f;
	
	// Remove the case when p is after q.
	if (epsilon_bigger(p.domain.left, q.domain.right)) return false;
	if (epsilon_bigger(p(q.domain.left), q(q.domain.left)) && epsilon_bigger(p(q.domain.right), q(q.domain.right))) return false;

	double t1 = -INFTY, t2 = INFTY;
	if (epsilon_equal(p.slope, q.slope))
	{
		if (epsilon_bigger(p(q.domain.left), q(q.domain.left)))
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
	if (epsilon_bigger(q.domain.right, p.domain.right) && epsilon_bigger(t2, p.domain.right))
	{
		// Intersection between identity function from max(dom(p)) and q.
		t2 = (p(p.domain.right) - p.domain.right - q.intercept) / (q.slope - 1.0);
	}
	
	t1 = max(t1, p.domain.left);
	t2 = min(t2, q.domain.right);

	// q is dominated in [t1, t2].
	// If all q is dominated, return true.
	if (epsilon_smaller_equal(t1, q.domain.left) && epsilon_bigger_equal(t2, q.domain.right))
	{
		return true;
	}
	else if (epsilon_smaller(t2, q.domain.left) || epsilon_bigger(t1, q.domain.right))
	{
		return false;
	}
	// A prefix is dominated.
	else if (epsilon_smaller_equal(t1, q.domain.left))
	{
		q.domain.left = t2+EPS;
		if (epsilon_equal(q.domain.left, q.domain.right))
		{
			q.intercept = q(t2+EPS);
			q.slope = 0.0;
		}
	}
	// A suffix is dominated
	else if (epsilon_bigger_equal(t2, q.domain.right))
	{
		q.domain.right = t1-EPS;
		if (epsilon_equal(q.domain.left, q.domain.right))
		{
			q.intercept = q(t1-EPS);
			q.slope = 0.0;
		}
	}
	return false;
}

GraphPath State::Piece::Path() const
{
	if (!prev) return {v};
	auto P = prev->Path();
	P.push_back(v);
	return P;
}

void State::Piece::Print(ostream& os) const
{
	os << f;
}

bool State::Merge(vector<Piece>& P)
{
	auto& P1 = F;
	auto& P2 = P;
	vector<Piece> F;

	// Sweep pieces.
	int i = 0, j = 0;
	int k = 0;
	bool any_p2 = false;
	while (i < P1.size() && j < P2.size())
	{
		if (k++ > 10000)
		{
			clog.precision(17);
			clog << "Merge: " << k << endl;
			clog << P1[i].f << " " << P2[j].f << endl;
			fail("Merge");
		}
		if (epsilon_smaller_equal(min(dom(P1[i].f)), min(dom(P2[j].f))) && P1[i].Dominate(P2[j])) { ++j; continue; }
		if (epsilon_smaller_equal(min(dom(P2[j].f)), min(dom(P1[i].f))) && P2[j].Dominate(P1[i])) { ++i; continue; }
		
		if (epsilon_smaller(min(dom(P1[i].f)), min(dom(P2[j].f))))
		{
			if (epsilon_smaller_equal(max(dom(P1[i].f)), min(dom(P2[j].f))))
			{
				F.push_back(P1[i]);
				++i;
			}
			else
			{
				F.emplace_back(Piece(P1[i].lb, P1[i].f.RestrictDomain({min(dom(P1[i].f)), min(dom(P2[j].f))}), P1[i].prev, P1[i].v));
				P1[i].f.domain.left = P2[j].f.domain.left+EPS;
			}
		}
		else
		{
			any_p2 = true;
			if (epsilon_smaller_equal(max(dom(P2[j].f)), min(dom(P1[i].f))))
			{
				F.push_back(P2[j]);
				++j;
			}
			else
			{
				F.emplace_back(Piece(P2[j].lb, P2[j].f.RestrictDomain({min(dom(P2[j].f)), min(dom(P1[i].f))}), P2[j].prev, P2[j].v));
				P2[j].f.domain.left = P1[i].f.domain.left+EPS;
			}
		}
	}
	any_p2 |= j < P2.size();
	// Add remaining pieces.
	for (; i < (int)P1.size(); ++i) F.push_back(P1[i]);
	for (; j < (int)P2.size(); ++j) F.push_back(P2[j]);
	
	this->F = F;
	return any_p2;
}
	
void State::DominateBy(const State& s2)
{
	auto P1 = s2.F;
	auto P2 = F;
	F.clear();
	
	// Sweep pieces.
	int i = 0, j = 0;
	int k = 0;
	for (; i < P1.size() && j < P2.size();)
	{
		if (k++ > 10000)
		{
			clog << "DominateBy: " << k << endl;
			fail("DominateBy");
		}
		if (epsilon_smaller_equal(min(dom(P1[i].f)), min(dom(P2[j].f))) && P1[i].Dominate(P2[j])) { ++j; continue; }
		if (epsilon_smaller_equal(min(dom(P2[j].f)), min(dom(P1[i].f))) && P2[j].Dominate(P1[i])) { ++i; continue; }
		
		if (epsilon_smaller(min(dom(P1[i].f)), min(dom(P2[j].f))))
		{
			if (epsilon_smaller_equal(max(dom(P1[i].f)), min(dom(P2[j].f))))
			{
				++i;
			}
			else
			{
				P1[i].f.domain.left = P2[j].f.domain.left+EPS;
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
				F.push_back(Piece(P2[j].lb, P2[j].f.RestrictDomain({min(dom(P2[j].f)), min(dom(P1[i].f))}), P2[j].prev, P2[j].v));
				P2[j].f.domain.left = P1[i].f.domain.left+EPS;
			}
		}
	}
	// Add remaining pieces.
	for (; j < P2.size(); ++j) F.push_back(P2[j]);
}

void State::Print(ostream& os) const
{
	os << F << endl;
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
	B[n-k+1][R-r-(!LSet.test(v))][v].push_back({S, DeltaR});
}

void Bounding::Bound(goc::Vertex v, const VertexSet& S, State& Delta)
{
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
			if (Delta.F[i].lb == -INFTY) { ++i; continue; }

			auto& p_i = Delta.F[i].f;
			auto& p_j = Delta2[j];
			
			if (epsilon_bigger(min(dom(p_i)), max(dom(p_j))))
			{
				// Case 1: p_i can not be extended with p_j, skip p_j.
				++j;
			}
			else if (epsilon_bigger(min(dom(p_j)), max(dom(p_i))))
			{
				// Case 2: p_j comes after p_i then p_i is bounded by waiting.
				double ti = p_i(max(dom(p_i))), tj = p_j(min(dom(p_j)));
				double wait = min(dom(p_j)) - max(dom(p_i));
				Delta.F[i].lb = min(Delta.F[i].lb, ti+wait+tj);
				++i;
			}
			else
			{
				// Case 3: p_i overlaps p_j, we must get the best bound.
				double t0 = max(min(dom(p_i)), min(dom(p_j)));
				double t1 = min(max(dom(p_i)), max(dom(p_j)));
				Delta.F[i].lb = min(Delta.F[i].lb, min(p_i(t0)+p_j(t0), p_i(t1)+p_j(t1)));
				if (epsilon_smaller_equal(max(dom(p_i)), max(dom(p_j)))) ++i;
				else ++j;
			}
		}
	}
}

Route run_dna(const VRPInstance& vrp, const VRPInstance& rvrp, NGStructure& NG, NGStructure& rNG, const vector<double>& lambda, CGExecutionLog* log, double& LB, Duration time_limit, bool bidirectional)
{
	int n = vrp.D.VertexCount();
	int max_iter = 200;
	int delta = 14;
	Stopwatch rolex(true);
	log->iteration_count = 0;
	log->iterations = vector<nlohmann::json>();
	for (int it = 0; it < max_iter; ++it)
	{
		clog << "\tIteration: " << it << "\tLB: " << LB << endl;
		
		// Run NGL.
		Route R;
		if (!bidirectional)
		{
			MLBExecutionLog iteration_log(true);
			R = run_ngltd(vrp, NG, lambda, &iteration_log, nullptr, LB, time_limit - rolex.Peek());
			log->iterations->push_back(iteration_log);
			if (iteration_log.status == MLBStatus::TimeLimitReached) { log->status = CGStatus::TimeLimitReached; break; }
		}
		else
		{
			BLBExecutionLog iteration_log(true);
			R = run_ngltd_bidirectional(vrp, rvrp, NG, rNG, lambda, &iteration_log, LB, time_limit - rolex.Peek());
			log->iterations->push_back(iteration_log);
			if (iteration_log.status == BLBStatus::TimeLimitReached) { log->status = CGStatus::TimeLimitReached; break; }
		}
		log->iteration_count++;

		// Find cycles.
		auto& P = R.path;
		vector<int> last(n, -1);
		vector<pair<int, int>> C;
		bool found_cycles = false;
		int k = 0; // all vertices from k+1 to i have |N(.)| < Delta.
		for (int i = 1; i < (int)P.size() - 1; ++i)
		{
			// If we have already visited vertex P[i], then we have a cycle.
			found_cycles |= last[P[i]] > -1;
			
			// If the last appearance of P[i] was after the last enhanceable
			if (last[P[i]] >= k) C.push_back({last[P[i]], i});
			
			// Set the last time we visited P[i].
			last[P[i]] = i;
			
			if (NG.N[P[i]].count() >= delta) k = i;
		}
		
		// If no cycles were found, we have the optimum.
		if (!found_cycles)
		{
			log->time = rolex.Peek();
			log->incumbent_value = LB;
			return R;
		}
		// If no breakable cycles are found, we must stop.
		if (C.empty())
		{
			clog << "\tReached NG neighbour size limit." << endl;
			break;
		}
		
		// Remove vertex-disjoint cycles from NG.
		VertexSet V;
		for (auto& c: C)
		{
			bool disjoint = true;
			for (int j = c.first + 1; j < c.second; ++j)
			{
				if (V.test(P[j]))
				{
					disjoint = false;
					break;
				}
			}
			if (!disjoint) continue;
			for (int j = c.first + 1; j < c.second; ++j)
			{
				NG.N[P[j]].set(P[c.first]);
				rNG.N[P[j]].set(P[c.first]);
				V.set(P[j]);
			}
		}
	}
	log->time = rolex.Peek();
	log->incumbent_value = LB;
	
	return Route({}, -1, INFTY);
}

struct Merger
{
	vector<double> lambda;
	GraphPath L;
	Matrix<vector<pair<VertexSet, State>>> F;
	
	Merger(int n, const vector<double>& lambda, const GraphPath& L) : lambda(lambda), L(L)
	{
		F = Matrix<vector<pair<VertexSet, State>>>(n, n);
	}
	
	void AddForward(const VertexSet& S, State& s, Vertex v, int r)
	{
		F[r][v].push_back({S, s});
	}
	
	void MergeBackward(const VertexSet& S, State& s, Vertex v, int r, double& min_cost, GraphPath& min_path)
	{
		if (s.F.empty()) return;
		int F_r = L.size() - 1 - r;
		if (L[F_r] != v) F_r--;
		vector<LinearFunction> Delta2;
		for (int i = s.F.size()-1; i >= 0; --i) Delta2.push_back(LinearFunction({-s.F[i].f.domain.right, s.F[i].f(s.F[i].f.domain.right)+lambda[v]}, {-s.F[i].f.domain.left, s.F[i].f(s.F[i].f.domain.left)+lambda[v]}));
		for (auto& S_s: F[F_r][v])
		{
			auto& Sf = S_s.first;
			auto& Delta = S_s.second;
			if ((Sf & S) != create_bitset<MAX_N>({v})) continue;
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
					// Case 2: p_j comes after p_i then p_i is bounded by waiting.
					double ti = p_i(max(dom(p_i))), tj = p_j(min(dom(p_j)));
					double wait = min(dom(p_j)) - max(dom(p_i));
					if (ti+wait+tj < min_cost)
					{
						min_cost = ti + wait + tj;
						min_path = Delta.F[i].Path();
						auto Pb = s.F[s.F.size()-j-1].Path();
						Pb.pop_back();
						Pb = reverse(Pb);
						min_path.insert(min_path.end(), Pb.begin(), Pb.end());
					}
					++i;
				}
				else
				{
					// Case 3: p_i overlaps p_j, we must get the best bound.
					double t0 = max(min(dom(p_i)), min(dom(p_j)));
					double t1 = min(max(dom(p_i)), max(dom(p_j)));
					if (min_cost > min(p_i(t0)+p_j(t0), p_i(t1)+p_j(t1)))
					{
						min_cost = min(p_i(t0)+p_j(t0), p_i(t1)+p_j(t1));
						min_path = Delta.F[i].Path();
						auto Pb = s.F[s.F.size()-j-1].Path();
						Pb.pop_back();
						Pb = reverse(Pb);
						min_path.insert(min_path.end(), Pb.begin(), Pb.end());
					}
					if (epsilon_smaller_equal(max(dom(p_i)), max(dom(p_j)))) ++i;
					else ++j;
				}
			}
		}
	}
};

Route run_ngltd_bidirectional(const VRPInstance& fvrp, const VRPInstance& bvrp, const NGStructure& fNG, const NGStructure& bNG, const vector<double>& lambda, BLBExecutionLog* blb_log, double& LB, Duration time_limit)
{
	Stopwatch rolex_blb(true);
	int n = fvrp.D.VertexCount();
	auto& N = fNG.N;
	vector<vector<vector<spp::sparse_hash_map<VertexSet, State>>>> DS[2] = {
		vector<vector<vector<spp::sparse_hash_map<VertexSet, State>>>>(n, vector<vector<spp::sparse_hash_map<VertexSet, State>>>(fNG.L.size(), vector<spp::sparse_hash_map<VertexSet, State>>(n))),
		vector<vector<vector<spp::sparse_hash_map<VertexSet, State>>>>(n, vector<vector<spp::sparse_hash_map<VertexSet, State>>>(fNG.L.size(), vector<spp::sparse_hash_map<VertexSet, State>>(n)))
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
		
		// Init structures.
		double OPT_end = INFTY, OPT_dur = INFTY;
		
		Stopwatch rolex(true), rolex_domination(false), rolex_extension(false), rolex_bounding(false), rolex_queuing(false);
		stretch_to_size(*log.count_by_length, n+1, 0);
		log.status = MLBStatus::Finished;
		
		State::Piece p0(-1, {{vrp.a[vrp.o], -lambda[vrp.o]}, {vrp.b[vrp.o], -lambda[vrp.o]}}, nullptr, vrp.o);
		vector<State::Piece> P0 = {p0};
		D[1][0][vrp.o].insert({create_bitset<MAX_N>({vrp.o}), State()}).first->second.Merge(P0);
		
		for (int k = 1; k <= nn[d]; ++k)
		{
			for (int r = 0; r < (int)L.size(); ++r)
			{
				for (int v = 0; v < n; ++v)
				{
					for (auto& S_Delta: D[k][r][v])
					{
						if (rolex.Peek() >= time_limit) { log.status = MLBStatus::TimeLimitReached; break; }
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
						
						log.processed_count++;
						log.extended_count += Delta.F.size();
						log.enumerated_count += Delta.F.size();
						(*log.count_by_length)[k] += Delta.F.size();
						
						// Extension of pieces.
						rolex_extension.Resume();
						for (Vertex w: vrp.D.Successors(v))
						{
							if (S.test(w)) continue;
							if (!V[r].test(w)) continue;
							if (vrp.prec_count[w] > k) continue;
							if (vrp.suc_count[w] > n - k - 1) continue;
							if (w != L[r + 1] && vrp.suc_count[L[r + 1]] > n - k - 2) continue;
							
							double LDTw_at_v = vrp.DepartureTime({v,w}, NG.tw[k+1][w].right);
							
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
									if (epsilon_bigger(min(dom(tau_j)), LDTw_at_v)) break;
									if (epsilon_smaller(max(dom(tau_j)), min(dom(p_i)))) continue;
									if (epsilon_bigger(min(dom(tau_j)), max(dom(p_i)))) break;
									
									// Find intersection of pieces.
									double t1 = min(LDTw_at_v, max(min(dom(tau_j)), min(dom(p_i))));
									double t2 = min(LDTw_at_v, min(max(dom(tau_j)), max(dom(p_i))));
									t2 = max(t1, t2); // We must avoid numerical errors where (t1 > t2 but t1 <=eps t2).
									
									// Extend.
									LinearFunction pp({t1 + tau_j(t1), p_i(t1) + tau_j(t1) - lambda[w]},
													  {t2 + tau_j(t2), p_i(t2) + tau_j(t2) - lambda[w]});
									
									if (k == n - 1) // If complete tour, add it to the solution.
									{
										if (epsilon_smaller(pp(min(dom(pp))), OPT_dur))
										{
											OPT_dur = pp(min(dom(pp)));
											OPT_end = min(dom(pp));
										}
										if (epsilon_smaller(pp(max(dom(pp))), OPT_dur))
										{
											OPT_dur = pp(max(dom(pp)));
											OPT_end = max(dom(pp));
										}
									}
									else // Otherwise, add it to the queue.
									{
										if (!EXT_P.empty() && epsilon_equal(min(dom(EXT_P.back().f)), min(dom(pp))))
											EXT_P.pop_back(); // Remove waiting time piece.
										EXT_P.push_back(State::Piece(-1, pp, &p, w));
									}
									
									// Stop if tau_vw exceeds the latest departure time.
									if (epsilon_bigger_equal(max(dom(tau_j)), max(dom(p_i)))) break;
								}
							}
							if (!EXT_P.empty())
							{
								rolex_extension.Pause();
								rolex_queuing.Resume();
								auto S_w = S & N[w];
								S_w.set(w);
								D[k + 1][r + (w == L[r + 1])][w].insert({S_w, State()}).first->second.Merge(EXT_P);
								rolex_queuing.Pause();
								rolex_extension.Resume();
							}
						}
						rolex_extension.Pause();
					}
				}
			}
		}
		
		log.extension_time = rolex_extension.Peek();
		log.bounding_time = rolex_bounding.Peek();
		log.domination_time = rolex_domination.Peek();
		log.queuing_time = rolex_queuing.Peek();
		log.time = rolex.Peek();
		if (d == 0) blb_log->forward_log = log;
		else blb_log->backward_log = log;
		if (rolex.Peek() >= time_limit) { blb_log->status = BLBStatus::TimeLimitReached; break; }
	}

	if (blb_log->status != BLBStatus::TimeLimitReached)
	{
		// Merge.
		Stopwatch rolex_merge(true);
		Merger M(n, lambda, fNG.L);
		int R = fNG.L.size();
		double min_cost = INFTY;
		for (int r = 0; r < R - 1; ++r)
			for (Vertex v: fvrp.D.Vertices())
				for (auto &S: DS[0][nn[0]][r][v])
					M.AddForward(S.first, S.second, v, r);

		GraphPath min_path;
		for (int r = 0; r < R - 1; ++r)
			for (Vertex v: fvrp.D.Vertices())
				for (auto &S: DS[1][nn[1]][r][v])
					M.MergeBackward(S.first, S.second, v, r, min_cost, min_path);

		VertexSet NGSet;
		for (auto &v: min_path)
		{
			NGSet = (NGSet & fNG.N[v]);
			NGSet.set(v);
		}
		blb_log->merge_time = rolex_merge.Peek();
		blb_log->time = rolex_blb.Peek();
		LB = max(LB, min_cost + sum(lambda));
		return fvrp.BestDurationRoute(min_path);
	}
	return Route({}, 0.0, INFTY);
}

Route run_ngltd(const VRPInstance& vrp, const NGStructure& NG, const vector<double>& lambda, MLBExecutionLog* log, Bounding* B, double& LB, Duration time_limit)
{
	int n = vrp.D.VertexCount();
	auto& N = NG.N;
	auto& L = NG.L;
	auto& V = NG.V;
	
	// Init structures.
	auto D = vector<vector<vector<spp::sparse_hash_map<VertexSet, State>>>>(n, vector<vector<spp::sparse_hash_map<VertexSet, State>>>(L.size(), vector<spp::sparse_hash_map<VertexSet, State>>(n)));
	double opt_cost = INFTY;
	GraphPath opt_path;
	
	Stopwatch rolex(true), rolex_domination(false), rolex_extension(false), rolex_bounding(false), rolex_queuing(false);
	stretch_to_size(*log->count_by_length, n+1, 0);
	log->status = MLBStatus::Finished;
	
	State::Piece p0(-1, {{vrp.a[vrp.o], -lambda[vrp.o]}, {vrp.b[vrp.o], -lambda[vrp.o]}}, nullptr, vrp.o);
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
					if (rolex.Peek() >= time_limit) { log->status = MLBStatus::TimeLimitReached; break; }
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
					if (B) B->AddBound(k, r, v, S, Delta);
					
					log->processed_count++;
					log->extended_count += Delta.F.size();
					log->enumerated_count += Delta.F.size();
					(*log->count_by_length)[k] += Delta.F.size();
					
					// Extension of pieces.
					rolex_extension.Resume();
					for (Vertex w: vrp.D.Successors(v))
					{
						if (S.test(w)) continue;
						if (!V[r].test(w)) continue;
						if (vrp.prec_count[w] > k) continue;
						if (vrp.suc_count[w] > n - k - 1) continue;
						if (w != L[r + 1] && vrp.suc_count[L[r + 1]] > n - k - 2) continue;
						
						double LDTw_at_v = vrp.DepartureTime({v,w}, NG.tw[k+1][w].right);
						
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
								if (epsilon_bigger(min(dom(tau_j)), LDTw_at_v)) break;
								if (epsilon_smaller(max(dom(tau_j)), min(dom(p_i)))) continue;
								if (epsilon_bigger(min(dom(tau_j)), max(dom(p_i)))) break;
								
								// Find intersection of pieces.
								double t1 = min(LDTw_at_v, max(min(dom(tau_j)), min(dom(p_i))));
								double t2 = min(LDTw_at_v, min(max(dom(tau_j)), max(dom(p_i))));
								t2 = max(t1, t2); // We must avoid numerical errors where (t1 > t2 but t1 <=eps t2).
								
								// Extend.
								LinearFunction pp({t1 + tau_j(t1), p_i(t1) + tau_j(t1) - lambda[w]},
												  {t2 + tau_j(t2), p_i(t2) + tau_j(t2) - lambda[w]});
								
								if (k == n - 1) // If complete tour, add it to the solution.
								{
									if (epsilon_smaller(pp(min(dom(pp))), opt_cost))
									{
										opt_cost = pp(min(dom(pp)));
										opt_path = p.Path();
										opt_path.push_back(vrp.d);
									}
									if (epsilon_smaller(pp(max(dom(pp))), opt_cost))
									{
										opt_cost = pp(max(dom(pp)));
										opt_path = p.Path();
										opt_path.push_back(vrp.d);
									}
								}
								else // Otherwise, add it to the queue.
								{
									if (!EXT_P.empty() && epsilon_equal(min(dom(EXT_P.back().f)), min(dom(pp))))
										EXT_P.pop_back(); // Remove waiting time piece.
									EXT_P.push_back(State::Piece(-1, pp, &p, w));
								}
								
								// Stop if tau_vw exceeds the latest departure time.
								if (epsilon_bigger_equal(max(dom(tau_j)), max(dom(p_i)))) break;
							}
						}
						if (!EXT_P.empty())
						{
							rolex_extension.Pause();
							rolex_queuing.Resume();
							auto S_w = S & N[w];
							S_w.set(w);
							D[k + 1][r + (w == L[r + 1])][w].insert({S_w, State()}).first->second.Merge(EXT_P);
							rolex_queuing.Pause();
							rolex_extension.Resume();
						}
					}
					rolex_extension.Pause();
				}
			}
		}
	}

	log->extension_time = rolex_extension.Peek();
	log->bounding_time = rolex_bounding.Peek();
	log->domination_time = rolex_domination.Peek();
	log->queuing_time = rolex_queuing.Peek();
	log->time = rolex.Peek();

	if (log->status == MLBStatus::TimeLimitReached) return Route({}, 0.0, INFTY);
	
	LB = max(LB, opt_cost + sum(lambda));
	return vrp.BestDurationRoute(opt_path);
}

Route run_exact_piecewise(const VRPInstance& vrp, const GraphPath& L, const vector<double>& lambda,
	double LB, double UB, MLBExecutionLog* log, Bounding* B, Duration time_limit)
{
	int n = vrp.D.VertexCount();
	int BASE = floor(LB), TOP = ceil(UB);
	
	// Init structures (D[k][v] -> S -> State)
	auto D = vector<vector<spp::sparse_hash_map<VertexSet, State>>>(n, vector<spp::sparse_hash_map<VertexSet, State>>(n));
	double OPT_end = INFTY, OPT_dur = INFTY;

	Stopwatch rolex(true), rolex_domination(false), rolex_extension(false), rolex_bounding(false), rolex_queuing(false);
	stretch_to_size(*log->count_by_length, n+1, 0);
	log->status = MLBStatus::Finished;
	
	State::Piece p0(LB, {{vrp.a[vrp.o], -lambda[vrp.o]}, {vrp.b[vrp.o], -lambda[vrp.o]}}, nullptr, vrp.o);
	vector<State::Piece> P0 = {p0};
	D[1][vrp.o].insert({create_bitset<MAX_N>({vrp.o}), State()}).first->second.Merge(P0);

	vector<vector<vector<spp::sparse_hash_set<VertexSet>>>> q(TOP-BASE+1, vector<vector<spp::sparse_hash_set<VertexSet>>>(n, vector<spp::sparse_hash_set<VertexSet>>(n)));
	q[0][1][vrp.o].insert(create_bitset<MAX_N>({vrp.o}));
	for (int lb = 0; lb <= TOP-BASE; ++lb)
	{
		for (int k = 1; k < n; ++k)
		{
			for (int v = 0; v < n; ++v)
			{
				for (auto& S: q[lb][k][v])
				{
					if (rolex.Peek() >= time_limit) { log->status = MLBStatus::TimeLimitReached; break; }
					rolex_extension.Resume();
					auto& Delta = D[k][v][S];

					// Get non dominated pieces.
					log->processed_count++;
					
					// Apply bounds.
					vector<LinearFunction> EXT; // Pieces with lb = lb that must be extended.
					rolex_extension.Pause();
					rolex_bounding.Resume();
					if (B) B->Bound(v, S, Delta);
					int next_lb = TOP-BASE+1;
					for (auto& p: Delta.F)
					{
						if (epsilon_bigger_equal(p.lb, UB)) continue;
						if (p.lb == -INFTY) continue;
						if (epsilon_smaller_equal((int)floor(p.lb+EPS), lb + BASE) || !B)
						{
							log->enumerated_count++;
							EXT.push_back(p.f);
							p.lb = -INFTY; // LB = -INFTY means it was already processed.
						}
						else
						{
							next_lb = min(next_lb, (int)floor(p.lb+EPS)-BASE);
						}
//						if (epsilon_smaller(floor(p.lb+EPS), lb + BASE) && p.lb != -INFTY)
//						{
//							clog.precision(17);
//							clog << p.lb+EPS << " " << lb + BASE << endl;
//							fail("Bounds are not increasing.");
//						}
					}
					if (next_lb <= TOP-BASE) q[next_lb][k][v].insert(S);
					rolex_bounding.Pause();
					if (EXT.empty()) continue;
					rolex_extension.Resume();
					log->extended_count += EXT.size();
					(*log->count_by_length)[k] += EXT.size();

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
								
								// Check that tau_j \cap p_i \cap [-INFTY, LDTw_at_v] \neq \emptyset.
								if (epsilon_bigger(min(dom(tau_j)), LDTw_at_v)) break;
								if (epsilon_smaller(max(dom(tau_j)), min(dom(p_i)))) continue;
								if (epsilon_bigger(min(dom(tau_j)), max(dom(p_i)))) break;
								
								// Find intersection of pieces.
								double t1 = min(LDTw_at_v, max(min(dom(tau_j)), min(dom(p_i))));
								double t2 = min(LDTw_at_v, min(max(dom(tau_j)), max(dom(p_i))));
								t2 = max(t1, t2); // We must avoid numerical errors where (t1 > t2 but t1 <=eps t2).
								
								// Extend.
								LinearFunction pp({t1+tau_j(t1), p_i(t1)+tau_j(t1)-lambda[w]}, {t2+tau_j(t2), p_i(t2)+tau_j(t2)-lambda[w]});
								if (k == n-1) // If complete tour, add it to the solution.
								{
									if (epsilon_smaller(pp(min(dom(pp))), OPT_dur))
									{
										OPT_dur = pp(min(dom(pp)));
										OPT_end = min(dom(pp));
									}
									if (epsilon_smaller(pp(max(dom(pp))), OPT_dur))
									{
										OPT_dur = pp(max(dom(pp)));
										OPT_end = max(dom(pp));
									}
								}
								else // Otherwise, add it to the queue.
								{
									if (!EXT_P.empty() && epsilon_equal(min(dom(EXT_P.back().f)) , min(dom(pp)))) EXT_P.pop_back(); // Remove waiting time piece.
									EXT_P.push_back(State::Piece(-1, pp, nullptr, w));
								}
								
								// Stop if tau_vw exceeds the latest departure time.
								if (epsilon_bigger_equal(max(dom(tau_j)), max(dom(p_i)))) break;
							}
						}
						if (!EXT_P.empty())
						{
							rolex_extension.Pause();
							rolex_queuing.Resume();
							*log->positive_domination_time += 1.0_sec;
							if (D[k + 1][w].insert({S_w, State()}).first->second.Merge(EXT_P))
							{
								q[lb][k + 1][w].insert(S_w);
								*log->negative_domination_time += 1.0_sec;
							}
							rolex_queuing.Pause();
							rolex_extension.Resume();
						}
					}
					rolex_extension.Pause();
				}
			}
			q[lb][k].clear();
		}
		q[lb].clear();
		if (OPT_dur < INFTY) break;
	}
	log->bounded_count = log->enumerated_count - log->extended_count;
	log->extension_time = rolex_extension.Peek();
	log->bounding_time = rolex_bounding.Peek();
	log->domination_time = rolex_domination.Peek();
	log->queuing_time = rolex_queuing.Peek();
	log->time = rolex.Peek();
	
	if (log->status == MLBStatus::Finished)
	{
		// Rebuild solution.
		GraphPath P;
		if (OPT_dur < INFTY)
		{
			P.push_back(vrp.d);
			double t = OPT_end, d = OPT_dur;
			int r = L.size() - 2;
			VertexSet E;
			E.set(vrp.d);
			int last_k = n+1;
			for (Vertex w = P.back(); P.size() < n; w = P.back())
			{
				int k = n - E.count();
				if (last_k == k)
				{
					fail("Could not rebuild solution.");
				}
				last_k = k;

				for (Vertex v: vrp.D.Predecessors(w))
				{
					if (E.test(v)) continue;
					double t_v = vrp.DepartureTime({v, w}, t);
					if (t_v == INFTY) continue;
					double d_v = d + lambda[w] - (t - t_v);
					bool found = false;

					for (auto& S_Delta: D[k][v])
					{
						auto& S = S_Delta.first;
						auto& Delta = S_Delta.second;
						if ((S & E) != 0) continue;
						for (auto& p: Delta.F)
						{
							found = (p.f.domain.Includes(t_v) && epsilon_smaller_equal(p.f(t_v)-EPS, d_v)) ||
									(epsilon_smaller(p.f.domain.right, t_v) &&
									 epsilon_smaller_equal(p.f(max(dom(p.f))) + t_v - max(dom(p.f))-EPS, d_v));
							if (found)
							{
								E.set(v);
								P.push_back(v);
								t = t_v;
								d = d_v;
								r = r - (v == L[r]);
								break;
							}
						}
						if (found) break;
					}
					if (found) break;
				}
			}
		}
		
		// Rebuild solution.
		return P.empty() ? Route({}, 0.0, UB) : vrp.BestDurationRoute(reverse(P));
	}
	return Route({}, 0.0, UB);
}
} // namespace tdtsptw