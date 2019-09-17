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
	Matrix<vector<vector<Label>>> q(UB.duration+1, n, vector<vector<Label>>(n, vector<Label>()));
	q[0][1][vrp.o].push_back(Label(nullptr, vrp.o, create_bitset<MAX_N>({vrp.o}), PWLFunction::ConstantFunction(0.0, vrp.tw[vrp.o]), vrp.a[vrp.o], lambda[vrp.o]));
	B.UB = best.duration;
	// Dominance structure: D[v][S] contains labels with v(l)=v, S(l)=S sorted by cost.
	vector<unordered_map<VertexSet, vector<Label>>> D(n);
	vector<Label*> pool;
	
	for (int LB = 0; LB < best.duration; ++LB)
	{
		for (int k = 1; k < n; ++k)
		{
			for (Vertex v: vrp.D.Vertices())
			{
				// Check if LB has reached best.
				if (epsilon_bigger_equal(LB, best.duration)) break;
				
				rolex_queuing.Resume();
				sort(q[LB][k][v].begin(), q[LB][k][v].end(), [](Label& l1, Label& l2) { return make_tuple(min(img(l1.Tdur)), l1.Ttime) < make_tuple(min(img(l2.Tdur)), l2.Ttime); });
				rolex_queuing.Pause();
				
				for (Label& l: q[LB][k][v])
				{
					// Check if LB has reached best.
					if (epsilon_bigger_equal(LB, best.duration)) break;
					
					// Domination.
					rolex_domination.Resume();
					bool is_dominated = false;
					PWLDominationFunction l_D = l.Tdur;
					for (Label& m: D[v][l.S])
					{
						if (epsilon_bigger(min(img(m.Tdur)), max(img(l.Tdur)))) break;
						if (!l_D.DominatePieces(m.Tdur)) continue;
						is_dominated = true;
						break;
					}
					rolex_domination.Pause();
					if (is_dominated)
					{
						log->dominated_count++;
						continue;
					}
					// Get non dominated pieces.
					l.Tdur = (PWLFunction)l_D;
					l.Ttime = min(dom(l.Tdur));
					insert_sorted(D[v][l.S], l, [] (const Label& l1, const Label& l2) { return min(img(l1.Tdur)) < min(img(l2.Tdur));});
					
					// Extension.
					(*log->count_by_length)[k]++;
					log->processed_count++;
					rolex_extension.Resume();
					pool.push_back(new Label(l));
					for (Vertex w: vrp.D.Successors(v))
					{
						// Feasibility check.
						if (l.S.test(w)) continue;
						VertexSet Sw = l.S;
						Sw.set(w);
						double Ttimew = vrp.ArrivalTime({v,w}, l.Ttime);
						bool any_unreachable = false;
						for (Vertex u: vrp.D.Vertices()) if (!Sw.test(u) && epsilon_smaller(vrp.LDT[w][u], Ttimew)) { any_unreachable = true; break; }
						if (any_unreachable) continue;
						
						// Extension.
						PWLFunction Tdurw = epsilon_smaller(max(dom(l.Tdur)), min(img(vrp.dep[v][w])))
										  ? PWLFunction::ConstantFunction(l.Tdur(max(dom(l.Tdur))) + vrp.a[w] - max(dom(l.Tdur)), {vrp.a[w], vrp.a[w]})
										  : (l.Tdur + vrp.tau[v][w]).Compose(vrp.dep[v][w]);
						if (Tdurw.Empty()) continue;
						Label lw(pool.back(), w, Sw, Tdurw, Ttimew, l.lambda + lambda[w]);
						
						// Process tour.
						if (w == vrp.d)
						{
							if (min(img(lw.Tdur)) < best.duration) best = Route(lw.Path(), max(dom(lw.Tdur))-lw.Tdur(max(dom(lw.Tdur))), min(img(lw.Tdur)));
							break;
						}
						
						// Bounding.
						auto LL = reverse(NG.L);
						int r = 0;
						for (r = 0; r < LL.size(); ++r)
							if (!lw.S.test(LL[r]))
							{
								r--;
								break;
							}
//						int r = 0;
//						for (r = 0; r < (int)NG.L.size()-1 && (!lw.S.test(NG.L[r+1]) || NG.L[r+1] == w); ++r) {}
						double LBw = B.CompletionBound(k+1, r, lw);
//						if (LBw < lb) { clog << LBw << " " << lb << endl; fail("Nope");}
//						clog << LBw << endl;
						
//						double LBw = LB;
						if (epsilon_bigger_equal(LBw, best.duration))
						{
							log->enumerated_count++;
							continue;
						}
//						if (epsilon_smaller(LBw, LB)) { clog << LBw << " " << LB << endl; fail("Nope"); }
//						clog << LBw << endl;
						LBw = max(LB, (int)floor(LBw));
						q[LBw][k+1][w].push_back(lw);
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
	log->time = rolex.Peek();
	
	return best;
}
} // namespace tdtsptw