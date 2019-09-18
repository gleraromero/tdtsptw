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
	int count = 0;
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
} // namespace tdtsptw