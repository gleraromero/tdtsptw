//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "labeling.h"
#include "pwl_domination_function.h"

using namespace std;
using namespace goc;

namespace tdtsptw
{
Route run_labeling(const VRPInstance& vrp, const Duration& time_limit, MLBExecutionLog* log)
{
	int n = vrp.D.VertexCount();
	
	*log = MLBExecutionLog(true);
	Stopwatch rolex(true), rolex2(false), rolex3(false);
	auto comp = [&] (Label* l, Label* m) { return l->D_min > m->D_min; };
	priority_queue<Label*, vector<Label*>, decltype(comp)> q(comp);
	q.push(new Label{nullptr, vrp.o, 1, create_bitset<MAX_N>({vrp.o}), PWLFunction::ConstantFunction(0.0, {vrp.tw[vrp.o].left, vrp.tw[vrp.o].left}),vrp.tw[vrp.o].left, 0.0});

	Matrix<unordered_map<VertexSet, vector<Label*>>> D(n, n+1);
	for (Vertex v: vrp.D.Vertices()) clog << v << ") " << vrp.tw[v] << endl;
	clog << vrp.LDT << endl;
	Route best({}, 0.0, INFTY);
	while (!q.empty())
	{
		if (rolex.Peek() >= time_limit) { log->status = MLBStatus::TimeLimitReached; break; }

		rolex2.Reset().Resume();
		Label* l = q.top();
		q.pop();
		*log->queuing_time += rolex2.Peek();

		// If l is a tour, then check if it is the best solution.
		if (l->k == n)
		{
			if (l->D_min < best.duration) best = Route(l->Path(), l->D.PreValue(l->D_min)-l->D_min, l->D_min);
			delete l;
			continue;
		}

		// Bounding step.
		if (epsilon_bigger_equal(l->D_min, best.duration)) { delete l; continue; }

		// Domination step.
		rolex2.Reset().Resume();
		PWLDominationFunction l_D = l->D;
		for (Label* m: D[l->v][l->k][l->S])
		{
			if (epsilon_smaller(max(img(l->D)), m->D_min)) break;
			if (l_D.DominatePieces(m->D)) break;
		}
		bool is_dominated = l_D.Empty();
		if (!is_dominated) { l->D = (PWLFunction)l_D; l->t = min(dom(l->D)); l->D_min = min(img(l->D)); }
		*log->domination_time += rolex2.Peek();
		if (is_dominated) { *log->positive_domination_time += rolex2.Peek(); log->dominated_count++; }
		else { *log->negative_domination_time += rolex2.Peek(); log->processed_count++; }
		if (is_dominated) { delete l; continue; } // If D(l) == \emptyset, then l is dominated.

		// l is not dominated, add it to the domination structure.
		rolex2.Reset().Resume();
		insert_sorted(D[l->v][l->k][l->S], l, [] (Label* l1, Label* l2) { return l1->D_min < l2->D_min; });
		*log->process_time += rolex2.Peek();

		// Extend label l.
		rolex2.Reset().Resume();
		for (Vertex w: vrp.D.Successors(l->v))
		{
			// Check feasibility.
			if (l->S.test(w)) continue;
			if (epsilon_bigger(l->t, vrp.LDT[l->v][w])) continue;
			TimeUnit t_w = vrp.ArrivalTime({l->v, w}, l->t);
			for (Vertex u: vrp.D.Vertices()) if (!l->S.test(u) && u != w && epsilon_smaller(vrp.LDT[w][u], t_w)) continue;
			PWLFunction D_w = epsilon_smaller(max(dom(l->D)), min(img(vrp.dep[l->v][w])))
						   ? PWLFunction::ConstantFunction(l->D(max(dom(l->D))) + min(vrp.tw[w]) - max(dom(l->D)), {min(vrp.tw[w]), min(vrp.tw[w])})
						   : (l->D + vrp.tau[l->v][w]).Compose(vrp.dep[l->v][w]);
			if (D_w.Empty()) continue;
			rolex2.Pause();
			rolex3.Reset().Resume();
			q.push(new Label(l, w, l->k+1, unite(l->S, {w}), D_w, min(dom(D_w)), min(img(D_w))));
			*log->queuing_time += rolex3.Peek();
			rolex2.Resume();
			log->extended_count++;
		}
		*log->extension_time += rolex2.Peek();
	}
	log->time = rolex.Peek();
	if (log->status == MLBStatus::DidNotStart) log->status = MLBStatus::Finished;

	return best;
}
} // namespace tdtsptw