//
// Created by Gonzalo Lera Romero on 07/05/2020.
//

#include "dynamic_neighbour_augmentation.h"
#include "relaxation_solver.h"
#include "label_sequence_td.h"

using namespace std;
using namespace goc;
using namespace nlohmann;

namespace tdtsptw
{
void dynamic_neighbour_augmentation(const VRPInstance& vrp_f, const VRPInstance& vrp_b, NGLInfo& ngl_info_f,
									NGLInfo& ngl_info_b, int delta, const vector<double>& penalties,
									const Duration& time_limit, Route* UB, double* lb, json* log)
{
	// Initialize log.
	Stopwatch rolex(true), rolex_temp(false);
	CGExecutionLog log_dna;
	log_dna.iteration_count = 0;
	log_dna.iterations = vector<json>();
	log_dna.status = CGStatus::Optimum;

	bool broke_cycles = true;
	while (broke_cycles)
	{
		broke_cycles = false;

		// Solve relaxation.
		Route opt;
		double opt_cost;
		json relaxation_log;
		rolex_temp.Reset().Resume();
		BLBStatus status = run_relaxation<LabelSequenceTD>(vrp_f, vrp_b, ngl_info_f, ngl_info_b, penalties, nullptr,
				time_limit - rolex.Peek(), &opt, &opt_cost, &relaxation_log);
		rolex_temp.Pause();
		log_dna.iteration_count++;
		log_dna.iterations->push_back(relaxation_log);
		if (status == BLBStatus::TimeLimitReached) { log_dna.status = CGStatus::TimeLimitReached; break; }
		if (status != BLBStatus::Finished) fail("Unexpected status from relaxation solver.");
		*lb = opt_cost + sum(penalties);
		clog << "> Iteration: " << log_dna.iteration_count << " - LB: " << *lb << " - Time: " << rolex_temp.Peek() << " - Total time: " << rolex.Peek() << endl;
		clog << ">\tRoute: " << opt.path << endl;
		// Check if lower bound reached UB.
		if (epsilon_equal(*lb, UB->duration))
		{
			clog << "> Lower bound reached upper bound." << endl;
			break;
		}
		// If solution is elementary finish algorithm, the optimum has been found.
		if (opt.IsElementary())
		{
			clog << "> Found elementary path while solving relaxation." << endl;
			*UB = opt;
			break;
		}

		// Find vertex-disjoint cycles to break.
		int n = vrp_f.D.VertexCount();
		auto& p = opt.path;
		auto& N = ngl_info_f.N;
		vector<bool> vertex_in_cycle(n, false); // indicates if a vertex was already in a broken cycle.
		for (int i = 0; i < n; ++i)
		{
			if (vertex_in_cycle[p[i]]) continue; // If i is already in a cycle, skip.
			for (int j = i+1; j < n; ++j)
			{
				if (vertex_in_cycle[p[j]]) break; // If j is in a cycle, then the cycle starting at i would not be disjoint.
				if (p[i] != p[j]) continue; // Not a cycle from i to j, skip.
				if (p[i] == p[j]) // Found a cycle disjoint with the previous one, modify neighbours.
				{
					broke_cycles = true;
					clog << ">\tBroke cycle: ";
					for (int k = i; k < j; ++k)
					{
						clog << p[k] << " ";
						vertex_in_cycle[p[k]] = true;
						ngl_info_f.N[p[k]].set(p[i]);
						ngl_info_b.N[p[k]].set(p[i]);
					}
					clog << p[i] << endl;
					break;
				}
				if (!N[p[j]].test(p[i]) && N[p[j]].count() >= delta) break; // vertex j already has a neighbourhood of size delta, then we can not break this cycle.
			}
		}
	}
}
} // namespace tdtsptw