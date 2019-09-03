//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "labeling.h"

using namespace std;
using namespace goc;

namespace tdtsptw
{
Route run_labeling(const VRPInstance& vrp, const Duration& time_limit, MLBExecutionLog* log)
{
	int n = vrp.D.VertexCount();
	
	*log = MLBExecutionLog(true);
	Stopwatch rolex(true);
	auto comp = [&] (Label* l, Label* m) { return min(img(l->D)) < min(img(m->D)); };
	priority_queue<Label*, vector<Label*>, decltype(comp)> q(comp);
	q.push(new Label{nullptr, vrp.o, 1, create_bitset<MAX_N>({vrp.o}), PWLFunction::ConstantFunction(0.0, vrp.tw[vrp.o])});
	
	Matrix<unordered_map<VertexSet, vector<Label*>>> D(n, n+1);
	
	Route best({}, 0.0, INFTY);
	while (!q.empty())
	{
		Label* l = q.top();
		q.pop();
		
		// Domination step.
		for ()
	}
}
} // namespace tdtsptw