//
// Created by Gonzalo Lera Romero on 12/05/2020.
//

#include "preprocess_time_precedence.h"

using namespace std;
using namespace goc;
using namespace nlohmann;

namespace tdtsptw
{
void preprocess_time_precedence(json& instance)
{
	int n = instance["digraph"]["vertex_count"];
	Vertex o = instance["start_depot"];
	Vertex d = instance["end_depot"];
	Matrix<double> LDT = instance["LDT"];
	Matrix<double> EAT = instance["EAT"];
	vector<double> a(n), b(n);
	for (int i = 0; i < n; ++i)
	{
		a[i] = instance["time_windows"][i][0];
		b[i] = instance["time_windows"][i][1];
	}

	// Calculate time_prec and time_succ.
	Matrix<pair<double, int>> time_prec(n, n+1); // time_prec[i][k] = max { t : |{j : LDT(i, j) <= t}| < k }.
	Matrix<pair<double, int>> time_succ(n, n+1); // time_succ[i][k] = min { t : |{j : EAT(j, i) >= t}| < k }.

	// Compute LDT times for the predecesors.
	for (Vertex i = 0; i < n; ++i)
	{
		vector<pair<double, Vertex> > limit_times;
		for (Vertex k = 0; k < n; ++k)
		{
			if (k == i) continue;

			limit_times.push_back(pair<double, Vertex>(LDT[i][k], k));
		}
		sort(limit_times.begin(), limit_times.end());
		for (int k = 1; k < n; k++)
		{
			time_prec[i][k] = pair<double, Vertex>(limit_times[k-1].first, limit_times[k-1].second);
		}
		time_prec[i][n] = {i == d ? b[i] : -INFTY, i};
	}

	// Compute EAT times for the successors
	for (Vertex i = 0; i < n; ++i)
	{
		vector<pair<double, Vertex> > limit_times;
		for (Vertex k = 0; k < n; ++k)
		{
			if (k == i) continue;
			limit_times.emplace_back(pair<double, Vertex>(EAT[k][i], k));
		}
		sort(limit_times.begin(), limit_times.end(), std::greater<>());
		for (int k = 1; k < n; k++)
		{
			time_succ[i][n-k+1] = pair<double, Vertex>(limit_times[k-1].first, limit_times[k-1].second);
		}
		time_succ[i][1] = {i == o ? a[i] : INFTY, i};
	}

	// Time windows with precedence. TWP[k][i] = feasible arrival times to i if k-th position in path.
	Matrix<Interval> TWP(n+1, n);
	for (int k = 1; k <= n; ++k)
		for (Vertex v = 0; v < n; ++v)
			TWP[k][v] = Interval(time_succ[v][k].first, time_prec[v][k].first);
//			TWP[k][v] = {a[v],b[v]};

	instance["time_windows_with_precedences"] = TWP;
}
}