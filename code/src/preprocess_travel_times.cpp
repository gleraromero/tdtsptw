//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "preprocess_travel_times.h"

using namespace std;
using namespace goc;
using namespace nlohmann;

namespace tdtsptw
{
namespace
{
// Calculates the time to depart to traverse arc e arriving at tf.
// Returns: INFTY if it is infeasible to depart inside the horizon.
double departing_time(const Matrix<int> C, const vector<Interval>& T, const Matrix<double>& S, const Matrix<double>& D, Arc e, double tf)
{
	double d = D[e.tail][e.head];
	auto& speed = S[C[e.tail][e.head]];
	double t = tf;
	for (int k = (int)T.size()-1; k >= 0; --k)
	{
		if (epsilon_equal(d, 0.0)) break;
		if (epsilon_bigger(T[k].left, tf)) continue;
		double remaining_time_in_k = min(T[k].right, tf) - T[k].left;
		double time_to_complete_d_in_k = d / speed[k];
		double time_in_k = min(remaining_time_in_k, time_to_complete_d_in_k);
		t -= time_in_k;
		d -= time_in_k * speed[k];
	}
	if (epsilon_bigger(d, 0.0)) return INFTY;
	return t;
}


// Calculates the travel time to traverse arc e departing at t0.
// Returns: INFTY if it is infeasible to arrive inside the horizon.
double travel_time(const Matrix<int> C, const vector<Interval>& T, const Matrix<double>& S, const Matrix<double>& D, Arc e, double t0)
{
	double d = D[e.tail][e.head];
	auto& speed = S[C[e.tail][e.head]];
	double t = t0;
	for (int k = 0; k < T.size(); ++k)
	{
		if (epsilon_equal(d, 0.0)) break;
		if (epsilon_smaller(T[k].right, t0)) continue;
		double remaining_time_in_k = T[k].right - max(T[k].left, t0);
		double time_to_complete_d_in_k = d / speed[k];
		double time_in_k = min(remaining_time_in_k, time_to_complete_d_in_k);
		t += time_in_k;
		d -= time_in_k * speed[k];
	}
	if (epsilon_bigger(d, 0.0)) return INFTY;
	return t-t0;
}

// Precondition: no speeds are 0.
PWLFunction compute_travel_time_function(const Matrix<int> C, const vector<Interval>& TZ, const Matrix<double>& S, const Matrix<double>& D, Arc e)
{
	// Calculate speed breakpoints.
	vector<double> speed_breakpoints;
	for (auto& z: TZ) speed_breakpoints.push_back(z.left);
	speed_breakpoints.push_back(TZ.back().right);

	// Travel time breakpoints are two sets
	// 	- B1: speed breakpoints which are feasible to depart
	// 	- B2: times t such that we arrive to head(e) at a speed breakpoint.
	vector<double> B1;
	B1.reserve(speed_breakpoints.size());
	for (double t: speed_breakpoints)
		if (travel_time(C, TZ, S, D, e, t) != INFTY)
			B1.push_back(t);

	vector<double> B2;
	B2.reserve(speed_breakpoints.size());
	for (double t: speed_breakpoints)
		if (departing_time(C, TZ, S, D, e, t) != INFTY)
			B2.push_back(departing_time(C, TZ, S, D, e, t));

	// Merge breakpoints in order in a set B.
	vector<double> B(B1.size()+B2.size());
	merge(B1.begin(), B1.end(), B2.begin(), B2.end(), B.begin());
	
	// Remove duplicates from B.
	B.resize(distance(B.begin(), unique(B.begin(), B.end())));
	
	// Calculate travel times for each t \in B.
	vector<double> T;
	T.reserve(B.size());
	for (double t: B) T.push_back(travel_time(C, TZ, S, D, e, t));
	
	// Create travel time function.
	PWLFunction tau;
	for (int i = 0; i < (int)B.size()-1; ++i)
		tau.AddPiece(LinearFunction(Point2D(B[i], T[i]), Point2D(B[i+1], T[i+1])));
	
	return tau;
}
}

void preprocess_travel_times(json& instance)
{
	int n = instance["digraph"]["vertex_count"];
	double T = instance["speed_zones"].back()[1];
	int K = instance["speed_zones"].size();
	int CC = instance["cluster_speeds"].size();
	
	Digraph D = instance["digraph"];
	Matrix<int> C = instance["clusters"];
	vector<Interval> TZ = instance["speed_zones"];
	Matrix<double> S = instance["cluster_speeds"];
	Matrix<double> Dist = instance["distances"];
	Matrix<PWLFunction> tau(D.VertexCount(), D.VertexCount());
	for (Arc e: D.Arcs())
		tau[e.tail][e.head] = compute_travel_time_function(C, TZ, S, Dist, e);
	instance["travel_times"] = tau;
}
} // namespace tdtsptw