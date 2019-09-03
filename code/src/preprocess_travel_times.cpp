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
double departing_time(const json& instance, Arc e, double tf)
{
	int c = instance["clusters"][e.tail][e.head]; // cluster of arc e.
	vector<Interval> T = instance["speed_zones"]; // T[k] = speed zone k.
	vector<double> speed = instance["cluster_speeds"][c]; // speed[k] = speed of traversing e in speed zone k.
	double d = instance["distances"][e.tail][e.head]; // distance of arc e.
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
double travel_time(const json& instance, Arc e, double t0)
{
	int c = instance["clusters"][e.tail][e.head]; // cluster of arc e.
	vector<Interval> T = instance["speed_zones"]; // T[k] = speed zone k.
	vector<double> speed = instance["cluster_speeds"][c]; // speed[k] = speed of traversing e in speed zone k.
	double d = instance["distances"][e.tail][e.head]; // distance of arc e.
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

// Returns the time when we arrive at the end of arc e if departing at t0.
double ready_time(const json& instance, Arc e, double t0)
{
	double tt = travel_time(instance, e, t0);
	return tt == INFTY ? tt : t0 + tt;
}

// Precondition: no speeds are 0.
PWLFunction compute_travel_time_function(const json& instance, Arc e)
{
	// Calculate speed breakpoints.
	vector<Interval> speed_zones = instance["speed_zones"];
	vector<double> speed_breakpoints;
	for (auto& z: speed_zones) speed_breakpoints.push_back(z.left);
	speed_breakpoints.push_back(speed_zones.back().right);
	
	// Travel time breakpoints are two sets
	// 	- B1: speed breakpoints which are feasible to depart
	// 	- B2: times t such that we arrive to head(e) at a speed breakpoint.
	vector<double> B1;
	for (double t: speed_breakpoints)
		if (travel_time(instance, e, t) != INFTY)
			B1.push_back(t);
		
	vector<double> B2;
	for (double t: speed_breakpoints)
		if (departing_time(instance, e, t) != INFTY)
			B2.push_back(departing_time(instance, e, t));
	
	// Merge breakpoints in order in a set B.
	vector<double> B(B1.size()+B2.size());
	merge(B1.begin(), B1.end(), B2.begin(), B2.end(), B.begin());
	
	// Remove duplicates from B.
	B.resize(distance(B.begin(), unique(B.begin(), B.end())));
	
	// Calculate travel times for each t \in B.
	vector<double> T;
	for (double t: B) T.push_back(travel_time(instance, e, t));
	
	// Create travel time function.
	PWLFunction tau;
	for (int i = 0; i < (int)B.size()-1; ++i)
		tau.AddPiece(LinearFunction(Point2D(B[i], T[i]), Point2D(B[i+1], T[i+1])));
	
	return tau;
}
}

void preprocess_travel_times(json& instance)
{
	Digraph D = instance["digraph"];
	Matrix<PWLFunction> tau(D.VertexCount(), D.VertexCount());
	for (Arc e: D.Arcs())
		tau[e.tail][e.head] = compute_travel_time_function(instance, e);
	instance["travel_times"] = tau;
}
} // namespace tdtsptw