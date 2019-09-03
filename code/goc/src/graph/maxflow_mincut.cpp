//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/graph/maxflow_mincut.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_selectors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>

using namespace std;

namespace goc
{
namespace
{
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, boost::no_property, boost::property<boost::edge_index_t, std::size_t>> BoostDigraph;
typedef boost::graph_traits<BoostDigraph>::vertex_descriptor BoostVertex;
typedef BoostDigraph::edge_descriptor BoostArc;
}

pair<double, STCut> maxflow_mincut(const Digraph& D, const function<double(int i, int j)>& c, int s, int t)
{
	// n = number of vertices.
	int n = D.VertexCount();
	
	// Build boost network B to work with.
	BoostDigraph B;
	vector<BoostArc> reverse_arcs;
	vector<float> capacities;
	
	for (int i = 0; i < n; ++i)
	{
		for (int j: D.Successors(i))
		{
			boost::add_edge(i, j, boost::num_edges(B), B);
			capacities.push_back(c(i,j));
		}
	}
	
	// Add boost reverse arcs.
	vector<BoostArc> reverse_reverses; // reverse arcs of the reverse arcs.
	for (int i = 0; i < n; ++i)
	{
		for (int j: D.Successors(i))
		{
			auto reverse_arc = boost::edge(j,i,B);
			if (!reverse_arc.second) // If the arc was not in the network.
			{
				reverse_arc = boost::add_edge(j, i, boost::num_edges(B), B);
				capacities.push_back(0.0);
				reverse_reverses.push_back(boost::edge(i,j,B).first);
			}
			reverse_arcs.push_back(reverse_arc.first);
		}
	}
	for (auto e: reverse_reverses) reverse_arcs.push_back(e);
	
	vector<int> color(n);
	vector<float> residual_capacity(num_edges(B), 0);
	
	auto capacity_map = boost::make_iterator_property_map(&capacities[0], boost::get(boost::edge_index, B));
	auto residual_capacity_map = boost::make_iterator_property_map(&residual_capacity[0], boost::get(boost::edge_index, B));
	auto reverse_arc_map = boost::make_iterator_property_map(&reverse_arcs[0], boost::get(boost::edge_index, B));
	auto color_map = boost::make_iterator_property_map(&color[0], boost::get(boost::vertex_index, B));
	
	// Solve max-flow with boykov_kolmogorov algorithm.
	BoostVertex source = s;
	BoostVertex sink = t;
	double max_flow = boost::boykov_kolmogorov_max_flow(B, capacity_map, residual_capacity_map, reverse_arc_map,
														color_map, boost::get(boost::vertex_index, B), source, sink);
	
	// Get min-cut.
	STCut min_cut;
	for (int i = 0; i < n; ++i)
	{
		if (color[i] == boost::black_color) min_cut.S.push_back(i);
		else min_cut.T.push_back(i);
	}
	
	return {max_flow, min_cut};
}
} // namespace goc.