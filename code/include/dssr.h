//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef TDTSPTW_DSSR_H
#define TDTSPTW_DSSR_H

#include <vector>
#include <tuple>

#include "goc/goc.h"

#include "vrp_instance.h"
#include "lbl_ng.h"

namespace tdtsptw
{
class TILabel : public goc::Printable
{
public:
	TILabel* prev;
	goc::Vertex v;
	VertexSet S;
	double Tdur;
	double Thelp;
	double lambda;
	
	TILabel(TILabel* prev, goc::Vertex v, VertexSet S, double Tdur, TimeUnit Thelp, double lambda);
	
	goc::GraphPath Path() const;
	
	virtual void Print(std::ostream& os) const;
};

class Label : public goc::Printable
{
public:
	Label* prev;
	goc::Vertex v;
	VertexSet S;
	goc::PWLFunction Tdur;
	TimeUnit Ttime;
	double lambda;
	
	Label(Label* prev, goc::Vertex v, VertexSet S, const goc::PWLFunction& Tdur, TimeUnit Ttime, double lambda);
	
	goc::GraphPath Path() const;
	
	virtual void Print(std::ostream& os) const;
	
	TILabel ToTI() const;
	
	inline double Thelp() const
	{
		return -(std::max(dom(Tdur))-Tdur(std::max(dom(Tdur))))-lambda;
	}
	
	inline double Tdurnum() const
	{
		return std::min(img(Tdur))-lambda;
	}
};

// Structure for obtaining lower bounds on the completion of all extensions.
class BoundingStructure
{
public:
	goc::Matrix<std::vector<goc::VectorMap<double, std::vector<Label>>>> S; // S[k][v][r][Thelp] sorted by min(img(Tdur)).
	VRPInstance* vrp;
	NGStructure* NG;
	std::vector<double> penalties;
	double UB;
	double penalties_sum;
	
	BoundingStructure(VRPInstance* vrp, NGStructure* NG, const std::vector<double>& penalties, double UB);
	
	void AddBound(int k, int r, const Label& l);
	
	double CompletionBound(int k, int r, const Label& l) const;
};

BoundingStructure run_dssr(const VRPInstance& vrp, const NGStructure& NG, int max_iter, const std::vector<double>& lambda,
	goc::Route& UB, double& LB, goc::CGExecutionLog* dssr_log, goc::MLBExecutionLog* exact_log);
} // namespace tdtsptw

#endif //TDTSPTW_DSSR_H
