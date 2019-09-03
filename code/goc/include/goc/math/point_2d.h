//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_MATH_POINT_2D_H
#define GOC_MATH_POINT_2D_H

#include <iostream>

#include "goc/lib/json.hpp"
#include "goc/print/printable.h"

namespace goc
{
// This class represents a point (x,y) in a 2d cartesian plane.
class Point2D : public Printable
{
public:
	double x, y;
	
	Point2D(double x=0.0, double y=0.0);
	
	virtual void Print(std::ostream& os) const;
};

// JSON format: [x, y].
void from_json(const nlohmann::json& j, Point2D& p);

void to_json(nlohmann::json& j, const Point2D& p);
} // namespace goc

#endif //GOC_MATH_POINT_2D_H
