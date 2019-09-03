//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/math/point_2d.h"

using namespace std;
using namespace nlohmann;

namespace goc
{
Point2D::Point2D(double x, double y) : x(x), y(y)
{

}

void Point2D::Print(ostream& os) const
{
	os << "(" << x << "," << y << ")";
}

void from_json(const json& j, Point2D& p)
{
	p.x = j[0];
	p.y = j[1];
}

void to_json(json& j, const Point2D& p)
{
	j = vector<double>({p.x, p.y});
}
} // namespace goc