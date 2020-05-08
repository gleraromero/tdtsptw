//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/vrp/route.h"

using namespace std;
using namespace nlohmann;

namespace goc
{
Route::Route() : t0(0), duration(0)
{ }

Route::Route(const GraphPath& path, double t0, double duration)
	: path(path), t0(t0), duration(duration)
{ }

bool Route::IsElementary() const
{
	vector<int> repetitions(*max_element(path.begin(), path.end())+1, 0);
	for (Vertex v: path) if (repetitions[v]++ > 0) return false;
	return true;
}

void Route::Print(ostream& os) const
{
	os << json(*this);
}

void to_json(json& j, const Route& r)
{
	j["path"] = r.path;
	j["t0"] = r.t0;
	j["duration"] = r.duration;
}

void from_json(const json& j, Route& r)
{
	vector<int> p = j["path"];
	r.path = p;
	r.t0 = j["t0"];
	r.duration = j["duration"];
}
} // namespace goc