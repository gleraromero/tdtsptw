//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_TIME_WATCH_H
#define GOC_TIME_WATCH_H

#include "goc/time/date.h"
#include "goc/time/point_in_time.h"

namespace goc
{
// Represents a device that indicates the time and date at a given moment.
class Watch
{
public:
	Watch() = default;
	
	Date Today() const;
	
	PointInTime Now() const;
	
	long NowTicks() const;
};
} // namespace goc

#endif //GOC_TIME_WATCH_H
