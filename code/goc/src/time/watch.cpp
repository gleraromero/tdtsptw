//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/time/watch.h"

#include <ctime>

namespace goc
{
Date Watch::Today() const
{
	time_t t = time(NULL);
	tm* timePtr = localtime(&t);
	int day = timePtr->tm_mday;
	int month = timePtr->tm_mon + 1;
	int year = timePtr->tm_year + 1900;
	return Date(day, month, year);
}

PointInTime Watch::Now() const
{
	Date today = Today();
	time_t t = time(NULL);
	tm* timePtr = localtime(&t);
	int hour = timePtr->tm_hour;
	int minute = timePtr->tm_min;
	int second = timePtr->tm_sec;
	return PointInTime(today, hour, minute, second);
}

long Watch::NowTicks() const
{
	return clock();
}
} // namespace goc