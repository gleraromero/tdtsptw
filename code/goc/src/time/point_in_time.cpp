//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/time/point_in_time.h"

#include "goc/exception/exception_utils.h"

namespace goc
{
PointInTime::PointInTime(const class Date& date, int hour, int minute, int second) : date_(date)
{
	if (hour < 0 || hour > 23) fail("The hour must be between 0 and 23");
	if (minute < 0 || hour > 59) fail("The minute must be between 0 and 59");
	if (second < 0 || second > 59) fail("The second must be between 0 and 59");
	hour_ = hour;
	minute_ = minute;
	second_ = second;
}

const goc::Date& PointInTime::Date() const
{
	return date_;
}

int PointInTime::Hour() const
{
	return hour_;
}

int PointInTime::Minute() const
{
	return minute_;
}

int PointInTime::Second() const
{
	return second_;
}

void PointInTime::Print(std::ostream& os) const
{
	os << Date() << " " << Hour() << ":" << Minute() << ":" << Second();
}
} // namespace goc