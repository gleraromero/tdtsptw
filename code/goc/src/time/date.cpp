//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/time/date.h"

#include <vector>

#include "goc/exception/exception_utils.h"

namespace goc
{
Date::Date()
{

}

Date::Date(int day, int month, int year)
{
	bool is_leap_year = (year % 4 == 0) && (year % 100 != 0 || year % 400 == 0);
	std::vector<int> days_in_month = {31, is_leap_year ? 29 : 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	if (year < 0) fail("Year must be at least 0");
	if (month < 1 || month > 12) fail("Month must be between 1 and 12");
	if (day < 1 || day > days_in_month[month - 1]) fail("Day must be in month's days");
	day_ = day;
	month_ = month;
	year_ = year;
}

int Date::Day() const
{
	return day_;
}

int Date::Month() const
{
	return month_;
}

int Date::Year() const
{
	return year_;
}

void Date::Print(std::ostream& os) const
{
	os << Day() << "-" << Month() << "-" << Year();
}
} // namespace goc