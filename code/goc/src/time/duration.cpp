//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/time/duration.h"

#include "goc/exception/exception_utils.h"

using namespace std;
using namespace nlohmann;

namespace goc
{
Duration Duration::None()
{
	return 0.0_sec;
}

Duration Duration::Max()
{
	return 1000.0_hr;
}

Duration::Duration()
{
	amount_in_ms_ = 0;
}

Duration::Duration(double amount, DurationUnit unit)
{
	amount_in_ms_ = ConvertToMilliseconds(amount, unit);
}

double Duration::Amount(DurationUnit unit) const
{
	return ConvertFromMilliseconds(amount_in_ms_, unit);
}

Duration Duration::operator+(const Duration& d) const
{
	return Duration(Amount(DurationUnit::Milliseconds) + d.Amount(DurationUnit::Milliseconds),
					DurationUnit::Milliseconds);
}

Duration Duration::operator-(const Duration& d) const
{
	return Duration(Amount(DurationUnit::Milliseconds) - d.Amount(DurationUnit::Milliseconds),
					DurationUnit::Milliseconds);
}

bool Duration::operator<(const Duration& d) const
{
	return Amount(DurationUnit::Milliseconds) < d.Amount(DurationUnit::Milliseconds);
}

bool Duration::operator>(const Duration& d) const
{
	return Amount(DurationUnit::Milliseconds) > d.Amount(DurationUnit::Milliseconds);
}

bool Duration::operator<=(const Duration& d) const
{
	return Amount(DurationUnit::Milliseconds) <= d.Amount(DurationUnit::Milliseconds);
}

bool Duration::operator>=(const Duration& d) const
{
	return Amount(DurationUnit::Milliseconds) >= d.Amount(DurationUnit::Milliseconds);
}

bool Duration::operator==(const Duration& d) const
{
	return Amount(DurationUnit::Milliseconds) == d.Amount(DurationUnit::Milliseconds);
}

bool Duration::operator!=(const Duration& d) const
{
	return Amount(DurationUnit::Milliseconds) != d.Amount(DurationUnit::Milliseconds);
}

Duration& Duration::operator+=(const Duration& d)
{
	amount_in_ms_ = amount_in_ms_ + d.amount_in_ms_;
	return *this;
}

Duration& Duration::operator-=(const Duration& d)
{
	amount_in_ms_ = amount_in_ms_ - d.amount_in_ms_;
	return *this;
}

void Duration::Print(ostream& os) const
{
	os << Amount(DurationUnit::Seconds);
}

double Duration::ConvertFromMilliseconds(double duration_ms, DurationUnit to_unit) const
{
	switch (to_unit)
	{
		case DurationUnit::Milliseconds:
			return duration_ms;
		case DurationUnit::Seconds:
			return duration_ms / 1000.0;
		case DurationUnit::Minutes:
			return duration_ms / 60000.0;
		case DurationUnit::Hours:
			return duration_ms / 3600000.0;
		default:
			fail("Unknown duration unit");
	}
	return -1;
}

double Duration::ConvertToMilliseconds(double duration, DurationUnit from_unit) const
{
	switch (from_unit)
	{
		case DurationUnit::Milliseconds:
			return duration;
		case DurationUnit::Seconds:
			return duration * 1000.0;
		case DurationUnit::Minutes:
			return duration * 60000.0;
		case DurationUnit::Hours:
			return duration * 3600000.0;
		default:
			fail("Unknown duration unit");
	}
	return -1;
}

Duration operator "" _ms(long double time)
{
	return Duration(time, DurationUnit::Milliseconds);
}

Duration operator "" _sec(long double time)
{
	return Duration(time, DurationUnit::Seconds);
}

Duration operator "" _min(long double time)
{
	return Duration(time, DurationUnit::Minutes);
}

Duration operator "" _hr(long double time)
{
	return Duration(time, DurationUnit::Hours);
}

void to_json(json& j, const Duration& d)
{
	json x = {};
	x["number"] = d.Amount(DurationUnit::Seconds);
	j = x["number"];
}

void from_json(const nlohmann::json& j, Duration& d)
{
	d = Duration(j, DurationUnit::Seconds);
}

} // namespace goc