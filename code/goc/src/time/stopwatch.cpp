//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/time/stopwatch.h"

#include <chrono>

namespace goc
{
Stopwatch::Stopwatch(bool started)
{
	partial_ = 0.0_ms;
	is_paused_ = !started;
	if (!is_paused_)
	{
		last_snapshot_ = std::chrono::duration_cast<std::chrono::milliseconds>(
			std::chrono::steady_clock::now().time_since_epoch()).count();
	}
}

Duration Stopwatch::Pause()
{
	if (!is_paused_)
	{
		is_paused_ = true;
		long long now = std::chrono::duration_cast<std::chrono::milliseconds>(
			std::chrono::steady_clock::now().time_since_epoch()).count();
		partial_ += Duration(now - last_snapshot_, DurationUnit::Milliseconds);
	}
	return Peek();
}

Duration Stopwatch::Resume()
{
	Duration peek = Peek();
	if (is_paused_)
	{
		is_paused_ = false;
		last_snapshot_ = std::chrono::duration_cast<std::chrono::milliseconds>(
			std::chrono::steady_clock::now().time_since_epoch()).count();;
	}
	return peek;
}

Stopwatch& Stopwatch::Reset()
{
	is_paused_ = true;
	partial_ = 0.0_ms;
	return *this;
}

Duration Stopwatch::Peek() const
{
	if (!is_paused_)
	{
		long long now = std::chrono::duration_cast<std::chrono::milliseconds>(
			std::chrono::steady_clock::now().time_since_epoch()).count();
		partial_ += Duration(now - last_snapshot_, DurationUnit::Milliseconds);
		last_snapshot_ = now;
	}
	return partial_;
}
} // namespace goc