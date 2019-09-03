//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_TIME_STOPWATCH_H
#define GOC_TIME_STOPWATCH_H

#include "goc/time/duration.h"

namespace goc
{
// Represents a device to measure how much time passed in an interval.
// - Uses CPP chrono steady clock.
class Stopwatch
{
public:
	// started indicates if the stopwatch starts turned on.
	Stopwatch(bool started = false);
	
	// Pauses the stopwatch but keeps track of the time spent.
	// Returns: the time spent so far.
	Duration Pause();
	
	// Resumes the stopwatch.
	// Returns: the time spent so far.
	Duration Resume();
	
	// Resets the time measured by the stopwatch to 0.
	// Returns: this object to concatenate calls.
	Stopwatch& Reset();
	
	// Returns: the time spent so far.
	Duration Peek() const;
	
private:
	mutable Duration partial_; // the time spent so far in the previously paused laps.
	mutable long long last_snapshot_; // the time the last lap was resumed.
	bool is_paused_; // indicates if the stopwatch is paused.
};
} // namespace goc

#endif //GOC_TIME_STOPWATCH_H
