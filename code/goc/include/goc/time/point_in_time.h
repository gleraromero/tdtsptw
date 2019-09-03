//
// Created by Gonzalo Lera Romero on 10/4/17.
//

#ifndef GOC_TIME_POINT_IN_TIME_H
#define GOC_TIME_POINT_IN_TIME_H

#include <iostream>

#include "goc/print/printable.h"
#include "goc/time/date.h"

namespace goc
{
// Represents a point in time.
// Example: 02-11-2019 03hs 04min 34secs.
class PointInTime : public Printable
{
public:
	PointInTime(const Date& date, int hour, int minute, int second);
	
	const goc::Date& Date() const;
	
	int Hour() const;
	
	int Minute() const;
	
	int Second() const;
	
	// Prints the point in time to the os stream.
	// Format: date hour:minute:seconds.
	virtual void Print(std::ostream& os) const;
	
private:
	goc::Date date_;
	int hour_, minute_, second_;
};
} // namespace goc

#endif //GOC_TIME_POINT_IN_TIME_H
