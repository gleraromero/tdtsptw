//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_TIME_DATE_H
#define GOC_TIME_DATE_H

#include <iostream>

#include "goc/print/printable.h"

namespace goc
{
// Represents a valid date in the gregorian calendar.
class Date : public Printable
{
public:
	Date();
	
	// Validates that the date is valid.
	Date(int day, int month, int year);
	
	int Day() const;
	
	int Month() const;
	
	int Year() const;
	
	// Prints the date in the dd-mm-yyyy format.
	virtual void Print(std::ostream& os) const;
	
private:
	int day_, month_, year_;
};
} // namespace goc

#endif //GOC_TIME_DATE_H
