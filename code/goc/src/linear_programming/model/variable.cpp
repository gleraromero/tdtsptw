//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/linear_programming/model/variable.h"

using namespace std;

namespace goc
{
Variable::Variable() : index_(nullptr)
{ }

int Variable::Index() const
{
	return index_ == nullptr ? -1 : *index_;
}

void Variable::Print(ostream& os) const
{
	os << name;
}

bool Variable::operator<(const Variable& v) const
{
	return index_ < v.index_;
}

bool Variable::operator==(const Variable& v) const
{
	return v.index_ == index_;
}

bool Variable::operator!=(const Variable& v) const
{
	return v.index_ != index_;
}

Variable::Variable(const string& name, int* index)
	: name(name), index_(index)
{ }
} // namespace goc