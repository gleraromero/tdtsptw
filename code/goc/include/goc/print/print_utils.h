//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

// In this file we include functions for printing objects using the << operator, so they do not have to be reimplemented
// over and over again.

#include <bitset>
#include <deque>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifndef GOC_PRINT_PRINT_UTILS_H
#define GOC_PRINT_PRINT_UTILS_H

namespace goc
{
// Prints a sequence of elements in the ostream os.
// The output format is "prefix element1, element2, ..., elementn suffix".
// Precondition: The iterated elements must have the <<(ostream&) operator implemented.
// Returns: the modified ostream.
template<typename Collection>
std::ostream& print_iterable(std::ostream& os, const Collection& collection, char prefix='{', char suffix='}')
{
	os << prefix;
	for (auto it = collection.begin(); it != collection.end(); ++it)
	{
		if (it != collection.begin()) os << ", ";
		os << *it;
	}
	os << suffix;
	return os;
}

// Prints the pair in the ostream os.
// The output format is "(first, second)".
// Precondition: The types T1, T2 must have the <<(ostream&) operator implemented.
// Returns: the modified ostream.
template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1, T2>& p)
{
	return os << "(" << p.first << ", " << p.second << ")";
}

// Prints the set in the ostream os.
// The output format is "{ element1, element2, ..., elementn }".
// Precondition: The type T must have the <<(ostream&) operator implemented.
// Returns: the modified ostream.
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::set<T>& s)
{
	return print_iterable(os, s);
}

// Prints the set in the ostream os.
// The output format is "{ element1, element2, ..., elementn }".
// Precondition: The type T must have the <<(ostream&) operator implemented.
// Returns: the modified ostream.
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::unordered_set<T>& s)
{
	return print_iterable(os, s);
}

// Prints the map in the ostream os.
// The output format is "{ (key1, value1), (key2, value2), ..., (keyn, valuen) }".
// Precondition: The types K, V must have the <<(ostream&) operator implemented.
// Returns: the modified ostream.
template<typename K, typename V>
std::ostream& operator<<(std::ostream& os, const std::map<K, V>& m)
{
	return print_iterable(os, m);
}

// Prints the map in the ostream os.
// The output format is "{ (key1, value1), (key2, value2), ..., (keyn, valuen) }".
// Precondition: The types K, V must have the <<(ostream&) operator implemented.
// Returns: the modified ostream.
template<typename K, typename V>
std::ostream& operator<<(std::ostream& os, const std::unordered_map<K, V>& m)
{
	return print_iterable(os, m);
}

// Prints the vector in the ostream os.
// The output format is "[ element1, element2, ..., elementn ]".
// Precondition: The type T must have the <<(ostream&) operator implemented.
// Returns: the modified ostream.
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
	os << '[';
	for (auto it = v.begin(); it != v.end(); ++it)
	{
		if (it != v.begin()) os << ", ";
		os << *it;
	}
	os << ']';
	return os;
}

// Prints the deque in the ostream os.
// The output format is "[ element1, element2, ..., elementn ]".
// Precondition: The type T must have the <<(ostream&) operator implemented.
// Returns: the modified ostream.
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::deque<T>& d)
{
	return print_iterable(os, d, '[', ']');
}

// Prints the list in the ostream os.
// The output format is "[ element1, element2, ..., elementn ]".
// Precondition: The type T must have the <<(ostream&) operator implemented.
// Returns: the modified ostream.
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::list<T>& l)
{
	return print_iterable(os, l, '[', ']');
}

// Prints the bitset in the ostream os.
// The output format is "( n1, n2, ..., nk )" where n1, ..., nk are set in the bitset.
// Returns: the modified ostream.
template<unsigned long N>
std::ostream& operator<<(std::ostream& os, const std::bitset<N>& b)
{
	os << "(";
	bool first_added = true;
	for (unsigned long i = 0; i < N; ++i)
	{
		if (b.test(i))
		{
			if (!first_added) os << ", ";
			first_added = false;
			os << i;
		}
	}
	return os << ")";
}

} // namespace goc

#endif //GOC_PRINT_PRINT_UTILS_H
