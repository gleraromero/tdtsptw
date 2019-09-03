//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_COLLECTION_COLLECTION_UTILS_H
#define GOC_COLLECTION_COLLECTION_UTILS_H

#include <algorithm>
#include <deque>
#include <functional>
#include <list>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>

namespace goc
{
// Returns: the collection without the elements in remove_set.
template<typename T>
std::vector<T> exclude(const std::vector<T>& collection, const std::set<T>& remove_set)
{
	std::vector<T> v;
	for (auto e: collection)
		if (remove_set.find(e) == remove_set.end())
			v.push_back(e);
	return v;
}

// Returns: a vector with all the numbers in [left, right).
// Observation: If left >= right then the returned vector is empty.
std::vector<int> range(int left, int right);

// Checks if the vector v has size at least 'size', otherwise it stretches the vector.
// It fills the new positions with the default element provided in the third argument 'default_element'.
template<typename T>
void stretch_to_size(std::vector<T>& v, int size, const T& default_element)
{
	while (v.size() < size) v.push_back(default_element);
}

// Returns: if the map 'm' contains the key: 'key'.
template<typename K, typename V>
bool includes_key(const std::map<K, V>& m, const K& key)
{
	return m.find(key) != m.end();
}

// Returns: if the map 'm' contains the key: 'key'.
template<typename K, typename V>
bool includes_key(const std::unordered_map<K, V>& m, const K& key)
{
	return m.find(key) != m.end();
}

// Returns: if the set 's' contains the element: 'element'.
template<typename T>
bool includes(const std::set<T>& s, const T& element)
{
	return s.find(element) != s.end();
}

// Returns: if the set 's' contains the element: 'element'.
template<typename T>
bool includes(const std::unordered_set<T>& s, const T& element)
{
	return s.find(element) != s.end();
}

// Returns: if the vector 'v' contains the element: 'element'.
template<typename T>
bool includes(const std::vector<T>& v, const T& element)
{
	return find(v.begin(), v.end(), element) != v.end();
}

// Returns: if the deque 'd' contains the element: 'element'.
template<typename T>
bool includes(const std::deque<T>& d, const T& element)
{
	return find(d.begin(), d.end(), element) != d.end();
}

// Returns: a new vector with the elements of 'v' in the opposite order.
template<typename T>
std::vector<T> reverse(const std::vector<T>& v)
{
	std::vector<T> r;
	r.reserve(v.size());
	for (int i = v.size()-1;i >= 0;--i) r.push_back(v[i]);
	return r;
}

// Returns: if the vector v1 is a prefix of the vector v2.
// We say (e1, e2, ..., ek) is a prefix of (e1, e2, ..., ek, ek+1, ...).
template<typename T>
bool is_prefix(const std::vector<T>& v1, const std::vector<T>& v)
{
	if (v.size() < v1.size()) return false;
	for (int i = 0; i < v1.size(); ++i) if (v[i] != v1[i]) return false;
	return true;
}

// Finds if 'element' is in 'v', if so, it returns a reference to that element in the vector, otherwise it adds it
// sorted.
// Precondition: T must have implemented the operator ==.
// ascendingly: indicates if elements in v are sorted ascendingly (e1 < e2 < ... < en).
template<typename T>
T& get_or_insert_sorted(const std::vector<T>& v, const T& element, bool ascendingly=true)
{
	// Find if the element exists in the vector, if so, return it.
	for (int i = 0; i < v.size(); ++i)
	{
		if (v[i] == element) return element;
		if (v[i] > element) break;
	}
	
	// Element does not exist, then add it.
	v.push_back(element);
	int i;
	for (i = (int)v.size()-1; i > 0; --i)
	{
		if (v[i] > v[i-1]) break;
		swap(v[i], v[i-1]);
	}
	return v[i];
}

// Adds item to the vector vec in its position using the predicate to compare.
// 	Pred(a, b): returns if item a should go before item b.
// Precondition: vec is sorted ascendingly.
template< typename T, typename Pred >
typename std::vector<T>::iterator insert_sorted( std::vector<T> & vec, T const& item, Pred pred )
{
	return vec.insert(std::upper_bound(vec.begin(), vec.end(), item, pred), item);
}

// Precondition: Type T must have the != operator implemented.
// Returns: if the two vectors have the same elements.
template<typename T>
bool operator==(const std::vector<T>& v1, const std::vector<T>& v2)
{
	if (v1.size() != v2.size()) return false;
	for (int i = 0; i < (int)v1.size(); ++i) if (v1[i] != v2[i]) return false;
	return true;
}

// Precondition: Type T must have the != operator implemented.
// Returns: if the two vectors have different elements.
template<typename T>
bool operator!=(const std::vector<T>& v1, const std::vector<T>& v2)
{
	return !(v1 == v2);
}

// Precondition: Type T must have the != operator implemented.
// Returns: if the two deques have the same elements.
template<typename T>
bool operator==(const std::deque<T>& d1, const std::deque<T>& d2)
{
	if (d1.size() != d2.size()) return false;
	for (int i = 0; i < (int)d1.size(); ++i) if (d1[i] != d2[i]) return false;
	return true;
}

// Precondition: Type T must have the != operator implemented.
// Returns: if the two deques have different elements.
template<typename T>
bool operator!=(const std::deque<T>& d1, const std::deque<T>& d2)
{
	return !(d1 == d2);
}

// Precondition: Type T must have the != operator implemented.
// Returns: if the two lists have the same elements.
template<typename T>
bool operator==(const std::list<T>& l1, const std::list<T>& l2)
{
	if (l1.size() != l2.size()) return false;
	for (auto it1 = l1.begin(), it2 = l2.begin(); it1 != l1.end(); ++it1, ++it2)
		if (*it1 != *it2)
			return false;
	return true;
}
} // namespace goc

#endif // GOC_COLLECTION_COLLECTION_UTILS_H