//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef NETWORKS2019_VECTOR_MAP_H
#define NETWORKS2019_VECTOR_MAP_H

#include <iostream>
#include <vector>

#include "goc/print/printable.h"

namespace goc
{
// This class implements a dictionary over a vector. It stores (key, value) pairs sorted by key.
// Invariant: there is only one entry per key.
// Invariant: the vector is sorted by K ascendingly.
template<typename K, typename V>
class VectorMap : public Printable
{
public:
	// Sets the value to the key if the key is not present.
	// Returns: a reference to the value associated to the key.
	V& Insert(const K& key, const V& value)
	{
		// Try to see if the element already existed.
		for (auto& e: S)
		{
			if (e.first == key) return e.second;
			else if (e.first > key) break;
		}
		
		// Element does not exist, then insert.
		S.push_back(make_pair(key, value));
		int i = S.size()-1;
		while (i > 0 && S[i-1].first > S[i].first)
		{
			swap(S[i-1], S[i]);
			--i;
		}
		return S[i].second;
	}
	
	// Returns: the (key, value) pair indexed with the index when sorted.
	std::pair<K,V>& operator[](int index)
	{
		return S[index];
	}
	
	// Returns: the (key, value) pair indexed with the index when sorted.
	const std::pair<K,V>& operator[](int index) const
	{
		return S[index];
	}
	
	// Removes the (key, value) pair at position index.
	// Returns: an index to the previous element if existent, otherwise to the first element.
	int Erase(int index)
	{
		S.erase(index);
		return std::max(0, index-1);
	}
	
	typename std::vector<std::pair<K, V>>::iterator begin()
	{
		return S.begin();
	}
	
	typename std::vector<std::pair<K, V>>::const_iterator begin() const
	{
		return S.begin();
	}
	
	typename std::vector<std::pair<K, V>>::iterator end()
	{
		return S.end();
	}
	
	typename std::vector<std::pair<K, V>>::const_iterator end() const
	{
		return S.end();
	}
	
	typename std::vector<std::pair<K, V>>::reverse_iterator rbegin()
	{
		return S.rbegin();
	}
	
	typename std::vector<std::pair<K, V>>::const_reverse_iterator rbegin() const
	{
		return S.rbegin();
	}
	
	typename std::vector<std::pair<K, V>>::reverse_iterator rend()
	{
		return S.rend();
	}
	
	typename std::vector<std::pair<K, V>>::const_reverse_iterator rend() const
	{
		return S.rend();
	}
	
	// Prints the (key, value) pair sequence sorted by key.
	virtual void Print(std::ostream& os) const
	{
		os << S;
	}
	
private:
	std::vector<std::pair<K, V>> S;
};
} // namespace goc

#endif //NETWORKS2019_VECTOR_MAP_H
