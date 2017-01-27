/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <stdlib.h>
#include <limits>
#include <map>
#include <string>
#include <sstream>

/**
 * C++ version char* style "itoa": Convert integer to string
 */
template<typename T>
char* itoa10(const T& value, char* result) {
	// Check that base is valid
	char* out = result;
	T quotient = value;
	if(std::numeric_limits<T>::is_signed) {
		if(quotient <= 0) quotient = -quotient;
	}
	// Now write each digit from most to least significant
	do {
		*out = "0123456789"[quotient % 10];
		++out;
		quotient /= 10;
	} while (quotient > 0);
	// Only apply negative sign for base 10
	if(std::numeric_limits<T>::is_signed) {
		// Avoid compiler warning in cases where T is unsigned
		if (value <= 0 && value != 0) *out++ = '-';
	}
	reverse( result, out );
	*out = 0; // terminator
	return out;
}

// extract numeric ID from the beginning of a string
inline
uint64_t extractIDFromRefName(const string& refName) {
    uint64_t id = 0;
    for (size_t ni = 0; ni < refName.length(); ni++) {
        if (refName[ni] < '0' || refName[ni] > '9')
            break;

        id *= 10;
        id += (refName[ni] - '0');
    }
    return id;
}

// Converts a numeric value to std::string (part of C++11)
template <typename T>
std::string to_string(T value) {
 ostringstream ss;
 ss << value;
 return ss.str();
}

// Addition for pairs
template <typename T,typename U>
std::pair<T,U> operator+(const std::pair<T,U>& l,const std::pair<T,U>& r) {
    return {l.first+r.first,l.second+r.second};
}

// Addition for pairs
template <typename T>
std::array<T,2> operator+(const std::array<T,2>& l,const std::array<T,2>& r) {
    return {l[0]+r[0],l[1]+r[1]};
}

// Addition for pairs
template <typename T>
std::array<T,2>& operator+=(std::array<T,2>& l,const std::array<T,2>& r) {
    l[0] += r[0];
    l[1] += r[1];
    return l;
}

// Add second elements of pairs (e.g. values of maps)
struct AddValues {
  template<class Value, class Pair>
  Value operator()(Value value, const Pair& pair) const {
    return value + pair.second;
  }
};



// Find an element in a map or return a default value
template<typename K,typename V>
inline
V find_or_use_default(const std::map<K, V>& my_map, const K& query, const V default_value) {
	typedef typename std::map<K,V>::const_iterator MapIterator;
	MapIterator itr = my_map.find(query);

	if (itr == my_map.end()) {
		return default_value;
	}

	return itr->second;
}

// get value from a map, but return default if value is not in map
//  adapted from http://stackoverflow.com/questions/2333728/stdmap-default-value
template<typename MAP>
const typename MAP::mapped_type& get_with_default(const MAP& m,
                                             const typename MAP::key_type& key,
                                             const typename MAP::mapped_type& defval)
{
    typename MAP::const_iterator it = m.find(key);
    if (it == m.end())
        return defval;

    return it->second;
}


#endif /*ifndef UTIL_H_*/
