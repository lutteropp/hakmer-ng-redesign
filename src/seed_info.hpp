/*
 * SeedInfo.hpp
 *
 *  Created on: Nov 27, 2018
 *      Author: sarah
 */

#pragma once

#include <stddef.h>

class SeedInfo {
public:
	size_t saPos;
	unsigned int k;
	unsigned int n;
	bool operator <(const SeedInfo& str) const {
		if (n == str.n) {
			return k < str.k;
		} else {
			return n < str.n;
		}
	}
	bool operator >(const SeedInfo& str) const {
		if (n == str.n) {
			return k > str.k;
		} else {
			return n > str.n;
		}
	}
	SeedInfo() {
		saPos = 0;
		k = 0;
		n = 0;
	}
	SeedInfo(size_t saPos, unsigned int k, unsigned int n) {
		this->saPos = saPos;
		this->k = k;
		this->n = n;
	}
};
