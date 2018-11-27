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
	std::vector<std::pair<size_t, size_t> > extraOccs;
	bool operator <(const SeedInfo& str) const {
		if (n == str.n) {
			if (extraOccs.size() == str.extraOccs.size()) {
				return k < str.k;
			} else {
				return extraOccs.size() < str.extraOccs.size();
			}
		} else {
			return n < str.n;
		}
	}
	bool operator >(const SeedInfo& str) const {
		if (n == str.n) {
			if (extraOccs.size() == str.extraOccs.size()) {
				return k > str.k;
			} else {
				return extraOccs.size() > str.extraOccs.size();
			}
		} else {
			return n > str.n;
		}
	}
	SeedInfo(size_t saPos, unsigned int k, unsigned int n) {
		this->saPos = saPos;
		this->k = k;
		this->n = n;
	}
	void addExtraOcc(const std::pair<size_t, size_t>& occ) {
		extraOccs.push_back(occ);
	}
};
