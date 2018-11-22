/*
 * simple_msa.hpp
 *
 *  Created on: Nov 14, 2018
 *      Author: sarah
 */

#pragma once

#include "msa_wrapper.hpp"
#include <string>

struct SimpleCoords {
public:
	size_t first = std::string::npos;
	size_t second = std::string::npos;
	size_t leftGapSize = 0;
	size_t rightGapSize = 0;

	size_t size() const {
		if (first == std::string::npos || second == std::string::npos || first > second) {
			return 0;
		} else {
			return second + 1 - first;
		}
	}
};

std::vector<std::string> computeMSA(const std::vector<std::string>& seqs);
std::vector<std::string> computeMSA(const std::vector<SimpleCoords>& seqCoords, const std::string& T);
