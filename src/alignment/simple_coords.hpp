/*
 * simple_coords.hpp
 *
 *  Created on: Nov 22, 2018
 *      Author: sarah
 */

#pragma once

#include <stddef.h>
#include <string>

class SimpleCoords {
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
