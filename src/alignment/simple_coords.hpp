/*
 * simple_coords.hpp
 *
 *  Created on: Nov 22, 2018
 *      Author: sarah
 */

#pragma once

#include <stddef.h>
#include <limits>

class SimpleCoords {
public:
	size_t first = std::numeric_limits<size_t>::max();
	size_t second = std::numeric_limits<size_t>::max();
	size_t leftGapSize = 0;
	size_t rightGapSize = 0;

	size_t size() const {
		if (first == std::numeric_limits<size_t>::max() || second == std::numeric_limits<size_t>::max() || first > second) {
			return 0;
		} else {
			return second + 1 - first;
		}
	}
};
