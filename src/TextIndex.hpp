/*
 * SuffixArray.hpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <vector>
#include <string>

class TextIndex {
public:
	TextIndex(const std::string& text);
	size_t countMatches(const std::string& pattern);
	std::vector<size_t> getMatches(const std::string& pattern);
};
