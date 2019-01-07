/*
 * dna_functions.hpp
 *
 *  Created on: Oct 25, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <stddef.h>
#include <cassert>
#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>

inline std::vector<char> createRevcompMapping() {
	std::vector<char> map(128, '#');
	map[static_cast<unsigned char>('A')] = 'T';
	map[static_cast<unsigned char>('T')] = 'A';
	map[static_cast<unsigned char>('U')] = 'U';
	map[static_cast<unsigned char>('C')] = 'G';
	map[static_cast<unsigned char>('G')] = 'C';
	map[static_cast<unsigned char>('N')] = 'N';
	// DNA ambiguity codes, see http://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
	map[static_cast<unsigned char>('Y')] = 'R';
	map[static_cast<unsigned char>('R')] = 'Y';
	map[static_cast<unsigned char>('W')] = 'W';
	map[static_cast<unsigned char>('S')] = 'S';
	map[static_cast<unsigned char>('K')] = 'M';
	map[static_cast<unsigned char>('M')] = 'K';
	map[static_cast<unsigned char>('D')] = 'H';
	map[static_cast<unsigned char>('V')] = 'B';
	map[static_cast<unsigned char>('H')] = 'D';
	map[static_cast<unsigned char>('B')] = 'V';
	map[static_cast<unsigned char>('X')] = 'X';
	// delimiters stay unchanged
	map[static_cast<unsigned char>('$')] = '$';
	// gaps stay unchanged
	map[static_cast<unsigned char>('-')] = '-';
	return map;
}

inline std::vector<unsigned short int> createAmbiguityMappings() {
	std::vector<unsigned short int> map(128, 0);
	map[static_cast<unsigned char>('A')] = 0b1000;
	map[static_cast<unsigned char>('T')] = 0b0001;
	map[static_cast<unsigned char>('U')] = 0b0001;
	map[static_cast<unsigned char>('C')] = 0b0100;
	map[static_cast<unsigned char>('G')] = 0b0010;
	map[static_cast<unsigned char>('Y')] = 0b0101;
	map[static_cast<unsigned char>('R')] = 0b1010;
	map[static_cast<unsigned char>('W')] = 0b1001;
	map[static_cast<unsigned char>('S')] = 0b0110;
	map[static_cast<unsigned char>('K')] = 0b0011;
	map[static_cast<unsigned char>('M')] = 0b1100;
	map[static_cast<unsigned char>('D')] = 0b1011;
	map[static_cast<unsigned char>('V')] = 0b1110;
	map[static_cast<unsigned char>('H')] = 0b1101;
	map[static_cast<unsigned char>('B')] = 0b0111;
	map[static_cast<unsigned char>('X')] = 0b1111;
	map[static_cast<unsigned char>('N')] = 0b1111;
	return map;
}

static std::vector<char> rcMapping = createRevcompMapping();
static std::vector<unsigned short int> ambiMapping = createAmbiguityMappings();

inline bool ambiguousMatch(char a, char b) {
	if (a == '$' || b == '$') {
		return false;
	}
	if (a == b) {
		return true;
	} else {
		return (ambiMapping[static_cast<unsigned char>(a)] & ambiMapping[static_cast<unsigned char>(b)] != 0);
	}
}

inline bool ambiguousEqual(const std::string& s1, const std::string& s2) {
	if (s1.size() != s2.size()) {
		return false;
	}
	for (size_t i = 0; i < s1.size(); ++i) {
		if (!ambiguousMatch(s1[i], s2[i])) {
			return false;
		}
	}
	return true;
}

inline std::string revComp(const std::string& str) {
	std::string res;
	for (int i = str.size() - 1; i >= 0; --i) {
		if (static_cast<unsigned char>(str[i]) > rcMapping.size()) {
			std::cout << str[i] << "\n";
			throw std::runtime_error("Problematic base encountered!");
		}

		if (rcMapping[static_cast<unsigned char>(str[i])] == '#') {
			std::cout << str[i] << "\n";
			throw std::runtime_error("Problematic base encountered 2!");
		}

		assert(rcMapping[static_cast<unsigned char>(str[i])] != '#');
		res += rcMapping[static_cast<unsigned char>(str[i])];
	}
	return res;
}
