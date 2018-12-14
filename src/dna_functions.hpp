/*
 * dna_functions.hpp
 *
 *  Created on: Oct 25, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <stddef.h>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

inline std::unordered_map<char, char> createRevcompMapping() {
	std::unordered_map<char, char> mapping; // TODO: Make this one static somehow
	mapping['A'] = 'T';
	mapping['T'] = 'A';
	mapping['U'] = 'U';
	mapping['C'] = 'G';
	mapping['G'] = 'C';
	mapping['N'] = 'N';
	// DNA ambiguity codes, see http://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
	mapping['Y'] = 'R';
	mapping['R'] = 'Y';
	mapping['W'] = 'W';
	mapping['S'] = 'S';
	mapping['K'] = 'M';
	mapping['M'] = 'K';
	mapping['D'] = 'H';
	mapping['V'] = 'B';
	mapping['H'] = 'D';
	mapping['B'] = 'V';
	mapping['X'] = 'X';
	// delimiters stay unchanged
	mapping['$'] = '$';
	// gaps stay unchanged
	mapping['-'] = '-';
	return mapping;
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

static std::unordered_map<char, char> rcMapping = createRevcompMapping();
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
		if (rcMapping.find(str[i]) != rcMapping.end()) {
			res += rcMapping[str[i]];
		} else {
			std::string c;
			c += str[i];
			std::cout << "i: " << i << "str.size(): " << str.size() << "\n";
			throw std::runtime_error("Cannot reverse-complement the following base: " + c);
		}
	}
	return res;
}
