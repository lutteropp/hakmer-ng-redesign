/*
 * dna_functions.hpp
 *
 *  Created on: Oct 25, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <stdexcept>

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

inline std::unordered_map<char, std::unordered_set<char> > createAmbiguityMappings() {
	std::unordered_map<char, std::unordered_set<char> > map;
	map['A'] = {'A', 'R', 'W', 'M', 'D', 'V', 'H', 'X', 'N', '-'};
	map['T'] = {'T', 'U', 'Y', 'W', 'K', 'D', 'H', 'B', 'X', 'N', '-'};
	map['U'] = {'T', 'U', 'Y', 'W', 'K', 'D', 'H', 'B', 'X', 'N', '-'};
	map['C'] = {'C', 'Y', 'S', 'M', 'V', 'H', 'B', 'X', 'N', '-'};
	map['G'] = {'G', 'R', 'S', 'K', 'D', 'V', 'B', 'X', 'N', '-'};
	map['Y'] = {'C', 'T', 'U', 'Y', 'W', 'S', 'K', 'D', 'V', 'H', 'B', 'X', 'N', '-'};
	map['R'] = {'A', 'G', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B', 'X', 'N', '-'};
	map['W'] = {'A', 'T', 'U', 'Y', 'R', 'W', 'K', 'M', 'D', 'V', 'H', 'B', 'X', 'N', '-'};
	map['S'] = {'G', 'C', 'Y', 'R', 'S', 'K', 'M', 'D', 'V', 'H', 'B', 'X', 'N', '-'};
	map['K'] = {'T', 'U', 'G', 'Y', 'R', 'W', 'S', 'K', 'D', 'V', 'H', 'B', 'X', 'N', '-'};
	map['M'] = {'C', 'A', 'Y', 'R', 'W', 'S', 'M', 'D', 'V', 'H', 'B', 'X', 'N', '-'};
	map['D'] = {'A', 'G', 'T', 'U', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B', 'X', 'N', '-'};
	map['V'] = {'A', 'C', 'G', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B', 'X', 'N', '-'};
	map['H'] = {'A', 'C', 'T', 'U', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B', 'X', 'N', '-'};
	map['B'] = {'G', 'C', 'T', 'U', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B', 'X', 'N', '-'};
	map['X'] = {'A', 'G', 'C', 'T', 'U', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B', 'X', 'N', '-'};
	map['N'] = {'A', 'C', 'G', 'T', 'U', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B', 'X', 'N', '-'};
	map['-'] = {'A', 'C', 'G', 'T', 'U', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B', 'X', 'N', '-'};
	return map;
}

static std::unordered_map<char, char> rcMapping = createRevcompMapping();
static std::unordered_map<char, std::unordered_set<char> > ambiMapping = createAmbiguityMappings();

inline bool ambiguousMatch(char a, char b) {
	if (a == b) {
		return true;
	} else {
		return (ambiMapping[a].find(b) != ambiMapping[a].end());
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
