/*
 * dist_config.hpp
 *
 *  Created on: May 15, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <unordered_map>
#include <vector>
#include <cctype>
#include <stdexcept>
#include <utility>
#include <iostream>

class DistConfig {
public:
	DistConfig(int sub, int ins, int del) {
		substitution_penalty = sub;
		insertion_penalty = ins;
		deletion_penalty = del;
	}

	DistConfig(bool dnaData, bool proteinData) {
		this->dnaData = dnaData;
		this->proteinData = proteinData;
		if (dnaData && proteinData) {
			throw std::runtime_error("The data can't be both DNA and protein data!");
		}
		if (dnaData) {
			buildMappingsDNA();
			insertion_penalty = 1;
			deletion_penalty = 1;
		}
		if (proteinData) {
			buildBlosum62();
			// 10 is randomly chosen, TODO: Find a better value
			insertion_penalty = 10;
			deletion_penalty = 10;
		}
		substitution_penalty = 1;
	}

	DistConfig(bool dnaData, bool proteinData, int sub, int ins, int del) {
		this->dnaData = dnaData;
		this->proteinData = proteinData;
		if (dnaData && proteinData) {
			throw std::runtime_error("The data can't be both DNA and protein data!");
		}
		if (dnaData) {
			buildMappingsDNA();
		}
		if (proteinData) {
			buildBlosum62();
			// 10 is randomly chosen, TODO: Find a better value
		}
		insertion_penalty = ins;
		deletion_penalty = del;
		substitution_penalty = sub;
	}

	DistConfig() {}

	int getSubstitutionPenalty(char a, char b) const {
		if (dnaData) {
			if (a == b) {
				return 0;
			}

			if (mappings.find(a) == mappings.end()) {
				std::string aStr = "";
				aStr += a;
				throw std::runtime_error("Unsupported base a: " + aStr);
			}

			if (mappings.find(b) == mappings.end()) {
				std::string bStr = "";
				bStr += b;
				std::cout << a << ", " << b << "\n";
				throw std::runtime_error("Unsupported base b: " + bStr);
			}

			for (size_t i = 0; i < mappings.at(a).size(); ++i) {
				for (size_t j = 0; j < mappings.at(b).size(); ++j) {
					if (mappings.at(a)[i] == mappings.at(b)[j]) {
						return 0;
					}
				}
			}
			return substitution_penalty;
		} else if (proteinData) {
			return blosumMat[blosumMapping.at(a)][blosumMapping.at(b)];

		} else { // normal text
			if (a == b) {
				return 0;
			} else {
				return substitution_penalty;
			}
		}
	}
	int getInsertionPenalty() const {
		return 1;
	}
	int getDeletionPenalty() const {
		return 1;
	}
private:
	int substitution_penalty = 1;
	int insertion_penalty = 1;
	int deletion_penalty = 1;
	bool dnaData = false;
	bool proteinData = false;
	std::unordered_map<char, std::vector<char> > mappings;
	std::vector<std::vector<int> > blosumMat;
	std::unordered_map<char, size_t> blosumMapping;

	void buildMappingsDNA() {
		mappings['A'] = {'A'};
		mappings['C'] = {'C'};
		mappings['G'] = {'G'};
		mappings['T'] = {'T'};
		mappings['Y'] = {'C', 'T'};
		mappings['R'] = {'A', 'G'};
		mappings['W'] = {'A', 'T'};
		mappings['S'] = {'G', 'C'};
		mappings['K'] = {'T', 'G'};
		mappings['M'] = {'C', 'A'};
		mappings['D'] = {'A', 'G', 'T'};
		mappings['V'] = {'A', 'C', 'G'};
		mappings['H'] = {'A', 'C', 'T'};
		mappings['B'] = {'C', 'G', 'T'};
		mappings['X'] = {'A', 'C', 'G', 'T'};
		mappings['N'] = {'A', 'C', 'G', 'T'};
		mappings['$'] = {};
	}

	void buildBlosum62() {
		std::string letters = "ARNDCQEGHILKMFPSTWYVBZX*";
		for (size_t i = 0; i < letters.size(); ++i) {
			blosumMapping[letters[i]] = i;
		}
		// from https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
		blosumMat = {
			{	4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,-2,-1,0,-4},
			{	-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3,-1,0,-1,-4},
			{	-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,3,0,-1,-4},
			{	-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -4},
			{	0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
			{	-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1, -4},
			{	-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4},
			{	0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1, -4},
			{	-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -4},
			{	-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -4},
			{	-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1, -4},
			{	-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1, -4},
			{	-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -4},
			{	-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -4},
			{	-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
			{	1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0, -4},
			{	0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0, -4},
			{	-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2, -4},
			{	-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1, -4},
			{	0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1, -4},
			{	-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1, -4},
			{	-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4},
			{	0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -4},
			{	-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1}};
		// We need to adapt this matrix.
		for (size_t i = 0; i < blosumMat.size(); ++i) {
			for (size_t j = 0; j < blosumMat[i].size(); ++j) {
				// First, since we are trying to minimize a distance instead of maximizing similarity, we need to multiply all values by -1.
				blosumMat[i][j] *= -1;
				// Second, since we can only work with non-negative weights, we add 8 (largest number from the original matrix) to all entries.
				blosumMat[i][j] += 8;
			}
		}

	}
};
