/*
 * pairwise_alignment.hpp
 *
 *  Created on: Oct 9, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <string>
#include "twodim_matrix.hpp"

const int GAP_PENALTY = 1;
const int MISMATCH_PENALTY = 1;

// Still an open question: How to select the alignment penalties (indels, substitutions)?
// Open thing to do: Improve runtime and space-requirements of the pairwise alignments...
// Maybe we could use edlib for very performant pairwise alignment (but it only supports Levenshtein distance)?
// TODO: Maybe also support affine gap penalties?

class PairwiseAlignment {
public:
	PairwiseAlignment();
	std::string getS1() const;
	std::string getS2() const;
	std::pair<std::string, std::string> extractAlignment();
	void printAlignment();
	void addChars(char a, char b);
	void addCharS1(char c);
	void addCharS2(char c);
	double pairwiseDistance();
	size_t getAlignmentWidth();
	void shrinkDownTo(size_t newS1Size, size_t newS2Size);
private:
	void backtrack(size_t i, size_t j, std::string& s1Aligned, std::string& s2Aligned);
	void update(size_t i, size_t j);
	std::string s1;
	std::string s2;
	TwoDimMatrix<int> matrix;
	std::pair<std::string, std::string> ali;
	bool aliValid;
};
