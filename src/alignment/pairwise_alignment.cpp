/*
 * pairwise_alignment.cpp
 *
 *  Created on: Oct 9, 2018
 *      Author: Sarah Lutteropp
 */

#include <iostream>
#include <algorithm>
#include "pairwise_alignment.hpp"
#include "../dna_functions.hpp"

PairwiseAlignment::PairwiseAlignment() :
		s1(""), s2(""), aliValid(false) {
	matrix.addRow();
	matrix.addColumn();
	// now we have a 1x1 matrix
	matrix.entryAt(0, 0) = 0;
}
std::string PairwiseAlignment::getS1() const {
	return s1;
}
std::string PairwiseAlignment::getS2() const {
	return s2;
}
std::pair<std::string, std::string> PairwiseAlignment::extractAlignment() {
	if (aliValid) {
		return ali;
	}
	std::string s1Aligned;
	std::string s2Aligned;
	s1Aligned.reserve(s1.size());
	s2Aligned.reserve(s2.size());
	backtrack(s1.size(), s2.size(), s1Aligned, s2Aligned);
	std::reverse(s1Aligned.begin(), s1Aligned.end());
	std::reverse(s2Aligned.begin(), s2Aligned.end());

	ali = std::make_pair(s1Aligned, s2Aligned);
	aliValid = true;
	return ali;
}
void PairwiseAlignment::printAlignment() {
	if (!aliValid) {
		extractAlignment();
	}
	std::cout << ali.first << "\n" << ali.second << "\n\n";
}
void PairwiseAlignment::addChars(char a, char b) {
	matrix.addRow();
	matrix.addColumn();
	s1 += a;
	s2 += b;
	// update the new row (except for last column)
	for (size_t j = 0; j < matrix.getM() - 1; ++j) {
		update(matrix.getN() - 1, j);
	}
	// update the new column (except for last row)
	for (size_t i = 0; i < matrix.getN() - 1; ++i) {
		update(i, matrix.getM() - 1);
	}
	// update entry at (last row, last column)
	update(matrix.getN() - 1, matrix.getM() - 1);
	aliValid = false;
}

void PairwiseAlignment::addCharS1(char c) {
	matrix.addRow();
	s1 += c;
	// update the new row
	for (size_t j = 0; j < matrix.getM(); ++j) {
		update(matrix.getN() - 1, j);
	}
	aliValid = false;
}
void PairwiseAlignment::addCharS2(char c) {
	matrix.addColumn();
	s2 += c;
	// update the new column
	for (size_t i = 0; i < matrix.getN(); ++i) {
		update(i, matrix.getM() - 1);
	}
	aliValid = false;
}

double PairwiseAlignment::pairwiseDistance() {
	return matrix.entryAt(s1.size(), s2.size());
}
size_t PairwiseAlignment::getAlignmentWidth() {
	std::pair<std::string, std::string> ali = extractAlignment();
	return ali.first.size();
}
void PairwiseAlignment::shrinkDownTo(size_t newS1Size, size_t newS2Size) {
	s1.resize(newS1Size);
	s1.shrink_to_fit();
	s2.resize(newS2Size);
	s2.shrink_to_fit();
	matrix.shrinkDownTo(newS1Size + 1, newS2Size + 1);
	aliValid = false;
}

void PairwiseAlignment::backtrack(size_t i, size_t j, std::string& s1Aligned, std::string& s2Aligned) {
	if (i == 0 || j == 0) { // backtracking has ended
	// add missing gaps
		std::string missingGaps1 = "";
		std::string missingText2 = "";
		for (size_t j1 = 0; j1 < j; ++j1) {
			missingGaps1 += "-";
			missingText2 += s2[j - j1 - 1];
		}
		std::string missingGaps2 = "";
		std::string missingText1 = "";
		for (size_t i1 = 0; i1 < i; ++i1) {
			missingGaps2 += "-";
			missingText1 += s1[i - i1 - 1];
		}
		s1Aligned += missingGaps1;
		s1Aligned += missingText1;
		s2Aligned += missingGaps2;
		s2Aligned += missingText2;
	} else { // backtracking continues
		int sub = (ambiguousMatch(s1[i - 1], s2[j - 1])) ? 0 : MISMATCH_PENALTY;
		int gap = GAP_PENALTY;
		int vertical = matrix.entryAt(i - 1, j) + gap;
		int horizontal = matrix.entryAt(i, j - 1) + gap;
		int diagonal = matrix.entryAt(i - 1, j - 1) + sub;
		int min = matrix.entryAt(i, j);
		if (diagonal == min) {
			s1Aligned += s1[i - 1];
			s2Aligned += s2[j - 1];
			backtrack(i - 1, j - 1, s1Aligned, s2Aligned);
		} else if (horizontal == min) {
			s1Aligned += "-";
			s2Aligned += s2[j - 1];
			backtrack(i, j - 1, s1Aligned, s2Aligned);
		} else { // vertical == min
			s1Aligned += s1[i - 1];
			s2Aligned += "-";
			backtrack(i - 1, j, s1Aligned, s2Aligned);
		}
	}
}
void PairwiseAlignment::update(size_t i, size_t j) {
	if (i == 0) {
		matrix.entryAt(i, j) = matrix.entryAt(i, j - 1) + GAP_PENALTY;
	} else if (j == 0) {
		matrix.entryAt(i, j) = matrix.entryAt(i - 1, j) + GAP_PENALTY;
	} else {
		int sub = (ambiguousMatch(s1[i - 1], s2[j - 1])) ? 0 : MISMATCH_PENALTY;
		int gap = GAP_PENALTY;
		int vertical = matrix.entryAt(i - 1, j) + gap;
		int horizontal = matrix.entryAt(i, j - 1) + gap;
		int diagonal = matrix.entryAt(i - 1, j - 1) + sub;
		int min = std::min(vertical, std::min(horizontal, diagonal));
		matrix.entryAt(i, j) = min;
	}
}
