/*
 * star_alignment.hpp
 *
 *  Created on: Oct 5, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <stdexcept>

#include "options.hpp"

// 1.) Build all O(n^2) pairwise alignments.
// 2.) Choose center sequence S that is closest to all other sequences (this is, sum of edit distances to the other sequences is minimal).
// 3.) Incrementally merge the pairwise alignments into one MSA (Progressive alignment). Gaps in S from the pairwise alignment stay, insert gaps to the other sequences as needed.

// Additional difficulty: We need to distinguish between left added flanks, seeds, and right added flanks.
// Still an open question: How to select the alignment penalties (indels, substitutions)?
// Open thing to do: Improve runtime and space-requirements of the pairwise alignments...
// Maybe we could use edlib for very performant pairwise alignment (but it only supports Levenshtein distance)?

// TODO: Maybe also support affine gap penalties?

const int GAP_PENALTY = 2;
const int MISMATCH_PENALTY = 1;

template<typename T>
class TwoDimMatrix {
public:
	TwoDimMatrix() :
			n(0), m(0) {
	}
	void init(size_t n, size_t m) {
		this->n = n;
		this->m = m;
		entries.resize(n);
		for (size_t i = 0; i < m; ++i) {
			entries[i].resize(m);
		}
	}
	void addRow() {
		n++;
		entries.resize(n);
		entries[n - 1].resize(m);
	}
	void addColumn() {
		m++;
		for (size_t i = 0; i < entries.size(); ++i) {
			entries[i].resize(m);
		}
	}
	T& entryAt(size_t i, size_t j) {
		return entries[i][j];
	}
	size_t getN() const {
		return n;
	}
	size_t getM() const {
		return m;
	}
private:
	size_t n;
	size_t m;
	std::vector<std::vector<T> > entries;
};

class PairwiseAlignment {
public:
	PairwiseAlignment() :
			s1(""), s2(""), aliValid(false) {
		matrix.addRow();
		matrix.addColumn();
		// now we have a 1x1 matrix
		matrix.entryAt(0, 0) = 0;
	}
	std::string getS1() const {
		return s1;
	}
	std::string getS2() const {
		return s2;
	}
	std::pair<std::string, std::string> extractAlignment() {
		if (aliValid) {
			return ali;
		}
		std::string s1Aligned;
		std::string s2Aligned;
		backtrack(s1.size(), s2.size(), s1Aligned, s2Aligned);
		ali = std::make_pair(s1Aligned, s2Aligned);
		aliValid = true;
		return ali;
	}
	void addChars(char a, char b) {
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
	double pairwiseDistance() {
		return matrix.entryAt(s1.size(), s2.size());
	}
	size_t getAlignmentWidth() {
		std::pair<std::string, std::string> ali = extractAlignment();
		return ali.first.size();
	}
private:
	void backtrack(size_t i, size_t j, std::string& s1Aligned, std::string& s2Aligned) {
		// TODO: This can be optimized by assembling the reverse aligned sequences and then later on reverse them
		if (i == 0 || j == 0) { // backtracking has ended
		// add missing gaps
			std::string missingGaps1 = "";
			std::string missingText2 = "";
			for (size_t j1 = 0; j1 < j; ++j1) {
				missingGaps1 += "-";
				missingText2 = s2[j - j1 - 1] + missingText2;
			}
			std::string missingGaps2 = "";
			std::string missingText1 = "";
			for (size_t i1 = 0; i1 < i; ++i1) {
				missingGaps2 += "-";
				missingText1 = s1[i - i1 - 1] + missingText1;
			}
			s1Aligned = missingText1 + missingGaps1 + s1Aligned;
			s2Aligned = missingText2 + missingGaps2 + s2Aligned;
		} else { // backtracking continues
			int sub = (s1[i - 1] == s2[j - 1]) ? 0 : MISMATCH_PENALTY;
			int gap = GAP_PENALTY;
			int vertical = matrix.entryAt(i - 1, j) + gap;
			int horizontal = matrix.entryAt(i, j - 1) + gap;
			int diagonal = matrix.entryAt(i - 1, j - 1) + sub;
			int min = matrix.entryAt(i, j);
			if (diagonal == min) {
				s1Aligned = s1[i - 1] + s1Aligned;
				s2Aligned = s2[j - 1] + s2Aligned;
				backtrack(i - 1, j - 1, s1Aligned, s2Aligned);
			} else if (horizontal == min) {
				s1Aligned = "-" + s1Aligned;
				s2Aligned = s2[j - 1] + s2Aligned;
				backtrack(i, j - 1, s1Aligned, s2Aligned);
			} else { // vertical == min
				s1Aligned = s1[i - 1] + s1Aligned;
				s2Aligned = "-" + s2Aligned;
				backtrack(i - 1, j, s1Aligned, s2Aligned);
			}
		}
	}
	void update(size_t i, size_t j) {
		if (i == 0) {
			matrix.entryAt(i, j) = matrix.entryAt(i, j - 1) + GAP_PENALTY;
		} else if (j == 0) {
			matrix.entryAt(i, j) = matrix.entryAt(i - 1, j) + GAP_PENALTY;
		} else {
			int sub = (s1[i - 1] == s2[j - 1]) ? 0 : MISMATCH_PENALTY;
			int gap = GAP_PENALTY;
			int vertical = matrix.entryAt(i - 1, j) + gap;
			int horizontal = matrix.entryAt(i, j - 1) + gap;
			int diagonal = matrix.entryAt(i - 1, j - 1) + sub;
			int min = std::min(vertical, std::min(horizontal, diagonal));
			matrix.entryAt(i, j) = min;
		}
	}
	std::string s1;
	std::string s2;
	TwoDimMatrix<int> matrix;
	std::pair<std::string, std::string> ali;
	bool aliValid;
};

class SplitPairwiseAlignment {
public:
	SplitPairwiseAlignment() {
	}
	void addCharsLeft(char a, char b) {
		leftFlank.addChars(a, b);
	}
	void addCharsRight(char a, char b) {
		rightFlank.addChars(a, b);
	}
	void setSeed(const std::string& seed1, const std::string& seed2) {
		if (seed1.size() != seed2.size()) {
			throw std::runtime_error("Different-sized seeds not supported yet");
		}
		for (size_t i = 0; i < seed1.size(); ++i) {
			seed.addChars(seed1[i], seed2[i]);
		}
	}
	std::pair<std::string, std::string> extractAlignment() {
		std::string s1Aligned;
		std::string s2Aligned;
		std::pair<std::string, std::string> leftFlankAlignment = leftFlank.extractAlignment();
		std::reverse(leftFlankAlignment.first.begin(), leftFlankAlignment.first.end());
		std::reverse(leftFlankAlignment.second.begin(), leftFlankAlignment.second.end());
		s1Aligned += leftFlankAlignment.first;
		s2Aligned += leftFlankAlignment.second;
		std::pair<std::string, std::string> seedAlignment = seed.extractAlignment();
		s1Aligned += seedAlignment.first;
		s2Aligned += seedAlignment.second;
		std::pair<std::string, std::string> rightFlankAlignment = rightFlank.extractAlignment();
		s1Aligned += rightFlankAlignment.first;
		s2Aligned += rightFlankAlignment.second;
		return std::make_pair(s1Aligned, s2Aligned);
	}
	size_t getAlignmentWidth() {
		return leftFlank.getAlignmentWidth() + seed.getAlignmentWidth() + rightFlank.getAlignmentWidth();
	}
	std::string getS1() const {
		std::string s1Left = leftFlank.getS1();
		std::reverse(s1Left.begin(), s1Left.end());
		return s1Left + seed.getS1() + rightFlank.getS1();
	}
	std::string getS2() const {
		std::string s2Left = leftFlank.getS2();
		std::reverse(s2Left.begin(), s2Left.end());
		return s2Left + seed.getS2() + rightFlank.getS2();
	}
	double pairwiseDistance() {
		return leftFlank.pairwiseDistance() + seed.pairwiseDistance() + rightFlank.pairwiseDistance();
	}
private:
	PairwiseAlignment leftFlank;
	PairwiseAlignment seed;
	PairwiseAlignment rightFlank;
};

class StarMSA {
public:
	StarMSA() :
			nTax(0) {
	}
	void init(size_t nTax) {
		this->nTax = nTax;
		pairwiseAlignments.init(nTax, nTax);
	}
	std::vector<std::string> assembleMSA() {
		std::vector<std::string> msa;
		msa.resize(nTax);
		// first, find the sequence with the smallest distance sum to the others
		std::vector<int> distanceSums(nTax, 0);
		for (size_t i = 0; i < nTax; ++i) {
			for (size_t j = 0; j < nTax; ++j) {
				if (i == j)
					continue;
				size_t firstIdx = std::min(i, j);
				size_t secondIdx = std::max(i, j);
				distanceSums[i] += pairwiseAlignments.entryAt(firstIdx, secondIdx).pairwiseDistance();
			}
		}
		size_t smallestIdx = 0;
		int smallestDist = distanceSums[0];
		for (size_t i = 1; i < nTax; ++i) {
			if (distanceSums[i] < smallestDist) {
				smallestIdx = i;
				smallestDist = distanceSums[i];
			}
		}
		for (size_t i = 0; i < nTax; ++i) {
			if (i == smallestIdx)
				continue;
			addToMSA(i, msa, smallestIdx);
		}
		return msa;
	}
	double pairwiseDistance(size_t idx1, size_t idx2) {
		size_t firstIdx = std::min(idx1, idx2);
		size_t secondIdx = std::max(idx1, idx2);
		return pairwiseAlignments.entryAt(firstIdx, secondIdx).pairwiseDistance();
	}
	double normalizedPairwiseDistance(size_t idx1, size_t idx2) {
		size_t firstIdx = std::min(idx1, idx2);
		size_t secondIdx = std::max(idx1, idx2);
		double editDist = pairwiseAlignments.entryAt(firstIdx, secondIdx).pairwiseDistance();
		size_t s1Size = pairwiseAlignments.entryAt(firstIdx, secondIdx).getS1().size();
		size_t s2Size = pairwiseAlignments.entryAt(firstIdx, secondIdx).getS2().size();
		double normEditDist = (2.0 * editDist) / (GAP_PENALTY * (s1Size + s2Size) + editDist);
		return normEditDist;
	}

	void addCharsLeft(const std::vector<char>& chars) {
		for (size_t i = 0; i < nTax - 1; ++i) {
			for (size_t j = i + 1; j < nTax; ++j) {
				pairwiseAlignments.entryAt(i, j).addCharsLeft(chars[i], chars[j]);
			}
		}
	}
	void addCharsRight(const std::vector<char>& chars) {
		for (size_t i = 0; i < nTax - 1; ++i) {
			for (size_t j = i + 1; j < nTax; ++j) {
				pairwiseAlignments.entryAt(i, j).addCharsRight(chars[i], chars[j]);
			}
		}
	}
	void setSeeds(const std::vector<std::string>& seeds) {
		for (size_t i = 0; i < nTax - 1; ++i) {
			for (size_t j = i + 1; j < nTax; ++j) {
				pairwiseAlignments.entryAt(i, j).setSeed(seeds[i], seeds[j]);
			}
		}
	}
	void setSeeds(const std::string& seed) {
		for (size_t i = 0; i < nTax - 1; ++i) {
			for (size_t j = i + 1; j < nTax; ++j) {
				pairwiseAlignments.entryAt(i, j).setSeed(seed, seed);
			}
		}
	}
private:
	void addToMSA(size_t taxonToAdd, std::vector<std::string>& msa, size_t centerSequenceIdx) {
		size_t firstIdx = std::min(taxonToAdd, centerSequenceIdx);
		size_t secondIdx = std::max(taxonToAdd, centerSequenceIdx);
		std::pair<std::string, std::string> alignmentToAdd = pairwiseAlignments.entryAt(firstIdx, secondIdx).extractAlignment();
		std::string aliSeqCenter;
		std::string aliSeqNewTaxon;
		if (centerSequenceIdx == firstIdx) {
			aliSeqCenter = alignmentToAdd.first;
			aliSeqNewTaxon = alignmentToAdd.second;
		} else {
			aliSeqCenter = alignmentToAdd.second;
			aliSeqNewTaxon = alignmentToAdd.first;
		}

		if (msa[centerSequenceIdx].empty()) {
			msa[centerSequenceIdx] += aliSeqCenter;
			msa[taxonToAdd] += aliSeqNewTaxon;
		} else {
			size_t i = 0;
			while (true) {
				if (msa[centerSequenceIdx][i] == aliSeqCenter[i]) {
					msa[taxonToAdd] += aliSeqNewTaxon[i];
				} else if (msa[centerSequenceIdx][i] == '-') {
					msa[taxonToAdd] += "-";
					msa[taxonToAdd] += aliSeqNewTaxon[i];
				} else if (aliSeqCenter[i] == '-') {
					// add gap to all sequences already added to the msa, at position i
					for (size_t j = 0; j < msa.size(); ++j) {
						if (msa[j].empty())
							continue;
						std::string left = msa[j].substr(0, i);
						std::string right = msa[j].substr(i, std::string::npos);
						msa[j] = left + "-" + right;
					}
					msa[taxonToAdd] += aliSeqNewTaxon[i];
				}
				++i;
			}
		}
	}

	size_t nTax;
	TwoDimMatrix<SplitPairwiseAlignment> pairwiseAlignments;
};

class NoGapsMSA {
public:
	NoGapsMSA() :
			nTax(0), width(0) {
	}

	void init(size_t nTax) {
		this->nTax = nTax;
		pairwiseHammingDistances.init(nTax, nTax);
		sequencesLeft.resize(nTax);
		sequencesMiddle.resize(nTax);
		sequencesRight.resize(nTax);
		width = 0;

		for (size_t i = 0; i < nTax; ++i) {
			for (size_t j = i + 1; j < nTax; ++j) {
				pairwiseHammingDistances.entryAt(i, j) = 0;
			}
		}
	}
	std::vector<std::string> assembleMSA() {
		std::vector<std::string> msa(nTax, "");
		for (size_t i = 0; i < nTax; ++i) {
			std::string reversedLeft = sequencesLeft[i];
			std::reverse(reversedLeft.begin(), reversedLeft.end());
			msa[i] = reversedLeft + sequencesMiddle[i] + sequencesRight[i];
		}
		return msa;
	}
	double pairwiseDistance(size_t idx1, size_t idx2) {
		size_t firstIdx = std::min(idx1, idx2);
		size_t secondIdx = std::max(idx1, idx2);
		return pairwiseHammingDistances.entryAt(firstIdx, secondIdx);
	}

	double normalizedPairwiseDistance(size_t idx1, size_t idx2) {
		size_t firstIdx = std::min(idx1, idx2);
		size_t secondIdx = std::max(idx1, idx2);
		return (double) pairwiseHammingDistances.entryAt(firstIdx, secondIdx) / width;
	}

	void addCharsLeft(const std::vector<char>& chars) {
		for (size_t i = 0; i < nTax; ++i) {
			sequencesLeft[i] += chars[i];
		}
		for (size_t i = 0; i < nTax; ++i) {
			for (size_t j = i + 1; j < nTax; ++j) {
				if (!isGapCharacter(chars[i]) && !isGapCharacter(chars[j]) && chars[i] != chars[j]) {
					pairwiseHammingDistances.entryAt(i, j) = pairwiseHammingDistances.entryAt(i, j) + 1;
				}
			}
		}
		width++;
	}
	void addCharsRight(const std::vector<char>& chars) {
		for (size_t i = 0; i < nTax; ++i) {
			sequencesRight[i] += chars[i];
		}
		for (size_t i = 0; i < nTax; ++i) {
			for (size_t j = i + 1; j < nTax; ++j) {
				if (!isGapCharacter(chars[i]) && !isGapCharacter(chars[j]) && chars[i] != chars[j]) {
					pairwiseHammingDistances.entryAt(i, j) = pairwiseHammingDistances.entryAt(i, j) + 1;
				}
			}
		}
		width++;
	}
	void setSeeds(const std::vector<std::string>& seeds) {
		for (size_t i = 0; i < nTax; ++i) {
			sequencesMiddle[i] = seeds[i];
		}

		for (size_t i = 0; i < nTax - 1; ++i) {
			for (size_t j = i + 1; j < nTax; ++j) {
				for (size_t k = 0; k < seeds[i].size(); ++k) {
					if ((!isGapCharacter(seeds[i][k])) && (!isGapCharacter(seeds[j][k])) && (seeds[i][k] != seeds[j][k])) {
						pairwiseHammingDistances.entryAt(i, j) = pairwiseHammingDistances.entryAt(i, j) + 1;
					}
				}
			}
		}
		width += seeds[0].size();
	}
	void setSeeds(const std::string& seed) {
		for (size_t i = 0; i < nTax; ++i) {
			sequencesMiddle[i] = seed;
		}
		width += seed.size();
	}
private:
	bool isGapCharacter(char c) {
		if (c == 'N' || c == 'n' || c == '-' || c == '?') {
			return true;
		} else {
			return false;
		}
	}
	size_t nTax;
	TwoDimMatrix<size_t> pairwiseHammingDistances;
	std::vector<std::string> sequencesLeft, sequencesMiddle, sequencesRight;
	size_t width;
};
