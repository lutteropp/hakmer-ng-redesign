/*
 * nogaps_msa.cpp
 *
 *  Created on: Oct 9, 2018
 *      Author: Sarah Lutteropp
 */

#include "nogaps_msa.hpp"

#include <algorithm>

NoGapsMSA::NoGapsMSA() :
		nTax(0), width(0) {
}

void NoGapsMSA::init(size_t nTax) {
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
std::vector<std::string> NoGapsMSA::assembleMSA() {
	std::vector<std::string> msa(nTax, "");
	for (size_t i = 0; i < nTax; ++i) {
		std::string reversedLeft = sequencesLeft[i];
		std::reverse(reversedLeft.begin(), reversedLeft.end());
		msa[i] = reversedLeft + sequencesMiddle[i] + sequencesRight[i];
	}
	return msa;
}
double NoGapsMSA::pairwiseDistance(size_t idx1, size_t idx2) {
	size_t firstIdx = std::min(idx1, idx2);
	size_t secondIdx = std::max(idx1, idx2);
	return pairwiseHammingDistances.entryAt(firstIdx, secondIdx);
}

double NoGapsMSA::normalizedPairwiseDistance(size_t idx1, size_t idx2) {
	size_t firstIdx = std::min(idx1, idx2);
	size_t secondIdx = std::max(idx1, idx2);
	return (double) pairwiseHammingDistances.entryAt(firstIdx, secondIdx) / width;
}

void NoGapsMSA::addCharsLeft(const std::vector<char>& chars) {
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
void NoGapsMSA::addCharsRight(const std::vector<char>& chars) {
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
void NoGapsMSA::setSeeds(const std::vector<std::string>& seeds) {
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
void NoGapsMSA::setSeeds(const std::string& seed) {
	for (size_t i = 0; i < nTax; ++i) {
		sequencesMiddle[i] = seed;
	}
	width += seed.size();
}

void NoGapsMSA::shrinkDownToLeftFlank(size_t newLeftFlankSize) {
	for (size_t i = 0; i < nTax - 1; ++i) {
		for (size_t j = i + 1; j < nTax; ++j) {
			for (size_t k = newLeftFlankSize - 1; k < sequencesLeft[0].size(); ++k) {
				if (!isGapCharacter(sequencesLeft[i][k]) && !isGapCharacter(sequencesLeft[j][k])
						&& sequencesLeft[i][k] != sequencesLeft[j][k]) {
					pairwiseHammingDistances.entryAt(i, j) = pairwiseHammingDistances.entryAt(i, j) - 1; // because this position will be deleted
				}
			}
		}
	}
	width -= sequencesLeft[0].size() - newLeftFlankSize;
	for (size_t i = 0; i < nTax; ++i) {
		sequencesLeft[i].resize(newLeftFlankSize);
		sequencesLeft[i].shrink_to_fit();
	}
}

void NoGapsMSA::shrinkDownToRightFlank(size_t newRightFlankSize) {
	for (size_t i = 0; i < nTax - 1; ++i) {
		for (size_t j = i + 1; j < nTax; ++j) {
			for (size_t k = newRightFlankSize - 1; k < sequencesRight[0].size(); ++k) {
				if (!isGapCharacter(sequencesRight[i][k]) && !isGapCharacter(sequencesRight[j][k])
						&& sequencesRight[i][k] != sequencesRight[j][k]) {
					pairwiseHammingDistances.entryAt(i, j) = pairwiseHammingDistances.entryAt(i, j) - 1; // because this position will be deleted
				}
			}
		}
	}
	width -= sequencesRight[0].size() - newRightFlankSize;
	for (size_t i = 0; i < nTax; ++i) {
		sequencesRight[i].resize(newRightFlankSize);
		sequencesRight[i].shrink_to_fit();
	}
}

bool NoGapsMSA::isGapCharacter(char c) {
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
