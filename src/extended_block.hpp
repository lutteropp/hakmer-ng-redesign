/*
 * extended_block.hpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <cmath>

#include "alignment/msa_wrapper.hpp"
#include "seeded_block.hpp"
#include "options.hpp"

inline double jukesCantorCorrection(double dist) {
	if (dist >= 0.75) { // Jukes Cantor Correction doesn't work if dist >= 0.75. In this case, it would return infinity. We change it to 1 here.
		return 1;
	}
	return -0.75 * std::log(1 - (4.0 / 3) * dist);
}

class ExtendedBlock {
public:
	ExtendedBlock(const SeededBlock& seededBlock, size_t nTax, bool noGaps);
	size_t getLeftFlankSize(size_t taxID) const;
	size_t getRightFlankSize(size_t taxID) const;
	void setLeftFlankSize(size_t taxID, size_t val);
	void setRightFlankSize(size_t taxID, size_t val);
	void incrementAllLeftFlanks();
	void incrementAllRightFlanks();
	void decrementAllLeftFlanks();
	void decrementAllRightFlanks();
	std::pair<size_t, size_t> getTaxonCoordsWithFlanks(size_t taxID) const;
	std::pair<size_t, size_t> getTaxonCoordsWithoutFlanks(size_t taxID) const;
	bool hasTaxon(size_t taxID) const;
	size_t getSeedSize() const;
	size_t getMaxLeftFlankSize() const;
	size_t getMaxRightFlankSize() const;
	size_t getNTaxInBlock() const;
	std::vector<size_t> getTaxonIDsInBlock() const;
	double getPairwiseNormalizedDistance(size_t idxInBlock1, size_t idxInBlock2, const Options& options);
	MSAWrapper msaWrapper;
private:
	SeededBlock mySeededBlock;
	std::vector<size_t> leftFlankSizes;
	std::vector<size_t> rightFlankSizes;
};
