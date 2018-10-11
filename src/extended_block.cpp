/*
 * extended_block.cpp
 *
 *  Created on: Sep 4, 2018
 *      Author: Sarah Lutteropp
 */

#include "extended_block.hpp"

ExtendedBlock::ExtendedBlock(const SeededBlock& seededBlock, size_t nTax, bool noGaps) :
		mySeededBlock(seededBlock), msaWrapper(noGaps) {
	leftFlankSizes.resize(nTax);
	rightFlankSizes.resize(nTax);
	for (size_t i = 0; i < nTax; ++i) {
		leftFlankSizes[i] = 0;
		rightFlankSizes[i] = 0;
	}
}

size_t ExtendedBlock::getLeftFlankSize(size_t taxID) const {
	return leftFlankSizes[taxID];
}

size_t ExtendedBlock::getRightFlankSize(size_t taxID) const {
	return rightFlankSizes[taxID];
}

void ExtendedBlock::setLeftFlankSize(size_t taxID, size_t val) {
	leftFlankSizes[taxID] = val;
}

void ExtendedBlock::setRightFlankSize(size_t taxID, size_t val) {
	rightFlankSizes[taxID] = val;
}

void ExtendedBlock::incrementAllLeftFlanks() {
	for (size_t i = 0; i < leftFlankSizes.size(); ++i) {
		if (mySeededBlock.hasTaxon(i)) {
			leftFlankSizes[i]++;
		}
	}
}

void ExtendedBlock::incrementAllRightFlanks() {
	for (size_t i = 0; i < rightFlankSizes.size(); ++i) {
		if (mySeededBlock.hasTaxon(i)) {
			rightFlankSizes[i]++;
		}
	}
}

void ExtendedBlock::decrementAllLeftFlanks() {
	for (size_t i = 0; i < leftFlankSizes.size(); ++i) {
		if (mySeededBlock.hasTaxon(i)) {
			leftFlankSizes[i]--;
		}
	}
}

void ExtendedBlock::decrementAllRightFlanks() {
	for (size_t i = 0; i < rightFlankSizes.size(); ++i) {
		if (mySeededBlock.hasTaxon(i)) {
			rightFlankSizes[i]--;
		}
	}
}

std::pair<size_t, size_t> ExtendedBlock::getTaxonCoordsWithFlanks(size_t taxID) const {
	return std::make_pair(mySeededBlock.getTaxonCoords(taxID).first - leftFlankSizes[taxID],
			mySeededBlock.getTaxonCoords(taxID).second + rightFlankSizes[taxID]);
}

std::pair<size_t, size_t> ExtendedBlock::getTaxonCoordsWithoutFlanks(size_t taxID) const {
	return mySeededBlock.getTaxonCoords(taxID);
}

bool ExtendedBlock::hasTaxon(size_t taxID) const {
	return mySeededBlock.hasTaxon(taxID);
}

size_t ExtendedBlock::getSeedSize() const {
	return mySeededBlock.getSeedSize();
}

size_t ExtendedBlock::getMaxLeftFlankSize() const {
	size_t max = 0;
	for (size_t i = 0; i < leftFlankSizes.size(); ++i) {
		if (leftFlankSizes[i] > max) {
			max = leftFlankSizes[i];
		}
	}
	return max;
}

size_t ExtendedBlock::getMaxRightFlankSize() const {
	size_t max = 0;
	for (size_t i = 0; i < rightFlankSizes.size(); ++i) {
		if (rightFlankSizes[i] > max) {
			max = rightFlankSizes[i];
		}
	}
	return max;
}

size_t ExtendedBlock::getNTaxInBlock() const {
	return mySeededBlock.getNTaxInBlock();
}

std::vector<size_t> ExtendedBlock::getTaxonIDsInBlock() const {
	return mySeededBlock.getTaxonIDsInBlock();
}

double ExtendedBlock::getPairwiseNormalizedDistance(size_t idxInBlock1, size_t idxInBlock2, const Options& options) {
	double dist = msaWrapper.normalizedPairwiseDistance(idxInBlock1, idxInBlock2);
	if (options.jukesCantor) {
		dist = jukesCantorCorrection(dist);
	}
	return dist;
}

std::vector<double> ExtendedBlock::getPairwiseNormalizedDistances(const Options& options) {
	std::vector<double> res;
	for (size_t i = 0; i < getNTaxInBlock(); ++i) {
		for (size_t j = i + 1; j < getNTaxInBlock(); ++j) {
			res.push_back(getPairwiseNormalizedDistance(i, j, options));
		}
	}
	return res;
}
