/*
 * extended_block.cpp
 *
 *  Created on: Sep 4, 2018
 *      Author: Sarah Lutteropp
 */

#include "extended_block.hpp"

ExtendedBlock::ExtendedBlock(const SeededBlock& seededBlock, size_t nTax, bool noGaps) :
		mySeededBlock(seededBlock), msaWrapper(noGaps) {
	leftFlankSize = 0;
	rightFlankSize = 0;
}

size_t ExtendedBlock::getLeftFlankSize() const {
	return leftFlankSize;
}

size_t ExtendedBlock::getRightFlankSize() const {
	return rightFlankSize;
}

void ExtendedBlock::setLeftFlankSize(size_t val) {
	leftFlankSize = val;
}

void ExtendedBlock::setRightFlankSize(size_t val) {
	rightFlankSize = val;
}

void ExtendedBlock::incrementLeftFlank() {
	leftFlankSize++;
}

void ExtendedBlock::incrementRightFlank() {
	rightFlankSize++;
}

void ExtendedBlock::decrementLeftFlank() {
	leftFlankSize--;
}

void ExtendedBlock::decrementRightFlank() {
	rightFlankSize--;
}

std::pair<size_t, size_t> ExtendedBlock::getTaxonCoordsWithFlanks(size_t taxID) const {
	return std::make_pair(mySeededBlock.getTaxonCoords(taxID).first - leftFlankSize,
			mySeededBlock.getTaxonCoords(taxID).second + rightFlankSize);
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
