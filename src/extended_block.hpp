/*
 * extended_block.hpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include "seeded_block.hpp"

class ExtendedBlock {
public:
	ExtendedBlock(const SeededBlock& seededBlock, size_t nTax);
	size_t getLeftFlankSize(size_t taxID);
	size_t getRightFlankSize(size_t taxID);
	void setLeftFlankSize(size_t taxID, size_t val);
	void setRightFlankSize(size_t taxID, size_t val);
	void incrementAllLeftFlanks();
	void incrementAllRightFlanks();
	void decrementAllLeftFlanks();
	void decrementAllRightFlanks();
	std::pair<size_t, size_t> getTaxonCoordsWithFlanks(size_t taxID);
	bool hasTaxon(size_t taxID);
private:
	SeededBlock mySeededBlock;
	std::vector<size_t> leftFlankSizes;
	std::vector<size_t> rightFlankSizes;
};

ExtendedBlock::ExtendedBlock(const SeededBlock& seededBlock, size_t nTax) :
		mySeededBlock(seededBlock) {
	leftFlankSizes.resize(nTax);
	rightFlankSizes.resize(nTax);
	for (size_t i = 0; i < nTax; ++i) {
		leftFlankSizes[i] = 0;
		rightFlankSizes[i] = 0;
	}
}

size_t ExtendedBlock::getLeftFlankSize(size_t taxID) {
	return leftFlankSizes[taxID];
}

size_t ExtendedBlock::getRightFlankSize(size_t taxID) {
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

std::pair<size_t, size_t> ExtendedBlock::getTaxonCoordsWithFlanks(size_t taxID) {
	return std::make_pair(mySeededBlock.getTaxonCoords()[taxID].first - leftFlankSizes[taxID],
			mySeededBlock.getTaxonCoords()[taxID].second + rightFlankSizes[taxID]);
}

bool ExtendedBlock::hasTaxon(size_t taxID) {
	return mySeededBlock.hasTaxon(taxID);
}
