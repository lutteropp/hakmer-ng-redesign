/*
 * extended_block.cpp
 *
 *  Created on: Sep 4, 2018
 *      Author: Sarah Lutteropp
 */

#include "extended_block.hpp"

ExtendedBlock::ExtendedBlock(const Seed& seededBlock, size_t nTax, bool noGaps) :
		mySeededBlock(seededBlock), msaWrapper(noGaps) {
	leftFlankSizes.resize(nTax);
	rightFlankSizes.resize(nTax);
}

double ExtendedBlock::getAverageLeftFlankSize() const {
	size_t leftFlankSum = 0;
	for (size_t i = 0; i < leftFlankSizes.size(); ++i) {
		if (this->hasTaxon(i)) {
			leftFlankSum += leftFlankSizes[i];
		}
	}
	return (double) leftFlankSum / this->getNTaxInBlock();
}

double ExtendedBlock::getAverageRightFlankSize() const {
	size_t rightFlankSum = 0;
	for (size_t i = 0; i < rightFlankSizes.size(); ++i) {
		if (this->hasTaxon(i)) {
			rightFlankSum += rightFlankSizes[i];
		}
	}
	return (double) rightFlankSum / this->getNTaxInBlock();
}

void ExtendedBlock::setLeftFlankSize(size_t val) {
	for (size_t i = 0; i < leftFlankSizes.size(); ++i) {
		if (this->hasTaxon(i)) {
			leftFlankSizes[i] = val;
		}
	}
}

void ExtendedBlock::setRightFlankSize(size_t val) {
	for (size_t i = 0; i < rightFlankSizes.size(); ++i) {
		if (this->hasTaxon(i)) {
			rightFlankSizes[i] = val;
		}
	}
}

void ExtendedBlock::decrementLeftFlank(size_t tID) {
	leftFlankSizes[tID]--;
}
void ExtendedBlock::incrementRightFlank(size_t tID) {
	rightFlankSizes[tID]++;
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

double ExtendedBlock::getAverageSeedSize() const {
	return mySeededBlock.getAverageSeedSize();
}

size_t ExtendedBlock::getNTaxInBlock() const {
	return mySeededBlock.getNTaxInBlock();
}

std::vector<size_t> ExtendedBlock::getTaxonIDsInBlock() const {
	return mySeededBlock.getTaxonIDsInBlock();
}

size_t ExtendedBlock::getTotalBasesUsed() const {
	size_t sum = 0;
	for (size_t tID : mySeededBlock.getTaxonIDsInBlock()) {
		sum += mySeededBlock.getSeedSize(tID) + leftFlankSizes[tID] + rightFlankSizes[tID];
	}
	return sum;
}
