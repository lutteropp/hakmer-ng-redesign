/*
 * extended_block.cpp
 *
 *  Created on: Sep 4, 2018
 *      Author: Sarah Lutteropp
 */

#include "extended_block.hpp"

ExtendedBlock::ExtendedBlock(const Seed& seededBlock, size_t nTax) :
		mySeededBlock(seededBlock) {
	leftFlankCoords.resize(nTax);
	rightFlankCoords.resize(nTax);
	for (size_t i = 0; i < nTax; ++i) {
		if (this->hasTaxon(i)) {
			leftFlankCoords[i].first = seededBlock.getSeedCoords(i).first;
			leftFlankCoords[i].second = seededBlock.getSeedCoords(i).first - 1;
			rightFlankCoords[i].first = seededBlock.getSeedCoords(i).second + 1;
			rightFlankCoords[i].second = seededBlock.getSeedCoords(i).second;
		}
	}
}

double ExtendedBlock::getAverageLeftFlankSize() const {
	size_t leftFlankSum = 0;
	for (size_t i = 0; i < leftFlankCoords.size(); ++i) {
		if (this->hasTaxon(i)) {
			leftFlankSum += leftFlankCoords[i].size();
		}
	}
	return (double) leftFlankSum / this->getNTaxInBlock();
}

double ExtendedBlock::getAverageRightFlankSize() const {
	size_t rightFlankSum = 0;
	for (size_t i = 0; i < rightFlankCoords.size(); ++i) {
		if (this->hasTaxon(i)) {
			rightFlankSum += rightFlankCoords[i].size();
		}
	}
	return (double) rightFlankSum / this->getNTaxInBlock();
}

void ExtendedBlock::setLeftFlankSize(size_t val) {
	for (size_t i = 0; i < leftFlankCoords.size(); ++i) {
		if (this->hasTaxon(i)) {
			// leftFlankCoords[i].second + 1 - leftFlankCoords[i].first === val
			// --> leftFlankCoords[i].second + 1 === val + leftFlankCoords[i].first
			// --> leftFlankCoords[i].second + 1 - val === leftFlankCoords[i].first
			leftFlankCoords[i].first = leftFlankCoords[i].second + 1 - val;
		}
	}
}

void ExtendedBlock::setRightFlankSize(size_t val) {
	for (size_t i = 0; i < rightFlankCoords.size(); ++i) {
		if (this->hasTaxon(i)) {
			// rightFlankCoords[i].second + 1 - rightFlankCoords[i].first === val
			// --> rightFlankCoords[i].second === val + rightFlankCoords[i].first - 1
			rightFlankCoords[i].second = val + rightFlankCoords[i].first - 1;
		}
	}
}

void ExtendedBlock::growLeftFlank(size_t tID) {
	leftFlankCoords[tID].first--;
}
void ExtendedBlock::growRightFlank(size_t tID) {
	rightFlankCoords[tID].second++;
}

std::pair<size_t, size_t> ExtendedBlock::getTaxonCoordsWithFlanks(size_t taxID) const {
	return std::make_pair(mySeededBlock.getSeedCoords(taxID).first - leftFlankCoords[taxID].size(),
			mySeededBlock.getSeedCoords(taxID).second + rightFlankCoords[taxID].size());
}

std::pair<size_t, size_t> ExtendedBlock::getTaxonCoordsWithoutFlanks(size_t taxID) const {
	return std::make_pair(mySeededBlock.getSeedCoords(taxID).first, mySeededBlock.getSeedCoords(taxID).second);
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
		sum += mySeededBlock.getSeedSize(tID) + leftFlankCoords[tID].size() + rightFlankCoords[tID].size();
	}
	return sum;
}

void ExtendedBlock::addLeftFlankGap(size_t tID) {
	leftFlankCoords[tID].leftGapSize++;
}
void ExtendedBlock::addRightFlankGap(size_t tID) {
	rightFlankCoords[tID].rightGapSize++;
}

std::vector<SimpleCoords> ExtendedBlock::getLeftFlankCoords() const {
	return leftFlankCoords;
}
std::vector<SimpleCoords> ExtendedBlock::getSeedCoords() const {
	return mySeededBlock.getSeedCoords();
}
std::vector<SimpleCoords> ExtendedBlock::getRightFlankCoords() const {
	return rightFlankCoords;
}
