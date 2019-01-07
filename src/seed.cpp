/*
 * seeded_block.cpp
 *
 *  Created on: Sep 4, 2018
 *      Author: Sarah Lutteropp
 */

#include "seed.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <utility>

Seed::Seed(size_t nTax) {
	seedCoords.resize(nTax);
	k = 0;
	originalK = 0;
	firstSAPos = 0;
	saPositions.resize(nTax);
	subRate = 0;
}

void Seed::addTaxon(size_t saPos, size_t taxID, size_t firstCoord, size_t lastCoord) {
	for (size_t i = 0; i < blockedTaxa.size(); ++i) {
		if (blockedTaxa[i] == taxID) {
			return;
		}
	}

	if (taxID >= seedCoords.size()) {
		throw std::runtime_error("Trying to add taxon ID that belongs to no taxon");
	}
	if (this->hasTaxon(taxID)) {
		std::cout << "taxID: " << taxID << "\n";
		std::cout << "firstCoord: " << firstCoord << "\n";
		std::cout << "second: " << lastCoord << "\n";
		std::cout << "old seedCoords[taxID].first: " << seedCoords[taxID].first << "\n";
		std::cout << "old seedCoords[taxID].second: " << seedCoords[taxID].second << "\n";
		throw std::runtime_error("This taxon is already present in the block");
	}
	seedCoords[taxID].first = firstCoord;
	seedCoords[taxID].second = lastCoord;
	k = lastCoord - firstCoord + 1;
	taxIDs.push_back(taxID);

	seedCoords[taxID].first = firstCoord;
	seedCoords[taxID].second = lastCoord;
	originalK = k;

	saPositions[taxID] = saPos;
}

size_t Seed::getNTaxInBlock() const {
	return taxIDs.size();
}

SimpleCoords Seed::getSeedCoords(size_t taxID) const {
	return seedCoords[taxID];
}

bool Seed::hasTaxon(size_t taxID) const {
	return (seedCoords[taxID].first != std::numeric_limits<size_t>::max()) && (seedCoords[taxID].second >= seedCoords[taxID].first);
}

void Seed::removeTaxon(size_t taxID) {
	seedCoords[taxID].first = std::numeric_limits<size_t>::max();
	seedCoords[taxID].second = std::numeric_limits<size_t>::max();
	seedCoords[taxID].leftGapSize = 0;
	seedCoords[taxID].rightGapSize = 0;
	taxIDs.erase(std::remove(taxIDs.begin(), taxIDs.end(), taxID), taxIDs.end());
}

double Seed::getAverageSeedSize() const {
	size_t seedSizeSum = 0;
	for (size_t tID : taxIDs) {
		if (!hasTaxon(tID)) {
			std::cout << "taxID: " << tID << "\n";
			std::cout << "old seedCoords[taxID].first: " << seedCoords[tID].first << "\n";
			std::cout << "old seedCoords[taxID].second: " << seedCoords[tID].second << "\n";
			throw std::runtime_error("This shouldn't happen");
		}
		seedSizeSum += getSeedSize(tID);
	}
	return (double) seedSizeSum / taxIDs.size();
}

std::vector<size_t> Seed::getTaxonIDsInBlock() const {
	return taxIDs;
}

void Seed::cleanTaxIDs() {
	std::vector<size_t> idsToRemove;
	for (size_t i = 0; i < taxIDs.size(); ++i) {
		if (!hasTaxon(taxIDs[i])) {
			idsToRemove.push_back(taxIDs[i]);
		}
	}
	for (size_t tID : idsToRemove) {
		removeTaxon(tID);
	}
}

void Seed::increaseAllTaxonCoordsRight() {
	for (size_t i = 0; i < taxIDs.size(); ++i) {
		seedCoords[taxIDs[i]].second++;
	}
	k++;
}

void Seed::decreaseAllTaxonCoordsLeft() {
	for (size_t i = 0; i < taxIDs.size(); ++i) {
		assert(seedCoords[taxIDs[i]].first > 0);
		seedCoords[taxIDs[i]].first--;
	}
	k++;
}

std::pair<size_t, size_t> computeForwardStrandCoordinates(const SimpleCoords& coords, size_t revCompStartPos) {
	size_t thisForwardFirst = coords.first;
	size_t thisForwardSecond = coords.second;
	if (thisForwardFirst >= revCompStartPos) {
		thisForwardFirst = revCompStartPos - (thisForwardFirst - revCompStartPos);
		thisForwardSecond = revCompStartPos - (thisForwardSecond - revCompStartPos);
	}
	return std::make_pair(std::min(thisForwardFirst, thisForwardSecond), std::max(thisForwardFirst, thisForwardSecond));
}

bool Seed::orderCompatible(const Seed& other, size_t revCompStartPos) const {
	bool orderSet = false;
	bool thisIsFirst;
	for (size_t i = 0; i < seedCoords.size(); ++i) {
		if (this->hasTaxon(i) && other.hasTaxon(i)) {
			// check if still order compatible
			bool thisIsFirstNow;
			size_t thisForwardFirst = seedCoords[i].first;
			if (thisForwardFirst >= revCompStartPos) {
				thisForwardFirst = revCompStartPos - (thisForwardFirst - revCompStartPos);
			}
			size_t otherForwardFirst = other.seedCoords[i].first;
			if (otherForwardFirst >= revCompStartPos) {
				otherForwardFirst = revCompStartPos - (otherForwardFirst - revCompStartPos);
			}
			thisIsFirstNow = (thisForwardFirst < otherForwardFirst);
			if (!orderSet) {
				thisIsFirst = thisIsFirstNow;
				orderSet = true;
			} else if (thisIsFirstNow != thisIsFirst) {
				return false;
			}
		}
	}
	return true;
}

bool Seed::overlap(const Seed& other, size_t revCompStartPos) const {
	for (size_t i = 0; i < seedCoords.size(); ++i) {
		if (this->hasTaxon(i) && other.hasTaxon(i)) {
			size_t thisForwardFirst = seedCoords[i].first;
			size_t thisForwardSecond = seedCoords[i].second;
			if (thisForwardFirst >= revCompStartPos) {
				thisForwardFirst = revCompStartPos - (thisForwardFirst - revCompStartPos);
				thisForwardSecond = revCompStartPos - (thisForwardSecond - revCompStartPos);
			}
			size_t otherForwardFirst = other.seedCoords[i].first;
			size_t otherForwardSecond = other.seedCoords[i].second;
			if (otherForwardFirst >= revCompStartPos) {
				otherForwardFirst = revCompStartPos - (otherForwardFirst - revCompStartPos);
				otherForwardSecond = revCompStartPos - (otherForwardSecond - revCompStartPos);
			}

			// check for overlap, faster version as in this blogpost: https://fgiesen.wordpress.com/2011/10/16/checking-for-interval-overlap/
			if (thisForwardFirst <= otherForwardSecond && otherForwardFirst <= thisForwardSecond) {
				return true;
			}
		}
	}
	return false;
}

size_t Seed::distance(const Seed& other, size_t revCompStartPos) const {
	if (overlap(other, revCompStartPos) || !orderCompatible(other, revCompStartPos)) {
		return std::numeric_limits<size_t>::infinity();
	}
	size_t dist = 0;
	for (size_t i = 0; i < seedCoords.size(); ++i) {
		if (this->hasTaxon(i) && other.hasTaxon(i)) {
			std::pair<size_t, size_t> thisForwardSorted = computeForwardStrandCoordinates(seedCoords[i], revCompStartPos);
			std::pair<size_t, size_t> otherForwardSorted = computeForwardStrandCoordinates(other.seedCoords[i], revCompStartPos);
			size_t actDist;
			if (thisForwardSorted.first > otherForwardSorted.second) {
				actDist = thisForwardSorted.first - otherForwardSorted.second;
			} else if (otherForwardSorted.first > thisForwardSorted.second) {
				actDist = otherForwardSorted.first - thisForwardSorted.second;
			} else {
				throw std::runtime_error("THIS SHOULD REALLY NOT HAPPEN");
			}
			dist = std::max(dist, actDist);
		}
	}
	return dist;
}

bool Seed::allSeedsSameSize() const {
	size_t size = 0;
	for (size_t i = 0; i < seedCoords.size(); ++i) {
		if (hasTaxon(i)) {
			size_t actSize = seedCoords[i].second + 1 - seedCoords[i].first;
			if (size != 0) {
				if (actSize != size) {
					return false;
				}
			} else {
				size = actSize;
			}
		}
	}
	return true;
}

size_t Seed::getSeedSize(size_t taxID) const {
	return seedCoords[taxID].size();
}

void Seed::increaseTaxonCoordsRight(size_t taxID) {
	seedCoords[taxID].second++;
}
void Seed::decreaseTaxonCoordsLeft(size_t taxID) {
	assert(seedCoords[taxID].first > 0);
	seedCoords[taxID].first--;
}
void Seed::decreaseTaxonCoordRight(size_t taxID) {
	seedCoords[taxID].second--;
}
void Seed::increaseTaxonCoordLeft(size_t taxID) {
	seedCoords[taxID].first++;
}

void Seed::addGapLeft(size_t taxonID) {
	seedCoords[taxonID].leftGapSize++;
}
void Seed::addGapRight(size_t taxonID) {
	seedCoords[taxonID].rightGapSize++;
}
void Seed::removeGapLeft(size_t taxonID) {
	seedCoords[taxonID].leftGapSize--;
}
void Seed::removeGapRight(size_t taxonID) {
	seedCoords[taxonID].rightGapSize--;
}

std::vector<SimpleCoords> Seed::getSeedCoords() const {
	return seedCoords;
}

size_t Seed::getOriginalK() const {
	return originalK;
}

size_t Seed::getFirstSAPos() const {
	return firstSAPos;
}

void Seed::setFirstSAPos(size_t pos) {
	firstSAPos = pos;
}

std::vector<size_t> Seed::getSAPositions() const {
	return saPositions;
}

void Seed::setSubRate(double rate) {
	subRate = rate;
}
double Seed::getSubRate() const {
	return subRate;
}

void Seed::blockTaxon(size_t taxID) {
	blockedTaxa.push_back(taxID);
}

const std::vector<size_t> Seed::getBlockedTaxa() const {
	return blockedTaxa;
}
