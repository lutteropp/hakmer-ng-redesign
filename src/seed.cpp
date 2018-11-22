/*
 * seeded_block.cpp
 *
 *  Created on: Sep 4, 2018
 *      Author: Sarah Lutteropp
 */

#include <stdexcept>
#include <iostream>
#include <limits>
#include <algorithm>
#include <cassert>
#include "seed.hpp"

Seed::Seed(size_t nTax) {
	taxonCoords.resize(nTax);
	for (size_t i = 0; i < nTax; ++i) {
		taxonCoords[i].first = std::string::npos;
		taxonCoords[i].second = std::string::npos;
	}
	n = 0;
	k = 0;
}

void Seed::addTaxon(size_t taxID, size_t firstCoord, size_t lastCoord) {
	if (taxID >= taxonCoords.size()) {
		throw std::runtime_error("Trying to add taxon ID that belongs to no taxon");
	}
	if (taxonCoords[taxID].first != std::string::npos) {
		std::cout << "taxID: " << taxID << "\n";
		std::cout << "firstCoord: " << firstCoord << "\n";
		std::cout << "lastCoord: " << lastCoord << "\n";
		std::cout << "old taxonCoords[taxID].first: " << taxonCoords[taxID].first << "\n";
		std::cout << "old taxonCoords[taxID].second: " << taxonCoords[taxID].second << "\n";
		throw std::runtime_error("This taxon is already present in the block");
	}
	n++;
	taxonCoords[taxID].first = firstCoord;
	taxonCoords[taxID].second = lastCoord;
	k = lastCoord - firstCoord + 1;
	taxIDs.push_back(taxID);
}

size_t Seed::getNTaxInBlock() const {
	return n;
}

std::pair<size_t, size_t> Seed::getTaxonCoords(size_t taxID) const {
	return taxonCoords[taxID];
}

bool Seed::hasTaxon(size_t taxID) const {
	return (taxonCoords[taxID].first != std::string::npos) && (taxonCoords[taxID].second > taxonCoords[taxID].first);
}

void Seed::removeTaxon(size_t taxID) {
	taxonCoords[taxID].first = std::string::npos;
	taxonCoords[taxID].second = std::string::npos;
	taxIDs.erase(std::remove(taxIDs.begin(), taxIDs.end(), taxID), taxIDs.end());
	n--;
}

double Seed::getAverageSeedSize() const {
	size_t seedSizeSum = 0;
	for (size_t tID : taxIDs) {
		seedSizeSum += getSeedSize(tID);
	}
	return (double) seedSizeSum / taxIDs.size();
}

std::vector<size_t> Seed::getTaxonIDsInBlock() const {
	return taxIDs;
}

void Seed::increaseAllTaxonCoordsRight() {
	for (size_t i = 0; i < taxIDs.size(); ++i) {
		taxonCoords[taxIDs[i]].second++;
	}
	k++;
}

void Seed::decreaseAllTaxonCoordsLeft() {
	for (size_t i = 0; i < taxIDs.size(); ++i) {
		assert(taxonCoords[taxIDs[i]].first > 0);
		taxonCoords[taxIDs[i]].first--;
	}
	k++;
}

void Seed::decreaseAllTaxonCoordsRight() {
	for (size_t i = 0; i < taxIDs.size(); ++i) {
		taxonCoords[taxIDs[i]].second--;
	}
	k--;
}
void Seed::increaseAllTaxonCoordsLeft() {
	for (size_t i = 0; i < taxIDs.size(); ++i) {
		taxonCoords[taxIDs[i]].first++;
	}
	k--;
}

std::pair<size_t, size_t> computeForwardStrandCoordinates(const std::pair<size_t, size_t>& coords, size_t revCompStartPos) {
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
	for (size_t i = 0; i < taxonCoords.size(); ++i) {
		if (this->hasTaxon(i) && other.hasTaxon(i)) {
			// check if still order compatible
			bool thisIsFirstNow;
			size_t thisForwardFirst = taxonCoords[i].first;
			if (thisForwardFirst >= revCompStartPos) {
				thisForwardFirst = revCompStartPos - (thisForwardFirst - revCompStartPos);
			}
			size_t otherForwardFirst = other.taxonCoords[i].first;
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
	for (size_t i = 0; i < taxonCoords.size(); ++i) {
		if (this->hasTaxon(i) && other.hasTaxon(i)) {
			size_t thisForwardFirst = taxonCoords[i].first;
			size_t thisForwardSecond = taxonCoords[i].second;
			if (thisForwardFirst >= revCompStartPos) {
				thisForwardFirst = revCompStartPos - (thisForwardFirst - revCompStartPos);
				thisForwardSecond = revCompStartPos - (thisForwardSecond - revCompStartPos);
			}
			size_t otherForwardFirst = other.taxonCoords[i].first;
			size_t otherForwardSecond = other.taxonCoords[i].second;
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
	for (size_t i = 0; i < taxonCoords.size(); ++i) {
		if (this->hasTaxon(i) && other.hasTaxon(i)) {
			std::pair<size_t, size_t> thisForwardSorted = computeForwardStrandCoordinates(taxonCoords[i], revCompStartPos);
			std::pair<size_t, size_t> otherForwardSorted = computeForwardStrandCoordinates(other.taxonCoords[i], revCompStartPos);
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
	for (size_t i = 0; i < taxonCoords.size(); ++i) {
		if (hasTaxon(i)) {
			size_t actSize = taxonCoords[i].second - taxonCoords[i].first;
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
	return taxonCoords[taxID].second - taxonCoords[taxID].first;
}

void Seed::increaseTaxonCoordsRight(size_t taxID) {
	taxonCoords[taxID].second++;
}
void Seed::decreaseTaxonCoordsLeft(size_t taxID) {
	assert(taxonCoords[taxID].first > 0);
	taxonCoords[taxID].first--;
}
void Seed::decreaseTaxonCoordRight(size_t taxID) {
	taxonCoords[taxID].second--;
}
void Seed::increaseTaxonCoordLeft(size_t taxID) {
	taxonCoords[taxID].first++;
}
