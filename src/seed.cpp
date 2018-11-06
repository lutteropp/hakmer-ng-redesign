/*
 * seeded_block.cpp
 *
 *  Created on: Sep 4, 2018
 *      Author: Sarah Lutteropp
 */

#include <stdexcept>
#include <iostream>
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
	return taxonCoords[taxID].first != std::string::npos;
}

size_t Seed::getSeedSize() const {
	return k;
}

std::vector<size_t> Seed::getTaxonIDsInBlock() const {
	return taxIDs;
}

void Seed::increaseTaxonCoordsRight() {
	for (size_t i = 0; i < taxIDs.size(); ++i) {
		taxonCoords[taxIDs[i]].second++;
	}
	k++;
}

void Seed::decreaseTaxonCoordsLeft() {
	for (size_t i = 0; i < taxIDs.size(); ++i) {
		taxonCoords[taxIDs[i]].first--;
	}
	k++;
}
