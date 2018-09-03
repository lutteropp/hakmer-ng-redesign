/*
 * seeded_block.hpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <vector>
#include <cstdlib>
#include <string>

class SeededBlock {
public:
	SeededBlock(size_t nTax);
	void addTaxon(size_t taxID, size_t firstCoord, size_t lastCoord);
	size_t getN();
	std::vector<std::pair<size_t, size_t> > getTaxonCoords();
	bool hasTaxon(size_t taxID);
private:
	std::vector<std::pair<size_t, size_t> > taxonCoords;
	size_t n;
};

SeededBlock::SeededBlock(size_t nTax) {
	taxonCoords.resize(nTax);
	for (size_t i = 0; i < nTax; ++i) {
		taxonCoords[i].first = std::string::npos;
		taxonCoords[i].second = std::string::npos;
	}
	n = 0;
}

void SeededBlock::addTaxon(size_t taxID, size_t firstCoord, size_t lastCoord) {
	n++;
	taxonCoords[taxID].first = firstCoord;
	taxonCoords[taxID].second = lastCoord;
}

size_t SeededBlock::getN() {
	return n;
}

std::vector<std::pair<size_t, size_t> > SeededBlock::getTaxonCoords() {
	return taxonCoords;
}

bool SeededBlock::hasTaxon(size_t taxID) {
	return taxonCoords[taxID].first != std::string::npos;
}
