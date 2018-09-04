/*
 * seeded_block.cpp
 *
 *  Created on: Sep 4, 2018
 *      Author: Sarah Lutteropp
 */

#include "seeded_block.hpp"

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

size_t SeededBlock::getN() const {
	return n;
}

std::vector<std::pair<size_t, size_t> > SeededBlock::getTaxonCoords() const {
	return taxonCoords;
}

bool SeededBlock::hasTaxon(size_t taxID) const {
	return taxonCoords[taxID].first != std::string::npos;
}
