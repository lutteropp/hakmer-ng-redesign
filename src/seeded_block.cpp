/*
 * seeded_block.cpp
 *
 *  Created on: Sep 4, 2018
 *      Author: Sarah Lutteropp
 */

#include "seeded_block.hpp"
#include <stdexcept>
#include <iostream>

SeededBlock::SeededBlock(size_t nTax) {
	taxonCoords.resize(nTax);
	for (size_t i = 0; i < nTax; ++i) {
		taxonCoords[i].first = std::string::npos;
		taxonCoords[i].second = std::string::npos;
	}
	n = 0;
	k = 0;
}

void SeededBlock::addTaxon(size_t taxID, size_t firstCoord, size_t lastCoord) {
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
}

size_t SeededBlock::getNTaxInBlock() const {
	return n;
}

std::pair<size_t, size_t> SeededBlock::getTaxonCoords(size_t taxID) const {
	return taxonCoords[taxID];
}

bool SeededBlock::hasTaxon(size_t taxID) const {
	return taxonCoords[taxID].first != std::string::npos;
}

size_t SeededBlock::getSeedSize() const {
	return k;
}

std::vector<size_t> SeededBlock::getTaxonIDsInBlock() const {
	std::vector<size_t> res;
	for (size_t i = 0; i < taxonCoords.size(); ++i) {
		if (hasTaxon(i)) {
			res.push_back(i);
		}
	}
	return res;
}
