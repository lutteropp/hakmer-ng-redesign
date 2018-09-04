/*
 * block_helper_functions.cpp
 *
 *  Created on: Sep 4, 2018
 *      Author: Sarah Lutteropp
 */

#include "block_helper_functions.hpp"

std::string extractTaxonSequence(const AlignedBlock& block, size_t taxID) {
	std::string res;
	if (block.hasTaxon(taxID)) {
		return block.getAlignment()[taxID];
	}
	return res;
}

std::string extractTaxonSequence(const ExtendedBlock& block, size_t taxID, const std::string& T) {
	std::string res;
	if (block.hasTaxon(taxID)) {
		std::pair<size_t, size_t> coord = block.getTaxonCoordsWithFlanks(taxID);
		return T.substr(coord.first, coord.second + 1 - coord.first);
	}
	return res;
}

std::string extractTaxonSequence(const SeededBlock& block, size_t taxID, const std::string& T) {
	std::string res;
	if (block.hasTaxon(taxID)) {
		std::pair<size_t, size_t> coord = block.getTaxonCoords(taxID);
		return T.substr(coord.first, coord.second + 1 - coord.first);
	}
	return res;
}

std::string createMissingString(size_t len) {
	std::string res;
	for (size_t i = 0; i < len; ++i) {
		res += '?';
	}
	return res;
}
