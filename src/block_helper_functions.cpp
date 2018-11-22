/*
 * block_helper_functions.cpp
 *
 *  Created on: Sep 4, 2018
 *      Author: Sarah Lutteropp
 */

#include "block_helper_functions.hpp"
#include "extended_block.hpp"
#include <iostream>
#include "seed.hpp"

std::string createMissingString(size_t len) {
	std::string res;
	for (size_t i = 0; i < len; ++i) {
		res += '?';
	}
	return res;
}

std::string extractTaxonSequence(ExtendedBlock& block, size_t taxID) {
	std::string res;
	size_t aliWidth = block.msaWrapper.getAlignmentWidth();
	std::vector<size_t> taxIDs = block.getTaxonIDsInBlock();
	if (block.hasTaxon(taxID)) {
		for (size_t i = 0; i < taxIDs.size(); ++i) {
			if (taxIDs[i] == taxID) {
				res = block.msaWrapper.assembleMSA()[i];
				break;
			}
		}
	} else {
		res = createMissingString(aliWidth);
	}
	return res;
}


std::string extractTaxonSequence(const Seed& block, size_t taxID, const std::string& T) {
	std::string res;
	if (block.hasTaxon(taxID)) {
		std::pair<size_t, size_t> coord = block.getTaxonCoords(taxID);
		res = T.substr(coord.first, coord.second + 1 - coord.first);
	} else {
		res = createMissingString(block.getAverageSeedSize());
	}
	return res;
}
