/*
 * block_helper_functions.cpp
 *
 *  Created on: Sep 4, 2018
 *      Author: Sarah Lutteropp
 */

#include "block_helper_functions.hpp"
#include "distance_estimator.hpp"

std::string createMissingString(size_t len) {
	std::string res;
	for (size_t i = 0; i < len; ++i) {
		res += '?';
	}
	return res;
}

std::string extractTaxonSequence(const AlignedBlock& block, size_t taxID) {
	std::string res;
	size_t aliWidth = block.getAlignmentWidth();
	if (block.hasTaxon(taxID)) {
		res = block.getAlignment()[taxID];
	} else {
		res = createMissingString(aliWidth);
	}
	return res;
}

std::string extractTaxonSequence(const ExtendedBlock& block, size_t taxID, const std::string& T) {
	std::string res;
	if (block.hasTaxon(taxID)) {
		std::pair<size_t, size_t> coord = block.getTaxonCoordsWithFlanks(taxID);
		res = T.substr(coord.first, coord.second + 1 - coord.first);
	} else {
		res = createMissingString(block.getMaxLeftFlankSize() + block.getSeedSize() + block.getMaxRightFlankSize());
	}
	return res;
}

std::string extractTaxonSequence(const SeededBlock& block, size_t taxID, const std::string& T) {
	std::string res;
	if (block.hasTaxon(taxID)) {
		std::pair<size_t, size_t> coord = block.getTaxonCoords(taxID);
		res = T.substr(coord.first, coord.second + 1 - coord.first);
	} else {
		res = createMissingString(block.getSeedSize());
	}
	return res;
}

std::vector<double> computePairwiseDistances(const AlignedBlock& block) {
	std::vector<double> res;
	std::vector<size_t> taxIDs = block.getTaxonIDsInBlock();
	for (size_t i = 0; i < taxIDs.size(); ++i) {
		for (size_t j = i + 1; j < taxIDs.size(); ++j) {
			HammingDistanceEstimator hamm(extractTaxonSequence(block, taxIDs[i]), extractTaxonSequence(block, taxIDs[j]));
			res.push_back(hamm.distance());
		}
	}
	return res;
}
