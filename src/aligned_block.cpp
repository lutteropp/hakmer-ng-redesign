/*
 * aligned_block.cpp
 *
 *  Created on: Sep 4, 2018
 *      Author: Sarah Lutteropp
 */

#include <stdexcept>

#include "aligned_block.hpp"
#include "block_helper_functions.hpp"

AlignedBlock::AlignedBlock(const ExtendedBlock& extendedBlock, size_t nTax) :
		aligned(false), myBlock(extendedBlock), nTax(nTax) {
}

bool AlignedBlock::isAligned() const {
	return aligned;
}

bool AlignedBlock::hasTaxon(size_t taxID) const {
	return myBlock.hasTaxon(taxID);
}

void AlignedBlock::align(const std::string& T, const Options& options) {
	if (aligned)
		return;
	std::vector<size_t> taxIDs = getTaxonIDsInBlock();
	for (size_t i = 0; i < taxIDs.size(); ++i) {
		alignment.push_back(extractTaxonSequence(myBlock, taxIDs[i], T));
	}
	if (!options.noIndels) {
		throw std::runtime_error("Actual alignment is not implemented yet.");
		// TODO: Implement the actual alignment, this is, insert gaps into the sequences just stored in the alignment variable
	}
	aligned = true;
}

std::vector<std::string> AlignedBlock::getAlignment() const {
	return alignment;
}

size_t AlignedBlock::getAlignmentWidth() const {
	for (size_t i = 0; i < alignment.size(); ++i) {
		if (alignment[i].size() > 0) {
			return alignment[i].size();
		}
	}
	return 0;
}

std::pair<size_t, size_t> AlignedBlock::getTaxonCoordsWithFlanks(size_t taxID) const {
	return myBlock.getTaxonCoordsWithFlanks(taxID);
}

std::vector<size_t> AlignedBlock::getTaxonIDsInBlock() const {
	return myBlock.getTaxonIDsInBlock();
}
