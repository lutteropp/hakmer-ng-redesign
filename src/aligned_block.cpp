/*
 * aligned_block.cpp
 *
 *  Created on: Sep 4, 2018
 *      Author: Sarah Lutteropp
 */

#include "aligned_block.hpp"

AlignedBlock::AlignedBlock(const ExtendedBlock& extendedBlock, size_t nTax) : aligned(false), myBlock(extendedBlock), nTax(nTax) {}

bool AlignedBlock::isAligned() const {
	return aligned;
}

bool AlignedBlock::hasTaxon(size_t taxID) const {
	return myBlock.hasTaxon(taxID);
}

void AlignedBlock::align() {
	if (aligned) return;
	// TODO: Implement the actual alignment
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
