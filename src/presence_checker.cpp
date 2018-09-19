/*
 * presence_checker.cpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#include "presence_checker.hpp"

PresenceChecker::PresenceChecker(const IndexedConcatenatedSequence& concat, bool revComp) {
	size_t nSites = concat.getConcatSize();
	if (revComp) {
		nSites /= 2;
	}
	freePos.resize(nSites);
	for (size_t i = 0; i < nSites; ++i) {
		freePos[i] = true;
	}
	this->nTax = concat.nTax();
	this->revComp = revComp;

	for (size_t i = 0; i < nSites; ++i) {
		if (concat.getConcatenatedSeq()[i] == '$') {
			setTaken(i);
		}
	}

}

bool PresenceChecker::isFree(size_t coord) const {
	if (coord >= freePos.size() && revComp) {
		coord = 2 * freePos.size() - 1 - coord;
	}
	if (coord >= freePos.size()) {
		return false;
	} else {
		return freePos[coord];
	}
}

bool PresenceChecker::isFree(size_t firstCoord, size_t lastCoord) const {
	for (size_t i = firstCoord; i <= lastCoord; ++i) {
		if (!isFree(i)) {
			return false;
		}
	}
	return true;
}

void PresenceChecker::setTaken(size_t coord) {
	if (coord >= freePos.size() && revComp) {
		coord = 2 * freePos.size() - 1 - coord;
	}
	freePos[coord] = false;
}

void PresenceChecker::setTaken(size_t firstCoord, size_t lastCoord) {
	for (size_t i = firstCoord; i <= lastCoord; ++i) {
		setTaken(i);
	}
}

void PresenceChecker::reserveSeededBlock(const SeededBlock& block) {
	for (size_t i = 0; i < nTax; ++i) {
		if (block.hasTaxon(i)) {
			std::pair<size_t, size_t> coords = block.getTaxonCoords(i);
			setTaken(coords.first, coords.second);
		}
	}
}

void PresenceChecker::reserveExtendedBlock(const ExtendedBlock& block) {
	for (size_t i = 0; i < nTax; ++i) {
		if (block.hasTaxon(i)) {
			std::pair<size_t, size_t> coords = block.getTaxonCoordsWithFlanks(i);
			setTaken(coords.first, coords.second);
		}
	}
}

bool PresenceChecker::isFine(const ExtendedBlock& block) {
	for (size_t i = 0; i < nTax; ++i) {
		if (block.hasTaxon(i)) {
			std::pair<size_t, size_t> coords = block.getTaxonCoordsWithFlanks(i);
			if (!isFree(coords.first, coords.second)) {
				return false;
			}
		}
	}
	return true;
}
