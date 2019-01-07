/*
 * presence_checker.cpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#include "presence_checker.hpp"

PresenceChecker::PresenceChecker(const IndexedConcatenatedSequence& concat, bool revComp) :
		revComp(revComp) {
	size_t nSites = concat.getConcatSize();
	if (revComp) {
		nSites /= 2;
	}
	freePos.resize(nSites);
	for (size_t i = 0; i < nSites; ++i) {
		freePos[i] = true;
	}
	this->nTax = concat.nTax();

	for (size_t i = 0; i < nSites; ++i) {
		if (concat.getConcatenatedSeq()[i] == '$') {
			setTaken(i);
		}
	}
	size = freePos.size();
}

PresenceChecker::PresenceChecker(const PresenceChecker& other) {
	this->revComp = other.revComp;
	this->freePos.resize(other.freePos.size());
	for (size_t i = 0; i < other.freePos.size(); ++i) {
		this->freePos[i] = other.freePos[i];
	}
	this->nTax = other.nTax;
	size = freePos.size();
}

bool PresenceChecker::isFree(size_t coord) const {
	if (revComp && coord >= size) {
		coord = 2*size - coord - 1; // TODO: This looks wrong.
	}
	if (coord >= size) {
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
	if (revComp && coord >= size) {
		coord = 2*size - coord - 1;
	}
	freePos[coord] = false;
}

void PresenceChecker::setFree(size_t coord) {
	if (revComp && coord >= size) {
		coord = 2*size - coord - 1;
	}
	freePos[coord] = true;
}

void PresenceChecker::setTaken(size_t firstCoord, size_t lastCoord) {
	for (size_t i = firstCoord; i <= lastCoord; ++i) {
		setTaken(i);
	}
}

void PresenceChecker::setFree(size_t firstCoord, size_t lastCoord) {
	for (size_t i = firstCoord; i <= lastCoord; ++i) {
		setFree(i);
	}
}

void PresenceChecker::reserveSeededBlock(const Seed& block) {
	for (size_t i = 0; i < nTax; ++i) {
		if (block.hasTaxon(i)) {
			SimpleCoords coords = block.getSeedCoords(i);
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

void PresenceChecker::freeExtendedBlock(const ExtendedBlock& block) {
	for (size_t i = 0; i < nTax; ++i) {
		if (block.hasTaxon(i)) {
			std::pair<size_t, size_t> coords = block.getTaxonCoordsWithFlanks(i);
			setFree(coords.first, coords.second);
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

bool PresenceChecker::isFineWithoutSeed(const ExtendedBlock& block) {
	for (size_t i = 0; i < nTax; ++i) {
		if (block.hasTaxon(i)) {
			std::pair<size_t, size_t> coords = block.getTaxonCoordsWithFlanks(i);
			std::pair<size_t, size_t> seedCoords = block.getTaxonCoordsWithoutFlanks(i);
			if (seedCoords.first - 1 >= coords.first && !isFree(coords.first, seedCoords.first - 1)) {
				return false;
			}
			if (coords.second >= seedCoords.second + 1 && !isFree(seedCoords.second + 1, coords.second)) {
				return false;
			}
		}
	}
	return true;
}

bool PresenceChecker::isFine(const Seed& block) {
	for (size_t i = 0; i < nTax; ++i) {
		if (block.hasTaxon(i)) {
			SimpleCoords coords = block.getSeedCoords(i);
			if (!isFree(coords.first, coords.second)) {
				return false;
			}
		}
	}
	return true;
}
