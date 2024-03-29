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
			alwaysTaken.push_back(i);
		}
	}
}

PresenceChecker::PresenceChecker(const PresenceChecker& other) {
	this->revComp = other.revComp;
	this->freePos = other.freePos;
	this->alwaysTaken = other.alwaysTaken;
	this->nTax = other.nTax;
}

size_t PresenceChecker::mirrorCoord(size_t coord) const {
	size_t coordBefore = coord;
	if (freePos.size() == 0) {
		throw std::runtime_error("Why is freePos.size() == 0???");
	}
	if (revComp && coord >= freePos.size()) {
		coord = 2 * freePos.size() - coord - 1;
	}
	if (coord >= freePos.size()) {
		std::cout << "coordBefore: " << coordBefore << "\n";
		std::cout << "coord: " << coord << "\n";
		std::cout << "freePos.size(): " << freePos.size() << "\n";
		throw std::runtime_error("Coord is too large");
	}
	return coord;
}

bool PresenceChecker::isFree(size_t coord) const {
	coord = mirrorCoord(coord);
	return freePos[coord];
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
	coord = mirrorCoord(coord);
	freePos[coord] = false;
}

void PresenceChecker::setFree(size_t coord) {
	coord = mirrorCoord(coord);
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
	for (size_t i = 0; i < alwaysTaken.size(); ++i) {
		setTaken(alwaysTaken[i]);
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
