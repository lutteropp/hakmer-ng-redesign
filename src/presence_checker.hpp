/*
 * presence_checker.hpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <stddef.h>
#include <vector>

#include "extended_block.hpp"
#include "indexed_concat.hpp"
#include "seed.hpp"

// TODO: Make this thing thread-safe

class PresenceChecker {
public:
	PresenceChecker(const IndexedConcatenatedSequence& concat, bool revComp);
	PresenceChecker(const PresenceChecker& other);
	bool isFree(size_t firstCoord, size_t lastCoord) const;
	bool isFree(size_t coord) const;
	void setTaken(size_t firstCoord, size_t lastCoord);
	void setTaken(size_t coord);
	void setFree(size_t firstCoord, size_t lastCoord);
	void setFree(size_t coord);
	void reserveSeededBlock(const Seed& block);
	void reserveExtendedBlock(const ExtendedBlock& block);
	void freeExtendedBlock(const ExtendedBlock& block);

	bool isFine(const ExtendedBlock& block);
	bool isFineWithoutSeed(const ExtendedBlock& block);
	bool isFine(const Seed& block);
private:
	std::vector<bool> freePos;
	size_t nTax;
	bool revComp;
	size_t size;
};
