/*
 * presence_checker.hpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <cstdlib>
#include <vector>

#include "seeded_block.hpp"
#include "extended_block.hpp"

// TODO: Make this thing thread-safe

class PresenceChecker {
public:
	PresenceChecker(size_t nSites, size_t nTax);
	bool isFree(size_t firstCoord, size_t lastCoord);
	bool isFree(size_t coord);
	void setTaken(size_t firstCoord, size_t lastCoord);
	void setTaken(size_t coord);
	void reserveSeededBlock(const SeededBlock& block);
	void reserveExtendedBlock(const ExtendedBlock& block);

	bool isFine(const ExtendedBlock& block);
private:
	std::vector<bool> freePos;
	size_t nTax;
};
