/*
 * block_extension.cpp
 *
 *  Created on: Nov 7, 2018
 *      Author: sarah
 */

#include "block_extension.hpp"

#include <stddef.h>
#include <algorithm>
#include <cassert>
#include <iterator>
#include <stdexcept>

#include "alignment/simple_coords.hpp"
#include "dna_functions.hpp"

bool canGoLeftAll(const Seed& block, const PresenceChecker& presenceChecker, size_t nTax, size_t offset = 1) {
	bool canGo = true;
	for (size_t i = 0; i < nTax; ++i) {
		if (block.hasTaxon(i)) {
			if (!presenceChecker.isFree(block.getSeedCoords(i).first - offset)) {
				canGo = false;
				break;
			}
		}
	}
	return canGo;
}

bool canGoRightAll(const Seed& block, const PresenceChecker& presenceChecker, size_t nTax, size_t offset = 1) {
	bool canGo = true;
	for (size_t i = 0; i < nTax; ++i) {
		if (block.hasTaxon(i)) {
			if (!presenceChecker.isFree(block.getSeedCoords(i).second + offset)) {
				canGo = false;
				break;
			}
		}
	}
	return canGo;
}

bool allLeftSame(const Seed& seededBlock, const std::string& T, const std::vector<size_t>& taxIDs) {
	char leftChar = T[seededBlock.getSeedCoords(taxIDs[0]).first - 1];
	for (size_t i = 1; i < taxIDs.size(); ++i) {
		char actChar = T[seededBlock.getSeedCoords(taxIDs[i]).first - 1];
		if (!ambiguousMatch(actChar, leftChar))
			return false;
	}
	return true;
}

bool allRightSame(const Seed& seededBlock, const std::string& T, const std::vector<size_t>& taxIDs) {
	char rightChar = T[seededBlock.getSeedCoords(taxIDs[0]).second + 1];
	for (size_t i = 1; i < taxIDs.size(); ++i) {
		char actChar = T[seededBlock.getSeedCoords(taxIDs[i]).second + 1];
		if (!ambiguousMatch(actChar, rightChar))
			return false;
	}
	return true;
}

bool allLeftSame(const Seed& seededBlock, const std::string& T, size_t nTax, size_t offset = 1) {
	std::vector<size_t> taxIDs = seededBlock.getTaxonIDsInBlock();
	char leftChar = T[seededBlock.getSeedCoords(taxIDs[0]).first - offset];
	for (size_t i = 1; i < taxIDs.size(); ++i) {
		char actChar = T[seededBlock.getSeedCoords(taxIDs[i]).first - offset];
		if (!ambiguousMatch(actChar, leftChar))
			return false;
	}
	return true;
}

bool allRightSame(const Seed& seededBlock, const std::string& T, size_t nTax, size_t offset = 1) {
	std::vector<size_t> taxIDs = seededBlock.getTaxonIDsInBlock();
	char rightChar = T[seededBlock.getSeedCoords(taxIDs[0]).second + offset];
	for (size_t i = 1; i < taxIDs.size(); ++i) {
		char actChar = T[seededBlock.getSeedCoords(taxIDs[i]).second + offset];
		if (!ambiguousMatch(actChar, rightChar))
			return false;
	}
	return true;
}

void trivialExtensionPartial(Seed& seededBlock, const std::string& T, PresenceChecker& presenceChecker, size_t nTax) {
	// we perform trivial extension as long as at least 4 taxa are still available in the current direction
	std::vector<size_t> taxIDsLeft = seededBlock.getTaxonIDsInBlock();
	std::vector<size_t> taxIDsRight = seededBlock.getTaxonIDsInBlock();
	size_t leftok = seededBlock.getNTaxInBlock();
	size_t rightok = seededBlock.getNTaxInBlock();
	while (true) { // partial trivial left extension
		leftok = 0;
		std::vector<size_t> taxIDsToRemove;
		for (size_t tID : taxIDsLeft) {
			if (seededBlock.getSeedCoords(tID).first > 0 && presenceChecker.isFree(seededBlock.getSeedCoords(tID).first - 1)) {
				leftok++;
			} else {
				taxIDsToRemove.push_back(tID);
			}
		}
		if (leftok < 4) {
			break;
		}
		for (size_t tID : taxIDsToRemove) {
			taxIDsLeft.erase(std::remove(taxIDsLeft.begin(), taxIDsLeft.end(), tID));
		}
		if (allLeftSame(seededBlock, T, taxIDsLeft)) {
			for (size_t tID : taxIDsLeft) {
				assert(seededBlock.getSeedCoords(tID).first <= 2 * T.size());
				assert(seededBlock.getSeedCoords(tID).second <= 2 * T.size());
				seededBlock.decreaseTaxonCoordsLeft(tID);
			}
			for (size_t tID : seededBlock.getTaxonIDsInBlock()) {
				if (std::find(taxIDsLeft.begin(), taxIDsLeft.end(), tID) == taxIDsLeft.end()) {
					seededBlock.addGapLeft(tID);
				}
			}
		} else {
			break;
		}
	}

	while (true) { // partial trivial right extension
		rightok = 0;
		std::vector<size_t> taxIDsToRemove;
		for (size_t tID : taxIDsRight) {
			if (presenceChecker.isFree(seededBlock.getSeedCoords(tID).second + 1)) {
				rightok++;
			} else {
				taxIDsToRemove.push_back(tID);
			}
		}
		if (rightok < 4) {
			break;
		}
		for (size_t tID : taxIDsToRemove) {
			taxIDsRight.erase(std::remove(taxIDsRight.begin(), taxIDsRight.end(), tID));
		}
		if (allRightSame(seededBlock, T, taxIDsRight)) {
			for (size_t tID : taxIDsRight) {
				seededBlock.increaseTaxonCoordsRight(tID);
			}
			for (size_t tID : seededBlock.getTaxonIDsInBlock()) {
				if (std::find(taxIDsRight.begin(), taxIDsRight.end(), tID) == taxIDsRight.end()) {
					seededBlock.addGapRight(tID);
				}
			}
		} else {
			break;
		}
	}
}

void trivialExtension(Seed& seededBlock, const std::string& T, PresenceChecker& presenceChecker, size_t nTax, const Options& options) {
	while (canGoLeftAll(seededBlock, presenceChecker, nTax) && allLeftSame(seededBlock, T, nTax)) {
		seededBlock.decreaseAllTaxonCoordsLeft();
	}
	while (canGoRightAll(seededBlock, presenceChecker, nTax) && allRightSame(seededBlock, T, nTax)) {
		seededBlock.increaseAllTaxonCoordsRight();
	}
	if (!options.simpleExtension) {
		return trivialExtensionPartial(seededBlock, T, presenceChecker, nTax);
	}
}

bool canGoLeftAll(const ExtendedBlock& block, const PresenceChecker& presenceChecker, size_t nTax, size_t offset) {
	bool canGo = true;
	for (size_t i = 0; i < nTax; ++i) {
		if (block.hasTaxon(i)) {
			if (!presenceChecker.isFree(block.getTaxonCoordsWithoutFlanks(i).first - offset)) {
				canGo = false;
				break;
			}
		}
	}
	return canGo;
}

bool canGoRightAll(const ExtendedBlock& block, const PresenceChecker& presenceChecker, size_t nTax, size_t offset) {
	bool canGo = true;
	for (size_t i = 0; i < nTax; ++i) {
		if (block.hasTaxon(i)) {
			if (!presenceChecker.isFree(block.getTaxonCoordsWithoutFlanks(i).second + offset)) {
				canGo = false;
				break;
			}
		}
	}
	return canGo;
}

size_t findPerfectFlankSize(ExtendedBlock& block, size_t nTax, const PresenceChecker& presenceChecker, const std::string& T,
		const Options& options, bool directionRight, size_t flankWidth) {
	size_t bestSize = 0;
	double bestScore = 1.0;
	size_t nTaxBlock = block.getNTaxInBlock();
	std::vector<size_t> taxIDs = block.getTaxonIDsInBlock();
	size_t finalFlankSize = 0;
	for (size_t i = 1; i <= flankWidth; ++i) {
		if (directionRight) {
			if (!canGoRightAll(block, presenceChecker, nTax, i)) {
				break;
			}
		} else {
			if (!canGoLeftAll(block, presenceChecker, nTax, i)) {
				break;
			}
		}
		for (size_t j = 0; j < taxIDs.size(); ++j) {
			size_t coord;
			if (directionRight) {
				coord = block.getTaxonCoordsWithoutFlanks(taxIDs[j]).second + i;
			} else {
				coord = block.getTaxonCoordsWithoutFlanks(taxIDs[j]).first - i;
			}
			if (coord >= T.size()) {
				throw std::runtime_error("This should not happen! Coord is too large.");
			}
		}
		finalFlankSize = i;
	}
	return finalFlankSize;
}

void extendBlockPartial(ExtendedBlock& block, const std::string& T, size_t nTax, PresenceChecker& presenceChecker, const Options& options,
		bool leftDirection, size_t flankWidth, size_t startingFlankSize = 0) {
	size_t minTaxaToBeOk = options.minTaxaPerBlock;

	std::vector<size_t> taxIDsLeft = block.getTaxonIDsInBlock();
	std::vector<bool> stillOk(nTax, false);
	for (size_t tID : taxIDsLeft) {
		stillOk[tID] = true;
	}

	size_t leftok = block.getNTaxInBlock();
	size_t flankSize = startingFlankSize;
	while (flankSize < flankWidth) {
		leftok = 0;
		flankSize++;
		std::vector<size_t> taxIDsToRemove;
		for (size_t tID : taxIDsLeft) {
			size_t coord;
			if (leftDirection) {
				coord = block.getTaxonCoordsWithFlanks(tID).first - 1;
			} else {
				coord = block.getTaxonCoordsWithFlanks(tID).second + 1;
			}
			if (presenceChecker.isFree(coord)) {
				leftok++;
			} else {
				taxIDsToRemove.push_back(tID);
				stillOk[tID] = false;
			}
		}
		if (leftok < minTaxaToBeOk) {
			break;
		}
		for (size_t tID : taxIDsToRemove) {
			taxIDsLeft.erase(std::remove(taxIDsLeft.begin(), taxIDsLeft.end(), tID));
		}

		for (size_t i = 0; i < nTax; ++i) {
			if (!block.hasTaxon(i)) {
				continue;
			}
			if (stillOk[i]) {
				if (leftDirection) {
					block.growLeftFlank(i);
				} else {
					block.growRightFlank(i);
				}
			} else {
				if (leftDirection) {
					block.addLeftFlankGap(i);
				} else {
					block.addRightFlankGap(i);
				}
			}
		}
	}
}

ExtendedBlock extendBlock(const Seed& seededBlock, const std::string& T, size_t nTax, PresenceChecker& presenceChecker, size_t flankWidth,
		const Options& options) {
	ExtendedBlock block(seededBlock, nTax);

	size_t bestLeft = findPerfectFlankSize(block, nTax, presenceChecker, T, options, false, flankWidth);
	size_t bestRight = findPerfectFlankSize(block, nTax, presenceChecker, T, options, true, flankWidth);

	block.setLeftFlankSize(bestLeft);
	block.setRightFlankSize(bestRight);

	if (!options.simpleExtension) {
		extendBlockPartial(block, T, nTax, presenceChecker, options, true, flankWidth, block.getAverageLeftFlankSize());
		extendBlockPartial(block, T, nTax, presenceChecker, options, false, flankWidth, block.getAverageRightFlankSize());
	}
	return block;
}
