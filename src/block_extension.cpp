/*
 * block_extension.cpp
 *
 *  Created on: Nov 7, 2018
 *      Author: sarah
 */

#include "block_extension.hpp"
#include "seed.hpp"
#include "presence_checker.hpp"

#include "extended_block.hpp"
#include "block_helper_functions.hpp"
#include <algorithm>
#include <cassert>

bool canGoLeftAll(const Seed& block, const PresenceChecker& presenceChecker, size_t nTax, size_t offset) {
	bool canGo = true;
	for (size_t i = 0; i < nTax; ++i) {
		if (block.hasTaxon(i)) {
			if (!presenceChecker.isFree(block.getTaxonCoords(i).first - offset)) {
				canGo = false;
				break;
			}
		}
	}
	return canGo;
}

bool canGoRightAll(const Seed& block, const PresenceChecker& presenceChecker, size_t nTax, size_t offset) {
	bool canGo = true;
	for (size_t i = 0; i < nTax; ++i) {
		if (block.hasTaxon(i)) {
			if (!presenceChecker.isFree(block.getTaxonCoords(i).second + offset)) {
				canGo = false;
				break;
			}
		}
	}
	return canGo;
}

bool allLeftSame(const Seed& seededBlock, const std::string& T, const std::vector<size_t>& taxIDs) {
	char leftChar = T[seededBlock.getTaxonCoords(taxIDs[0]).first - 1];
	for (size_t i = 1; i < taxIDs.size(); ++i) {
		char actChar = T[seededBlock.getTaxonCoords(taxIDs[i]).first - 1];
		if (!ambiguousMatch(actChar, leftChar))
			return false;
	}
	return true;
}

bool allRightSame(const Seed& seededBlock, const std::string& T, const std::vector<size_t>& taxIDs) {
	char rightChar = T[seededBlock.getTaxonCoords(taxIDs[0]).second + 1];
	for (size_t i = 1; i < taxIDs.size(); ++i) {
		char actChar = T[seededBlock.getTaxonCoords(taxIDs[i]).second + 1];
		if (!ambiguousMatch(actChar, rightChar))
			return false;
	}
	return true;
}

bool allLeftSame(const Seed& seededBlock, const std::string& T, size_t nTax, size_t offset = 1) {
	std::vector<size_t> taxIDs = seededBlock.getTaxonIDsInBlock();
	char leftChar = T[seededBlock.getTaxonCoords(taxIDs[0]).first - offset];
	for (size_t i = 1; i < taxIDs.size(); ++i) {
		char actChar = T[seededBlock.getTaxonCoords(taxIDs[i]).first - offset];
		if (!ambiguousMatch(actChar, leftChar))
			return false;
	}
	return true;
}

bool allRightSame(const Seed& seededBlock, const std::string& T, size_t nTax, size_t offset = 1) {
	std::vector<size_t> taxIDs = seededBlock.getTaxonIDsInBlock();
	char rightChar = T[seededBlock.getTaxonCoords(taxIDs[0]).second + offset];
	for (size_t i = 1; i < taxIDs.size(); ++i) {
		char actChar = T[seededBlock.getTaxonCoords(taxIDs[i]).second + offset];
		if (!ambiguousMatch(actChar, rightChar))
			return false;
	}
	return true;
}

std::pair<size_t, size_t> computeBestCaseMaxSizes(Seed& seededBlock, const std::string& T, PresenceChecker& presenceChecker, size_t nTax) {
	size_t maxSizeLeft = 0;
	size_t maxSizeRight = 0;
	bool canLeft = true;
	bool canRight = true;
	size_t leftOffset = 1;
	size_t rightOffset = 1;
	while (canLeft || canRight) {
		if (!canGoLeftAll(seededBlock, presenceChecker, nTax, leftOffset)) {
			canLeft = false;
		}
		if (!canGoRightAll(seededBlock, presenceChecker, nTax, rightOffset)) {
			canRight = false;
		}
		if (canLeft && !allLeftSame(seededBlock, T, nTax, leftOffset)) {
			canLeft = false;
		}
		if (canRight && !allRightSame(seededBlock, T, nTax, rightOffset)) {
			canRight = false;
		}
		if (canLeft) {
			maxSizeLeft++;
			leftOffset++;
		}
		if (canRight) {
			maxSizeRight++;
			rightOffset++;
		}
	}
	return std::make_pair(maxSizeLeft, maxSizeRight);
}

void trivialExtensionPartial(Seed& seededBlock, const std::string& T, PresenceChecker& presenceChecker, size_t nTax) {
	// we perform trivial extension as long as at least 4 taxa are still availabe in the current direction
	std::vector<size_t> taxIDsLeft = seededBlock.getTaxonIDsInBlock();
	std::vector<size_t> taxIDsRight = seededBlock.getTaxonIDsInBlock();
	size_t leftok = seededBlock.getNTaxInBlock();
	size_t rightok = seededBlock.getNTaxInBlock();
	while (true) { // partial trivial left extension
		leftok = 0;
		std::vector<size_t> taxIDsToRemove;
		for (size_t tID : taxIDsLeft) {
			if (seededBlock.getTaxonCoords(tID).first > 0 && presenceChecker.isFree(seededBlock.getTaxonCoords(tID).first - 1)) {
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
				assert(seededBlock.getTaxonCoords(tID).first <= 2 * T.size());
				assert(seededBlock.getTaxonCoords(tID).second <= 2 * T.size());
				seededBlock.decreaseTaxonCoordsLeft(tID);
				presenceChecker.setTaken(seededBlock.getTaxonCoords(tID).first);
			}
		} else {
			break;
		}
	}

	while (true) { // partial trivial right extension
		rightok = 0;
		std::vector<size_t> taxIDsToRemove;
		for (size_t tID : taxIDsRight) {
			if (presenceChecker.isFree(seededBlock.getTaxonCoords(tID).second + 1)) {
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
				presenceChecker.setTaken(seededBlock.getTaxonCoords(tID).second);
			}
		} else {
			break;
		}
	}
}

void trivialExtension(Seed& seededBlock, const std::string& T, PresenceChecker& presenceChecker, size_t nTax, const Options& options) {
	if (!options.simpleExtension) {
		return trivialExtensionPartial(seededBlock, T, presenceChecker, nTax);
	}

	std::pair<size_t, size_t> bestCaseMaxSize = computeBestCaseMaxSizes(seededBlock, T, presenceChecker, nTax);
	for (size_t i = 1; i <= bestCaseMaxSize.first; ++i) {
		/*if (!canGoLeftAll(seededBlock, presenceChecker, nTax, 1)) {
		 break;
		 }*/
		seededBlock.decreaseAllTaxonCoordsLeft();
	}

	for (size_t i = 1; i <= bestCaseMaxSize.second; ++i) {
		/*if (!canGoRightAll(seededBlock, presenceChecker, nTax, 1)) {
		 break;
		 }*/
		seededBlock.increaseAllTaxonCoordsRight();
	}
}

void swap(std::vector<double>& vec, size_t i, size_t j) {
	double tmp = vec[i];
	vec[i] = vec[j];
	vec[j] = tmp;
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
		const Options& options, bool directionRight) {
	size_t bestSize = 0;
	double bestScore = 1.0;
	size_t nTaxBlock = block.getNTaxInBlock();
	std::vector<size_t> taxIDs = block.getTaxonIDsInBlock();
	size_t finalFlankSize = 0;
	for (size_t i = 1; i <= options.flankWidth; ++i) {
		if (directionRight) {
			if (!canGoRightAll(block, presenceChecker, nTax, i)) {
				break;
			}
		} else {
			if (!canGoLeftAll(block, presenceChecker, nTax, i)) {
				break;
			}
		}
		std::vector<char> charsToAdd(taxIDs.size());
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
			charsToAdd[j] = T[coord];
		}

		if (directionRight) {
			block.msaWrapper.addCharsRight(charsToAdd);
		} else {
			block.msaWrapper.addCharsLeft(charsToAdd);
		}
		finalFlankSize = i;
	}
	return finalFlankSize;
}

std::vector<char> findCharsToAdd(ExtendedBlock& block, const std::vector<size_t>& taxIDsPresent, const std::vector<size_t>& allTaxIDs,
		const std::string& T, bool leftDirection) {
	std::vector<char> charsToAdd;
	charsToAdd.resize(allTaxIDs.size());
	std::vector<bool> presence(allTaxIDs.size(), false);
	for (size_t tID : taxIDsPresent) {
		presence[tID] = true;
	}

	for (size_t i = 0; i < allTaxIDs.size(); ++i) {
		char c = '-';
		if (presence[i]) {
			size_t coord;
			if (leftDirection) {
				coord = block.getTaxonCoordsWithFlanks(allTaxIDs[i]).first - 1;
			} else {
				coord = block.getTaxonCoordsWithFlanks(allTaxIDs[i]).second + 1;
			}
			c = T[coord];
		}
		charsToAdd[i] = c;
	}
	return charsToAdd;
}

ExtendedBlock extendBlockPartial(ExtendedBlock& block, const std::string& T, size_t nTax, PresenceChecker& presenceChecker,
		const Options& options, size_t startingLeftFlankSize = 0, size_t startingRightFlankSize = 0) {
	size_t minTaxaToBeOk = options.minTaxaPerBlock;

	std::vector<size_t> taxIDsLeft = block.getTaxonIDsInBlock();
	std::vector<size_t> taxIDsRight = block.getTaxonIDsInBlock();
	size_t leftok = block.getNTaxInBlock();
	size_t rightok = block.getNTaxInBlock();
	size_t leftFlankSize = startingLeftFlankSize;
	size_t rightFlankSize = startingRightFlankSize;
	while (leftFlankSize < options.flankWidth) { // partial left extension
		leftok = 0;
		leftFlankSize++;
		std::vector<size_t> taxIDsToRemove;
		for (size_t tID : taxIDsLeft) {
			if (presenceChecker.isFree(block.getTaxonCoordsWithFlanks(tID).first - 1)) {
				leftok++;
			} else {
				taxIDsToRemove.push_back(tID);
			}
		}
		if (leftok < minTaxaToBeOk) {
			break;
		}
		for (size_t tID : taxIDsToRemove) {
			taxIDsLeft.erase(std::remove(taxIDsLeft.begin(), taxIDsLeft.end(), tID));
		}
		for (size_t tID : taxIDsLeft) {
			block.decrementLeftFlank(tID);
		}

		// we still need to add the chars to the MSA
		std::vector<char> charsToAdd = findCharsToAdd(block, taxIDsLeft, block.getTaxonIDsInBlock(), T, true);
		block.msaWrapper.addCharsLeft(charsToAdd);
	}

	while (rightFlankSize < options.flankWidth) { // partial trivial right extension
		rightFlankSize++;
		rightok = 0;
		std::vector<size_t> taxIDsToRemove;
		for (size_t tID : taxIDsRight) {
			if (presenceChecker.isFree(block.getTaxonCoordsWithFlanks(tID).second + 1)) {
				rightok++;
			} else {
				taxIDsToRemove.push_back(tID);
			}
		}
		if (rightok < minTaxaToBeOk) {
			break;
		}
		for (size_t tID : taxIDsToRemove) {
			taxIDsRight.erase(std::remove(taxIDsRight.begin(), taxIDsRight.end(), tID));
		}
		for (size_t tID : taxIDsRight) {
			block.incrementRightFlank(tID);
		}

		// we still need to add the chars to the MSA
		std::vector<char> charsToAdd = findCharsToAdd(block, taxIDsRight, block.getTaxonIDsInBlock(), T, false);
		block.msaWrapper.addCharsRight(charsToAdd);
	}

	return block;
}

ExtendedBlock extendBlock(const Seed& seededBlock, const std::string& T, size_t nTax, PresenceChecker& presenceChecker,
		const Options& options) {
	ExtendedBlock block(seededBlock, nTax, options.noIndels);

	if (options.maxMismatches == 0 && (!options.trimSeeds || !options.simpleTrimming) && options.simpleExtension) {
		std::string seedSequence = extractTaxonSequence(seededBlock, seededBlock.getTaxonIDsInBlock()[0], T);
		block.msaWrapper.init(block.getNTaxInBlock());
		block.msaWrapper.setSeeds(seedSequence);
	} else {
		std::vector<std::string> seedSequences;
		for (size_t i = 0; i < block.getTaxonIDsInBlock().size(); ++i) {
			seedSequences.push_back(extractTaxonSequence(seededBlock, seededBlock.getTaxonIDsInBlock()[i], T));
		}
		block.msaWrapper.init(block.getNTaxInBlock());
		block.msaWrapper.setSeeds(seedSequences);
	}

	size_t bestLeft = findPerfectFlankSize(block, nTax, presenceChecker, T, options, false);
	size_t bestRight = findPerfectFlankSize(block, nTax, presenceChecker, T, options, true);

	block.setLeftFlankSize(bestLeft);
	block.setRightFlankSize(bestRight);

	block.msaWrapper.shrinkDownToLeftFlank(bestLeft);
	block.msaWrapper.shrinkDownToRightFlank(bestRight);

	if (!options.simpleExtension) {
		return extendBlockPartial(block, T, nTax, presenceChecker, options, block.getAverageLeftFlankSize(),
				block.getAverageRightFlankSize());
	}

	block.msaWrapper.clearMSADataStructures();

	return block;
}
