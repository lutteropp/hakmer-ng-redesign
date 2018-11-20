/*
 * block_extension.cpp
 *
 *  Created on: Nov 7, 2018
 *      Author: sarah
 */

#include "block_extension.hpp"
#include "seed.hpp"
#include "presence_checker.hpp"

#include "hmm_classification.hpp"
#include "block_helper_functions.hpp"

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

void trivialExtension(Seed& seededBlock, const std::string& T, PresenceChecker& presenceChecker, size_t nTax) {
	// first, perform trivial extension of the seeded block
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

double deltaScore(const std::array<double, 6>& pairwiseDist) {
	// This happens if two sequences are too divergent, this is, relative uncorrected distance > 0.75
	for (size_t i = 0; i < pairwiseDist.size(); ++i) {
		if (std::isnan(pairwiseDist[i]) || std::isinf(pairwiseDist[i])) {
			return 1;
		}
	}

	// Check for 3 or more identical sequences, as this leads to star topology
	if (pairwiseDist[0] == 0 && pairwiseDist[1] == 0) { // d(a,b) and d(a,c) -> a,b,c are equal
		return 1;
	} else if (pairwiseDist[0] == 0 && pairwiseDist[2] == 0) { // d(a,b) and d(a,d) -> a,b,d are equal
		return 1;
	} else if (pairwiseDist[1] == 0 && pairwiseDist[2] == 0) { // d(a,c) and d(a,d) -> a,c,d are equal
		return 1;
	} else if (pairwiseDist[3] == 0 && pairwiseDist[4] == 0) { // d(b,c) and d(b,d) -> b,c,d are equal
		return 1;
	}

	double d12_34 = pairwiseDist[0] + pairwiseDist[5];
	double d13_24 = pairwiseDist[1] + pairwiseDist[4];
	double d14_23 = pairwiseDist[2] + pairwiseDist[3];

	std::vector<double> distances = { d12_34, d13_24, d14_23 };
	if (distances[0] > distances[1])
		swap(distances, 0, 1);
	if (distances[0] > distances[2])
		swap(distances, 0, 2);
	if (distances[1] > distances[2])
		swap(distances, 1, 2);

	double score;
	if (distances[0] == distances[2]) { // star topology
		score = 1.0;
	} else {
		score = (distances[2] - distances[1]) / (distances[2] - distances[0]);
	}
	if (std::isnan(score)) {
		throw std::runtime_error("The delta score is nan!!!\n");
	}
	return score;
}

double averageDeltaScore(ExtendedBlock& block, size_t nTaxBlock, const Options& options) {
	double sum = 0.0;
	size_t numQuartets;

	if (options.quickDelta) {
		numQuartets = nTaxBlock - 3;
		for (size_t i = 0; i < nTaxBlock - 3; ++i) {
			size_t j = i + 1;
			size_t k = i + 2;
			size_t l = i + 3;
			std::array<double, 6> pairwiseDist = { block.getPairwiseNormalizedDistance(i, j, options), block.getPairwiseNormalizedDistance(
					i, k, options), block.getPairwiseNormalizedDistance(i, l, options), block.getPairwiseNormalizedDistance(j, k, options),
					block.getPairwiseNormalizedDistance(j, l, options), block.getPairwiseNormalizedDistance(k, l, options) };
			sum += deltaScore(pairwiseDist);
		}
	} else {
		numQuartets = nTaxBlock * (nTaxBlock - 1) * (nTaxBlock - 2) * (nTaxBlock - 3) / 24;
		for (size_t i = 0; i < nTaxBlock - 3; ++i) {
			for (size_t j = i + 1; j < nTaxBlock - 2; ++j) {
				for (size_t k = j + 1; k < nTaxBlock - 1; ++k) {
					for (size_t l = k + 1; l < nTaxBlock; ++l) {
						std::array<double, 6> pairwiseDist = { block.getPairwiseNormalizedDistance(i, j, options),
								block.getPairwiseNormalizedDistance(i, k, options), block.getPairwiseNormalizedDistance(i, l, options),
								block.getPairwiseNormalizedDistance(j, k, options), block.getPairwiseNormalizedDistance(j, l, options),
								block.getPairwiseNormalizedDistance(k, l, options) };
						sum += deltaScore(pairwiseDist);
					}
				}
			}
		}
	}
	return sum / numQuartets;
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

std::pair<size_t, double> findPerfectFlankSizeStatic(ExtendedBlock& block, size_t nTax, const PresenceChecker& presenceChecker,
		const std::string& T, const Options& options, bool directionRight) {
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
	double score = averageDeltaScore(block, nTaxBlock, options);
	return std::make_pair(finalFlankSize, score);
}

std::pair<size_t, double> findPerfectFlankSizeDelta(ExtendedBlock& block, size_t nTax, const PresenceChecker& presenceChecker,
		const std::string& T, const Options& options, bool directionRight) {
	size_t bestSize = 0;
	double bestScore = 1.0;
	size_t nTaxBlock = block.getNTaxInBlock();
	std::vector<size_t> taxIDs = block.getTaxonIDsInBlock();

	for (size_t i = 1; i <= options.maximumExtensionWidth; ++i) {
		if (bestScore <= options.maxDelta) {
			break;
		}
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

		double score = averageDeltaScore(block, nTaxBlock, options);
		if (score < bestScore) {
			bestScore = score;
			bestSize = i;
			//std::cout << "found a better score: " << score << "\n" << " with flank size: " << i << "\n";
		} else {
			if (score > bestScore && i - bestSize >= options.earlyStopCount) {
				break;
			}
		}
	}

	return std::make_pair(bestSize, bestScore);
}

// TODO: FIX THIS!!!
std::pair<size_t, double> findPerfectFlankSizeHMM(ExtendedBlock& block, size_t nTax, const PresenceChecker& presenceChecker,
		const std::string& T, const Options& options, bool directionRight) {
	size_t bestSize = 0;
	size_t nTaxBlock = block.getNTaxInBlock();
	std::vector<size_t> taxIDs = block.getTaxonIDsInBlock();

	Params hmmParams = prepareHmmParams(T, options);

	for (size_t i = 1; i <= options.maximumExtensionWidth; ++i) {
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
		bestSize++;
	}

	// Run the HomologyHMM on the pairwise alignments to check if we should already stop the extension...
	if (bestSize > 0) {
		if (isEntireMSAHomologous(block.msaWrapper, hmmParams)) {
			std::cout << "Homo\n";
		} else {
			std::cout << "No homo ;-P\n";
		}

		/*for (size_t t1 = 0; t1 < block.getNTaxInBlock(); t1++) {
		 for (size_t t2 = t1 + 1; t2 < block.getNTaxInBlock(); t2++) {
		 size_t goodSites;
		 if (directionRight) {
		 goodSites = findNumGoodSites(block.msaWrapper.getRightFlankAlignment(t1, t2).first,
		 block.msaWrapper.getRightFlankAlignment(t1, t2).second, hmmParams);
		 } else {
		 goodSites = findNumGoodSites(block.msaWrapper.getReversedLeftFlankAlignment(t1, t2).first,
		 block.msaWrapper.getReversedLeftFlankAlignment(t1, t2).second, hmmParams);
		 }
		 bestSize = std::min(bestSize, goodSites);
		 }
		 }*/
	}

	/*std::pair<size_t, size_t> seedCoords; // TODO: This currently only works with exactly matching seeds I think...
	 std::string seed;
	 for (size_t i = 0; i < nTax; ++i) {
	 if (block.hasTaxon(i)) { // TODO: ALso, this is super sloooow!!!
	 seed = T.substr(block.getTaxonCoordsWithoutFlanks(i).first, block.getSeedSize());
	 seedCoords.first = block.msaWrapper.assembleMSA()[0].find(seed);
	 seedCoords.second = seedCoords.first + block.getSeedSize() - 1;
	 break;
	 }
	 }
	 if (!directionRight && seedCoords.first == 0) {
	 bestSize = 0;
	 } else {
	 bestSize = findNumGoodSitesMSA(block.msaWrapper, directionRight, seedCoords, hmmParams);
	 }*/

	block.msaWrapper.disassembleMSA();

	return std::make_pair(bestSize, 1);
}

std::pair<size_t, double> findPerfectFlankSizeBigGaps(ExtendedBlock& block, size_t nTax, const PresenceChecker& presenceChecker,
		const std::string& T, const Options& options, bool directionRight) {
	throw std::runtime_error("Not implemented yet");
}

std::pair<size_t, double> findPerfectFlankSize(ExtendedBlock& block, size_t nTax, const PresenceChecker& presenceChecker,
		const std::string& T, const Options& options, bool directionRight) {
	if (!options.dynamicFlanks) {
		return findPerfectFlankSizeStatic(block, nTax, presenceChecker, T, options, directionRight);
	} else if (!options.useHMM) {
		return findPerfectFlankSizeDelta(block, nTax, presenceChecker, T, options, directionRight);
	} else if (!options.useBigGapsCriterion) {
		return findPerfectFlankSizeHMM(block, nTax, presenceChecker, T, options, directionRight);
	} else {
		return findPerfectFlankSizeBigGaps(block, nTax, presenceChecker, T, options, directionRight);
	}
}

ExtendedBlock extendBlock(const Seed& seededBlock, const std::string& T, size_t nTax, PresenceChecker& presenceChecker,
		const Options& options) {
	ExtendedBlock block(seededBlock, nTax, options.noIndels);

	if (options.maxMismatches == 0 && (!options.trimSeeds || !options.simpleTrimming)) {
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

	std::pair<size_t, double> bestLeft;
	std::pair<size_t, double> bestRight;

	if (!options.dynamicFlanks && options.fixedFlanks) {
		bestLeft.first = options.flankWidth;
		bestRight.first = options.flankWidth;
		// we still need to add the chars to the MSA
		std::vector<size_t> taxIDs = block.getTaxonIDsInBlock();
		for (size_t i = 1; i <= options.flankWidth; ++i) {
			std::vector<char> charsToAddLeft(taxIDs.size());
			std::vector<char> charsToAddRight(taxIDs.size());
			for (size_t j = 0; j < taxIDs.size(); ++j) {
				size_t coordRight = 0;
				size_t coordLeft = 0;
				coordRight = block.getTaxonCoordsWithoutFlanks(taxIDs[j]).second + i;
				coordLeft = block.getTaxonCoordsWithoutFlanks(taxIDs[j]).first - i;
				if (coordRight >= T.size() || coordLeft >= T.size()) {
					throw std::runtime_error("This should not happen! Coord is too large.");
				}
				charsToAddLeft[j] = T[coordLeft];
				charsToAddRight[j] = T[coordRight];
			}
			block.msaWrapper.addCharsRight(charsToAddRight);
			block.msaWrapper.addCharsLeft(charsToAddLeft);
		}
	} else {
		bestLeft = findPerfectFlankSize(block, nTax, presenceChecker, T, options, false);
		bestRight = findPerfectFlankSize(block, nTax, presenceChecker, T, options, true);
	}

	block.setLeftFlankSize(bestLeft.first);
	block.setRightFlankSize(bestRight.first);

	block.msaWrapper.shrinkDownToLeftFlank(bestLeft.first);
	block.msaWrapper.shrinkDownToRightFlank(bestRight.first);

	block.msaWrapper.clearMSADataStructures();

	return block;
}
