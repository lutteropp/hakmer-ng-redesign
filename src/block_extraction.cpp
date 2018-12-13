/*
 * block_extraction.cpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#include "block_extraction.hpp"

#include <stddef.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "alignment/simple_coords.hpp"
#include "alignment/simple_msa.hpp"
#include "block_extension.hpp"
#include "dna_functions.hpp"
#include "indexing/approx_matching.hpp"
#include "seed_info.hpp"

size_t posToTaxon(size_t pos, const std::vector<IndexedTaxonCoords>& taxonCoords, size_t concatSize, bool revComp) {
	if (revComp && pos >= concatSize / 2) {
		pos = concatSize - pos - 1;
	}
	for (size_t i = 0; i < taxonCoords.size(); ++i) {
		if (taxonCoords[i].contains(pos)) {
			return i;
		}
	}
	return std::string::npos;
}

/**
 * Compute the entropy of a given k-mer.
 * @param S The entire sequence data.
 * @param startPos Index of the first base in the k-mer.
 * @param k Size of the k-mer.
 * @param NCount Variable to store the number of 'N' bases in the k-mer.
 */
float info_content(const std::string& seed, unsigned int *NCount) {
	char histo[64] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	float x = 0.0;
	*NCount = 0;
	if (seed.size() == 3) {
		return 0.0; // degenerate case
	}

	size_t p = 0, q = 0;

	for (unsigned int i = 0; i <= seed.size() - 3; i++) {
		unsigned int I = 0;
		q = p;
		for (int j = 0; j < 3; j++) {
			I = (I << 2);
			switch (seed[q]) {
			case 'A':
				break;
			case 'C':
				I |= 1;
				break;
			case 'G':
				I |= 2;
				break;
			case 'T':
				I |= 3;
				break;
			case 'U': // Added by Sarah
				I |= 3;
				break;
			case 'N':
				++*NCount;
				break; // implicitly treating an N as an A...
			default:
				break; // implicitly treating all other chars as A also... :)
			}
			++q;
		}
		++p;
		++histo[I];
	}
	for (unsigned int i = 0; i < 64; i++) {
		x += histo[i] * (histo[i] - 1) / 2.0;
	}
	x /= (seed.size() - 3);
	return x;
}

/**
 * Checks whether a k-mer is of high enough complexity. For example, k-mers like AAAAAAAAAAA are bad and will fail the test.
 * @param S The input sequence data.
 * @param Position of the first base of the k-mer in the sequence data.
 * @param kmerSize Size of the k-mer.
 * @param Ncutoff Maximum number of 'N' bases allowed in a k-mer.
 * @param lowComplexityCutoff The k-mer will only be accepted if its complexity is strictly greater than this value.
 * @param gLowComplexityKmers Variable to keep track of the number of rejected k-mers due to low complexity.
 */
bool highComplexity(const std::string& seed, unsigned int Ncutoff, double lowComplexityCutoff) {
	unsigned int NCount;
	float lowComplexity = info_content(seed, &NCount);

	if (NCount > Ncutoff || lowComplexity > lowComplexityCutoff) {
		return false;
	} else {
		return true;
	}
}

bool acceptSeedComplexity(size_t actSAPos, size_t matchCount, size_t k, size_t nTax, const std::vector<size_t>& SA,
		PresenceChecker& presenceChecker, const std::vector<IndexedTaxonCoords>& taxonCoords, const std::string& T,
		const Options& options) {
	if (!options.lowComplexity) {
		std::string seed = T.substr(SA[actSAPos], k);
		if (!highComplexity(seed, k / 8, 1.0)) {
			return false;
		}
	}
	return true;
}

void trimSeededBlockExtra(Seed& block, PresenceChecker& presenceChecker, const Options& options) {
	// Further trimming: Prune seed nucleotide occurrences which are kept in too few taxa.
	std::vector<size_t> taxIDs = block.getTaxonIDsInBlock();
	if (taxIDs.size() < options.minTaxaPerBlock) {
		return;
	}
	std::vector<unsigned int> nFine(
			block.getSeedSize(taxIDs[0]) + block.getSeedCoords(taxIDs[0]).leftGapSize + block.getSeedCoords(taxIDs[0]).rightGapSize, 0);
	for (size_t i = 0; i < nFine.size(); ++i) {
		for (size_t tID : taxIDs) {
			if (block.getSeedCoords(tID).leftGapSize <= i && block.getSeedCoords(tID).leftGapSize + block.getSeedSize(tID) >= i) {
				nFine[i]++;
			}
		}
	}
	// left side...
	for (size_t i = 0; i < nFine.size(); ++i) {
		if (nFine[i] >= options.minSeedTaxInColumn) {
			break;
		}
		// now we must remove the leftmost characters
		for (size_t tID : taxIDs) {
			if (block.getSeedCoords(tID).leftGapSize <= i) {
				block.increaseTaxonCoordLeft(tID);
				block.addGapLeft(tID);
			}
		}
	}
	// right side...
	for (size_t i = 0; i < nFine.size(); ++i) {
		if (nFine[nFine.size() - i - 1] >= options.minSeedTaxInColumn) {
			break;
		}
		// now we must remove the rightmost characters
		for (size_t tID : taxIDs) {
			if (block.getSeedCoords(tID).leftGapSize + block.getSeedSize(tID) <= nFine.size() - i - 1) {
				block.decreaseTaxonCoordRight(tID);
				block.addGapRight(tID);
			}
		}
	}
	for (size_t tID : taxIDs) {
		if (!block.hasTaxon(tID)) {
			block.removeTaxon(tID);
		}
	}

	// Remove taxa with not enough kept seed sequence data from the seed
	for (size_t i = 0; i < taxIDs.size(); ++i) {
		if (block.hasTaxon(taxIDs[i]) && block.getSeedCoords(taxIDs[i]).size() < options.minSeedSitesKept) {
			block.removeTaxon(taxIDs[i]);
		}
	}

	// reduce gap sizes to reasonable amount
	// TODO: This, of course, changes the k value
	taxIDs = block.getTaxonIDsInBlock();

	if (taxIDs.empty()) {
		return;
	}

	size_t minLeftGapSize = std::numeric_limits<size_t>::max();
	size_t minRightGapSize = std::numeric_limits<size_t>::max();
	for (size_t tID : taxIDs) {
		minLeftGapSize = std::min(minLeftGapSize, block.getSeedCoords(tID).leftGapSize);
		minRightGapSize = std::min(minRightGapSize, block.getSeedCoords(tID).rightGapSize);
	}
	while (minLeftGapSize > 0) {
		for (size_t tID : taxIDs) {
			block.removeGapLeft(tID);
		}
		minLeftGapSize--;
	}
	while (minRightGapSize > 0) {
		for (size_t tID : taxIDs) {
			block.removeGapRight(tID);
		}
		minRightGapSize--;
	}
}

void trimSeededBlock(Seed& block, PresenceChecker& presenceChecker, const Options& options) {
	std::vector<size_t> taxIDs = block.getTaxonIDsInBlock();
	// prune positions which are overlapping with other, already taken blocks
	for (size_t i = 0; i < taxIDs.size(); ++i) {
		size_t tID = taxIDs[i];
		while (block.hasTaxon(tID)) {
			size_t leftCoord = block.getSeedCoords(tID).first;
			if (!presenceChecker.isFree(leftCoord)) {
				block.increaseTaxonCoordLeft(tID);
				block.addGapLeft(tID);
			} else {
				break;
			}
		}
		while (block.hasTaxon(tID)) {
			size_t rightCoord = block.getSeedCoords(tID).second;
			if (!presenceChecker.isFree(rightCoord)) {
				block.decreaseTaxonCoordRight(tID);
				block.addGapRight(tID);
			} else {
				break;
			}
		}
		if (!block.hasTaxon(tID)) {
			block.removeTaxon(tID);
		}
	}
	assert(block.getTaxonIDsInBlock().size() == block.getNTaxInBlock());
}

ExtendedBlock trimmedExtendedBlock(const SeedInfo& seedInfo, const std::vector<std::pair<size_t, size_t>>& extraOccs,
		const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker, const Options& options) {
	// make a new seed out of the extended block.
	std::vector<size_t> taxonUsage(concat.nTax(), 0);
	Seed seed(concat.nTax(), seedInfo);
	for (size_t i = 0; i < seedInfo.n; ++i) {
		size_t firstCoord = concat.getSuffixArray()[seedInfo.saPos + i];
		size_t taxID = concat.posToTaxon(firstCoord);
		size_t lastCoord = firstCoord + seedInfo.k - 1;
		seed.addTaxon(taxID, firstCoord, lastCoord);
		taxonUsage[taxID]++;
	}
	for (size_t i = 0; i < extraOccs.size(); ++i) {
		size_t taxID = concat.posToTaxon(extraOccs[i].first);
		if (!seed.hasTaxon(taxID)) {
			seed.addTaxon(taxID, extraOccs[i].first, extraOccs[i].second);
		}
		taxonUsage[taxID]++;
	}
	if (options.discardParalogMismatches) {
		for (size_t i = 0; i < concat.nTax(); ++i) {
			if (taxonUsage[i] > 1) {
				seed.removeTaxon(i);
			}
		}
	}
	trimSeededBlock(seed, presenceChecker, options);
	trimSeededBlockExtra(seed, presenceChecker, options);
	trivialExtensionPartial(seed, concat.getConcatenatedSeq(), presenceChecker, concat.nTax(), options);
	return extendBlock(seed, concat.getConcatenatedSeq(), concat.nTax(), presenceChecker, seedInfo.k, options);
}

double jukesCantorCorrection(double dist) {
	// TODO: Jukes Cantor Correction doesn't work if dist >= 0.75. In this case, it will return infinity.
	return -0.75 * log(1 - (4.0 / 3) * dist);
}

double averageSubstitutionRate(const std::vector<std::string>& msa) {
	std::vector<size_t> validSeqIndices;
	for (size_t i = 0; i < msa.size(); ++i) {
		for (size_t j = 0; j < msa[i].size(); ++j) {
			if (msa[i][j] != '?' && msa[i][j] != '-') {
				validSeqIndices.push_back(i);
				break;
			}
		}
	}
	//size_t nPairs = msa[0].size() * (validSeqIndices.size() * (validSeqIndices.size() - 1) / 2);
	size_t nPairs = 0;
	size_t nDiff = 0;
	for (size_t i = 0; i < msa[0].size(); ++i) {
		for (size_t j = 0; j < validSeqIndices.size(); ++j) {
			for (size_t k = j + 1; k < validSeqIndices.size(); ++k) {
				/*if (msa[validSeqIndices[j]][i] != '-' && msa[validSeqIndices[k]][i] != '-' && msa[validSeqIndices[j]][i] != '?'
				 && msa[validSeqIndices[k]][i] != '?' && !ambiguousMatch(msa[validSeqIndices[j]][i], msa[validSeqIndices[k]][i])) {
				 nDiff++;
				 }*/
				if (msa[validSeqIndices[j]][i] != '-' && msa[validSeqIndices[k]][i] != '-' && msa[validSeqIndices[j]][i] != '?'
						&& msa[validSeqIndices[k]][i] != '?') {
					nPairs++;
					if (!ambiguousMatch(msa[validSeqIndices[j]][i], msa[validSeqIndices[k]][i])) {
						nDiff++;
					}
				}
			}
		}
	}
	double subRate = (double) nDiff / nPairs;
	return jukesCantorCorrection(subRate);
}

void addValidExtraOccs(const std::vector<std::pair<size_t, size_t> >& extraOccs, const IndexedConcatenatedSequence& concat,
		const Options& options, Seed& newSeed, const ExtendedBlock& block) {
	std::unordered_set<size_t> taxIDs;
	for (size_t tID : block.getTaxonIDsInBlock()) {
		taxIDs.insert(tID);
	}
	std::vector<size_t> taxUsage(concat.nTax(), 0);
	for (size_t i = 0; i < extraOccs.size(); ++i) {
		size_t taxID = posToTaxon(extraOccs[i].first, concat.getTaxonCoords(), concat.getConcatenatedSeq().size(),
				options.reverseComplement);
		taxUsage[taxID]++;
		if (taxID < concat.nTax() && taxIDs.find(taxID) == taxIDs.end()) {
			// if we are here, the occurrence is fine. This means work for us.
			newSeed.addTaxon(taxID, extraOccs[i].first, extraOccs[i].second);
			taxIDs.insert(taxID);
		}
	}
	if (options.discardParalogMismatches) { // remove paralog mismatches again
		for (size_t i = 0; i < concat.nTax(); ++i) {
			if (taxUsage[i] > 1) {
				newSeed.removeTaxon(i);
			}
		}
	}
}

size_t processExtendedBlockBuffer(std::vector<ExtendedBlock>& extendedBlockBuffer, const Options& options, SummaryStatistics& stats,
		BlockWriter& writer, const IndexedConcatenatedSequence& concat, ApproximateMatcher& approxMatcher,
		PresenceChecker& presenceChecker) {
	size_t newlyAddedBlocks = 0;
#pragma omp parallel for reduction(+:newlyAddedBlocks)
	for (size_t i = 0; i < extendedBlockBuffer.size(); ++i) {
		ExtendedBlock block = extendedBlockBuffer[i];
		std::vector<std::pair<size_t, size_t> > extraOccs;

		if (block.getNTaxInBlock() < concat.nTax()) {
			double subRate = block.getSubRate();
//#pragma omp critical
			//		std::cout << "subRate: " << subRate << "\n";
			size_t k = block.getMySeededBlock().mySeedInfo.k;
			size_t maxMismatches = k * subRate;

			if (maxMismatches > 0) { // augment the block with approximate matches
				std::string pattern = concat.getConcatenatedSeq().substr(concat.getSuffixArray()[block.getMySeededBlock().mySeedInfo.saPos],
						k);
				std::vector<size_t> taxPresence(concat.nTax(), 0);
				for (size_t tID : block.getTaxonIDsInBlock()) {
					taxPresence[tID] += 2; // in order not to search for approx matches in the taxa we already have with exact matches
				}
				// do it on a per-taxon basis

				for (size_t taxonID = 0; taxonID < concat.nTax(); taxonID++) {
					if (taxPresence[taxonID] >= 2) {
						continue;
					}
					std::vector<std::pair<size_t, size_t> > eOccs = approxMatcher.findFewOccurrences(concat.getConcatenatedSeq(),
							concat.getSuffixArray(taxonID), concat.getLcpArray(taxonID), presenceChecker, pattern, maxMismatches, 1, false,
							concat, taxPresence);
					for (size_t eIdx = 0; eIdx < eOccs.size(); ++eIdx) {
						extraOccs.push_back(eOccs[eIdx]);
					}
				}
			}
		}

		if (block.getNTaxInBlock() >= options.minTaxaPerBlock) {
#pragma omp critical
			{
				presenceChecker.freeExtendedBlock(block);
				block = trimmedExtendedBlock(block.getMySeededBlock().mySeedInfo, extraOccs, concat, presenceChecker, options);
				if (block.getNTaxInBlock() >= options.minTaxaPerBlock) {
					presenceChecker.reserveExtendedBlock(block);
				}
			}

			if (block.getNTaxInBlock() >= options.minTaxaPerBlock) {
				std::vector<std::string> msa = computeMSA(block, concat.getConcatenatedSeq(), concat.nTax(), options);
#pragma omp critical
				{
					newlyAddedBlocks++;
					stats.updateSummaryStatistics(block, concat.nTax());
					if (!options.outpath.empty()) {
						writer.writeTemporaryBlockMSA(msa, concat.nTax());
					}
				}
			}
		}
	}
	return newlyAddedBlocks;
}

std::vector<ExtendedBlock> processSeedInfoBuffer(std::vector<SeedInfo>& seedInfoBuffer, const IndexedConcatenatedSequence& concat,
		PresenceChecker& presenceChecker, const Options& options, ApproximateMatcher& approxMatcher) {
	std::vector<ExtendedBlock> extendedBlocks;

#pragma omp parallel for
	for (size_t i = 0; i < seedInfoBuffer.size(); ++i) {
		SeedInfo seedInfo = seedInfoBuffer[i];
		// make a SeededBlock out of the seedInfo
		Seed block(concat.nTax(), seedInfo);
		size_t sIdx = seedInfo.saPos;
		size_t k = seedInfo.k;
		size_t matchCount = seedInfo.n;
		for (size_t i = sIdx; i < sIdx + matchCount; ++i) {
			size_t firstPos = concat.getSuffixArray()[i];
			size_t tID = concat.posToTaxon(firstPos);
			block.addTaxon(tID, firstPos, firstPos + k - 1);
		}

		// first check if the current seed candidate is hopeless
		trimSeededBlock(block, presenceChecker, options);
		if (block.getNTaxInBlock() < options.minTaxaPerBlock) {
			continue;
		}

		k = block.getSeedCoords(block.getTaxonIDsInBlock()[0]).size() + block.getSeedCoords(block.getTaxonIDsInBlock()[0]).leftGapSize
				+ block.getSeedCoords(block.getTaxonIDsInBlock()[0]).rightGapSize;

		// check some little assertion...
		for (size_t tID : block.getTaxonIDsInBlock()) {
			assert(k == block.getSeedCoords(tID).size() + block.getSeedCoords(tID).leftGapSize + block.getSeedCoords(tID).rightGapSize);
		}

		// Final trimming and adding
		size_t oldSeedSize = block.getAverageSeedSize();
		trimSeededBlock(block, presenceChecker, options);
		trimSeededBlockExtra(block, presenceChecker, options);

		if (block.getNTaxInBlock() >= options.minTaxaPerBlock) { // we need this extra check because in parallel mode, our taxa in the block could get invalidated all the time
			// Partially extend the seed
			trivialExtensionPartial(block, concat.getConcatenatedSeq(), presenceChecker, concat.nTax(), options);
			ExtendedBlock extendedBlock = extendBlock(block, concat.getConcatenatedSeq(), concat.nTax(), presenceChecker, oldSeedSize,
					options);

			bool discardMe = (options.discardUninformativeBlocks && extendedBlock.getAverageLeftFlankSize() == 0
					&& extendedBlock.getAverageRightFlankSize() == 0);
			if (!discardMe) {
				std::vector<std::string> msa = computeMSA(extendedBlock, concat.getConcatenatedSeq(), concat.nTax(), options);
				if (block.getNTaxInBlock() < concat.nTax()) {
					double subRate = averageSubstitutionRate(msa);
					if (subRate > options.maxAvgSubstitutionRate) {
						discardMe = true;
					} else {
						extendedBlock.setSubRate(subRate);
					}
				}
			}

			if (!discardMe) {
#pragma omp critical
				{
					presenceChecker.reserveExtendedBlock(extendedBlock);
					extendedBlocks.push_back(extendedBlock);
				}
			}
		}
	}
	return extendedBlocks;
}

std::pair<size_t, size_t> findKAndNumExactMatches(size_t saPos, const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker,
		size_t minK, const Options& options) {
	std::pair<size_t, size_t> emptyPair = std::make_pair(0, 0);
	size_t kStart = std::max(minK, concat.getLcpArray()[saPos] + 1);
	if ((concat.getSuffixArray()[saPos] + kStart >= concat.getConcatenatedSeq().size())
			|| (!presenceChecker.isFree(concat.getSuffixArray()[saPos], concat.getSuffixArray()[saPos] + kStart - 1))) {
		return emptyPair;
	}

// at least nMin next entries have to be included
	std::unordered_set<size_t> taxIDs;
	size_t minNeighborK = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < options.minTaxaPerBlock; ++i) {
		if (i > 0 && concat.getLcpArray()[saPos + i] < kStart) {
			return emptyPair;
		}
		size_t tID = concat.posToTaxon(concat.getSuffixArray()[saPos + i]);
		taxIDs.insert(tID);
		if (i > 0) {
			minNeighborK = std::min(minNeighborK, concat.getLcpArray()[saPos + i]);
		}
	}
	if (taxIDs.size() < options.minTaxaPerBlock) {
		return emptyPair;
	}
	size_t actMatches = options.minTaxaPerBlock;
// now, we are at a position where at least the minimum amount of taxa has been found. Let's check how much more we can go to the right without encountering any paralogs.
	size_t newMinK = kStart;
	for (size_t i = options.minTaxaPerBlock; i < concat.getSuffixArray().size() - saPos; ++i) {
		if (concat.getLcpArray()[saPos + i] < kStart) {
			break;
		}
		size_t tID = concat.posToTaxon(concat.getSuffixArray()[saPos + i]);
		taxIDs.insert(tID);
		if (taxIDs.size() != actMatches + 1) {
			newMinK = concat.getLcpArray()[saPos + i] + 1;
			break;
		}
		actMatches++;
	}

	if (minNeighborK < newMinK) {
		return emptyPair;
	}
	size_t k = newMinK; // now, we just need to check how many matches we still have and what's their min value

	size_t minValNewMatches = std::numeric_limits<size_t>::max();
	size_t nMatches = 1;
	for (size_t i = 1; i < concat.getSuffixArray().size() - saPos; ++i) {
		if (concat.getLcpArray()[saPos + i] < k) {
			break;
		}
		nMatches++;
		minValNewMatches = std::min(minValNewMatches, concat.getLcpArray()[saPos + i]);
	}
	k = minValNewMatches;

	return std::make_pair(k, nMatches);
}

std::vector<SeedInfo> extractSeedInfos(const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker, size_t minK,
		size_t maxK, const Options& options) {
	std::vector<SeedInfo> res;
	size_t actSAPos = 0;
	double lastP = 0;

	ApproximateMatcher approxMatcher(true);

#pragma omp parallel for schedule(dynamic)
	for (size_t sIdx = 0; sIdx < concat.getSuffixArray().size() - options.minTaxaPerBlock; ++sIdx) {
		bool skip = false;
		std::pair<size_t, size_t> p = findKAndNumExactMatches(sIdx, concat, presenceChecker, minK, options);
		size_t k = p.first;
		size_t matchCount = p.second;
		if (k == 0) {
			skip = true;
		}

		if (!skip
				&& acceptSeedComplexity(sIdx, matchCount, k, concat.nTax(), concat.getSuffixArray(), presenceChecker,
						concat.getTaxonCoords(), concat.getConcatenatedSeq(), options)) {

			// perform simple extension of the seed
			Seed block(concat.nTax());
			for (size_t idx = 0; idx < matchCount; ++idx) {
				size_t firstPos = concat.getSuffixArray()[sIdx + idx];
				size_t tID = concat.posToTaxon(firstPos);
				block.addTaxon(tID, firstPos, firstPos + k - 1);
			}

			if (!canGoLeftAll(block, presenceChecker, concat.nTax())
					|| !allLeftSame(block, concat.getConcatenatedSeq(), block.getTaxonIDsInBlock())) {
				SeedInfo info(sIdx, k, matchCount);
#pragma omp critical
				res.push_back(info);
				if (options.verboseDebug) {
#pragma omp critical
					std::cout << "Pushing back a seeded occurrence with " << info.n << " taxa and seed size " << info.k << "\n";
				}
			}
		}
		double progress = (double) 100 * sIdx / concat.getSuffixArray().size(); // TODO: Fix this, this looks kinda wrong in parallel mode
		if (progress > lastP + 1) {
#pragma omp critical
			{
				if (progress > lastP + 1) {
					std::cout << progress << " %\n";
					lastP = progress;
				}
			}
		}
	}
	return res;
}

size_t elbowMethod(const std::vector<std::pair<size_t, size_t> >& seedSizeCounts, const Options& options) {
// see https://www.linkedin.com/pulse/finding-optimal-number-clusters-k-means-through-elbow-asanka-perera/
// see https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
// we need to find the point with the largest distance to the line from the first to the last point; this point corresponds to our chosen minK value.
	size_t minReasonableCount = 1000;
	size_t lastIdx = seedSizeCounts.size() - 1;
	while (seedSizeCounts[lastIdx].second < minReasonableCount && lastIdx > 0) {
		lastIdx--;
	}
	int maxDist = 0;
	size_t maxDistIdx = 0;
	int x1 = seedSizeCounts[0].first;
	int y1 = seedSizeCounts[0].second;
	int x2 = seedSizeCounts[lastIdx].first;
	int y2 = seedSizeCounts[lastIdx].second;
	for (size_t i = 1; i <= lastIdx; ++i) { // because the endpoints trivially have distance 0
		int x0 = seedSizeCounts[i].first;
		int y0 = seedSizeCounts[i].second;
		int d = std::abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1);
		if (d > maxDist) {
			maxDist = d;
			maxDistIdx = i;
		}
	}
	return seedSizeCounts[maxDistIdx].first;
}

void printSeedSizeHistogram(const std::vector<std::pair<size_t, size_t> >& seedSizes) {
	std::cout << "Seed size histogram:\n";
	for (size_t i = 0; i < seedSizes.size(); ++i) {
		std::cout << seedSizes[i].first << ": " << seedSizes[i].second << "\n";
	}
}

std::vector<std::pair<size_t, size_t> > countSeedSizes(const std::vector<SeedInfo>& seedInfos, const Options& options) {
	std::vector<std::pair<size_t, size_t> > res;
	std::vector<size_t> seedSizes(1000, 0);
	for (size_t i = 0; i < seedInfos.size(); ++i) {
		if (seedInfos[i].k >= seedSizes.size()) {
			seedSizes.resize(seedInfos[i].k + 1);
		}
		seedSizes[seedInfos[i].k]++;
	}
	for (size_t i = 0; i < seedSizes.size(); ++i) {
		if (seedSizes[i] > 0) {
			std::pair<size_t, size_t> p = std::make_pair(i, seedSizes[i]);
			res.push_back(p);
		}
	}
	return res;
}

size_t selectAndProcessSeedInfos(const std::vector<SeedInfo>& seededBlockInfos, ApproximateMatcher& approxMatcher,
		const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker, BlockWriter& writer, SummaryStatistics& stats,
		const Options& options, size_t minK, size_t maxK) {
	size_t newlyAddedBlocks = 0;
	std::vector<SeedInfo> buffer;
	size_t lastN = seededBlockInfos[0].n;
	buffer.push_back(seededBlockInfos[0]);
	for (size_t i = 1; i < seededBlockInfos.size(); ++i) {
		if (seededBlockInfos[i].n == lastN) {
			if (seededBlockInfos[i].k >= minK && seededBlockInfos[i].k <= maxK) {
				buffer.push_back(seededBlockInfos[i]);
			}
		} else {
			std::cout << "Number of seeded blocks with " << lastN << " exact matches: " << buffer.size() << "\n";
			std::vector<ExtendedBlock> extendedBlockBuffer = processSeedInfoBuffer(buffer, concat, presenceChecker, options, approxMatcher);
			std::cout << "Assembling " << extendedBlockBuffer.size() << " extended blocks with " << lastN << " exact taxa...\n";
			newlyAddedBlocks += processExtendedBlockBuffer(extendedBlockBuffer, options, stats, writer, concat, approxMatcher,
					presenceChecker);
			buffer.clear();
			lastN = seededBlockInfos[i].n;
			if (seededBlockInfos[i].k >= minK && seededBlockInfos[i].k <= maxK) {
				buffer.push_back(seededBlockInfos[i]);
			}
		}
	}
// process the last buffer
	std::cout << "Number of seeded blocks with " << lastN << " exact matches: " << buffer.size() << "\n";
	std::vector<ExtendedBlock> extendedBlockBuffer = processSeedInfoBuffer(buffer, concat, presenceChecker, options, approxMatcher);
	std::cout << "Assembling " << extendedBlockBuffer.size() << " extended blocks with " << lastN << " exact taxa...\n";
	newlyAddedBlocks += processExtendedBlockBuffer(extendedBlockBuffer, options, stats, writer, concat, approxMatcher, presenceChecker);
	std::cout << "Newly added blocks: " << newlyAddedBlocks << "\n";
	return newlyAddedBlocks;
}

void printHypotheticalBestCaseTaxonCoverage(const IndexedConcatenatedSequence& concat, const std::vector<SeedInfo>& seedInfos) {
// if we would ignore overlaps and stuff like this
	std::vector<size_t> taxUsage(concat.nTax(), 0);
	for (size_t i = 0; i < seedInfos.size(); ++i) {
		for (size_t j = 0; j < seedInfos[i].n; ++j) {
			size_t tID = concat.posToTaxon(concat.getSuffixArray()[seedInfos[i].saPos + j]);
			taxUsage[tID] += seedInfos[i].k;
		}
	}

	std::cout << "\nHypothetical best-case taxon usages if we wouldn't care about overlaps and stuff like this:\n";
	for (size_t i = 0; i < concat.nTax(); ++i) {
		double p = ((double) taxUsage[i] * 100) / concat.getTaxonCoords(i).getTotalLength();
		std::cout << concat.getTaxonLabels()[i] << ": " << p << " %\n";
	}
	std::cout << "\n";
}

void extractExtendedBlocks(const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker, BlockWriter& writer,
		SummaryStatistics& stats, const Options& options, size_t minK, size_t maxK) {
	PresenceChecker seedingPresenceChecker(presenceChecker);
	size_t initialMinK = minK;

	std::cout << "Extracting seeded blocks...\n";
	std::vector<SeedInfo> seededBlockInfos = extractSeedInfos(concat, seedingPresenceChecker, minK, maxK, options);
	std::cout << "seeded block infos.size(): " << seededBlockInfos.size() << "\n";
	std::cout << "Processing seeds...\n";
	std::sort(seededBlockInfos.begin(), seededBlockInfos.end(), std::greater<SeedInfo>());

	std::vector<std::pair<size_t, size_t> > seedSizes;
	seedSizes = countSeedSizes(seededBlockInfos, options);
	printSeedSizeHistogram(seedSizes);
	printHypotheticalBestCaseTaxonCoverage(concat, seededBlockInfos);

	size_t newMinK = elbowMethod(seedSizes, options);
	std::cout << "New chosen minK by using the elbow method: " << newMinK << ". Ignoring all seeds with smaller k than this value.\n";
	minK = newMinK;
	ApproximateMatcher approxMatcher(options.mismatchesOnly);

	size_t newlyAddedBlocks = selectAndProcessSeedInfos(seededBlockInfos, approxMatcher, concat, presenceChecker, writer, stats, options,
			minK, maxK);
	double seqDataUsed = stats.getAmountSeqDataUsed(concat.getSequenceDataSize());
	while (seqDataUsed < options.minSeqDataUsage) {
		std::cout << "Current percentage of sequence data used: " << seqDataUsed * 100
				<< " %. This is too low. Trying to find more blocks with lower minimum kmer size.\n";
		size_t newMinK = std::max((size_t) initialMinK, (size_t) minK - 1);
		if (newMinK != minK) {
			std::cout << "Using new value for minK: " << newMinK << "\n";
			selectAndProcessSeedInfos(seededBlockInfos, approxMatcher, concat, presenceChecker, writer, stats, options, newMinK, minK - 1);
			minK = newMinK;
		} else {
			std::cout << "Unfortunately, a smaller minK value makes not much sense. :-(\n";
			break;
		}
		seqDataUsed = stats.getAmountSeqDataUsed(concat.getSequenceDataSize());
	}
	std::cout << "seqDataUsed: " << seqDataUsed << "\n";

// TODO: Remove me again, this is just out of curiosity
//std::vector<Superseed> superseeds = buildSuperseeds(seededBlocks, concat.getConcatenatedSeq(), presenceChecker, concat.nTax(), options);
//std::cout << "We have made " << superseeds.size() << " superseeds out of " << seededBlocks.size() << " seeded blocks.\n";
}
