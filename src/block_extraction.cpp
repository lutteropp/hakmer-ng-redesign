/*
 * block_extraction.cpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#include "block_extraction.hpp"
#include "dna_functions.hpp"
#include "indexing/approx_matching.hpp"

#include <unordered_set>
#include <cmath>
#include <queue>
#include <algorithm>

#include "indexing/suffix_array_classic.hpp"

#include "build_superseeds.hpp"
#include "block_extension.hpp"
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

bool acceptSeed(size_t actSAPos, size_t matchCount, size_t k, size_t nTax, const std::vector<size_t>& SA, PresenceChecker& presenceChecker,
		const std::vector<IndexedTaxonCoords>& taxonCoords, const std::string& T, const Options& options) {
	if (matchCount > nTax /*|| matchCount > options.maxTaxaPerBlock*/) { // easy test for paralogy
		return false;
	}
	if (matchCount < options.minTaxaPerBlock) { // easy test for not enough taxa
		return false;
	}

	if (!options.lowComplexity) {
		std::string seed = T.substr(SA[actSAPos], k);
		if (!highComplexity(seed, k / 8, 1.0)) {
			//std::cout << "rejecting seed complexity: " << seed << "\n";
			return false;
		}
		//std::cout << "accepting seed complexity: " << seed << "\n";
	}

	size_t concatSize = T.size();

	// more complicated check for paralogy
	std::unordered_set<size_t> takenTaxa;
	for (size_t i = actSAPos; i < actSAPos + matchCount; ++i) {
		size_t taxID = posToTaxon(SA[i], taxonCoords, concatSize, options.reverseComplement);
		if (takenTaxa.find(taxID) != takenTaxa.end()) { // multiple match in taxon
			return false;
		}
		takenTaxa.insert(taxID);
	}

	return true;
}

size_t countMatches(size_t actSAPos, const std::vector<size_t>& lcp, size_t k) {
	size_t matchCount = 1;
	for (size_t i = actSAPos + 1; i < lcp.size(); ++i) {
		if (lcp[i] >= k) {
			matchCount++;
		} else {
			break;
		}
	}
	return matchCount;
}

/*void trimExtendedBlockSimple(ExtendedBlock& block, PresenceChecker& presenceChecker) {
 std::vector<size_t> taxIDs = block.getTaxonIDsInBlock();

 // first, trim the left side
 while (block.getAverageSeedSize() > 0) {
 bool ok = true;
 for (size_t i = 0; i < taxIDs.size(); ++i) {
 size_t coord = block.getTaxonCoordsWithFlanks(taxIDs[i]).first;
 if (!presenceChecker.isFree(coord)) {
 ok = false;
 break;
 }
 }
 if (ok) {
 break;
 } else {
 // trim back the left side by 1
 block.increaseAllTaxonCoordsLeft();
 }
 }
 // then, trim the right side (if seed.size() is not already zero)
 while (block.getAverageSeedSize() > 0) {
 bool ok = true;
 for (size_t i = 0; i < taxIDs.size(); ++i) {
 size_t coord = block.getSeedCoords(taxIDs[i]).second;
 if (!presenceChecker.isFree(coord)) {
 ok = false;
 break;
 }
 }
 if (ok) {
 break;
 } else {
 // trim back the right side by 1
 block.decreaseAllTaxonCoordsRight();
 }
 }
 for (size_t i = 0; i < taxIDs.size(); ++i) {
 if (!block.hasTaxon(taxIDs[i])) {
 block.removeTaxon(taxIDs[i]);
 }
 }
 }

 void trimExtendedBlock(ExtendedBlock& block, PresenceChecker& presenceChecker, Options& options) {
 if (options.simpleTrimming) {
 trimExtendedBlockSimple(block, presenceChecker);
 } else {
 std::vector<size_t> taxIDs = block.getTaxonIDsInBlock();
 for (size_t i = 0; i < taxIDs.size(); ++i) {
 size_t tID = taxIDs[i];
 while (block.hasTaxon(tID)) {
 size_t leftCoord = block.getTaxonCoordsWithFlanks(tID).first;
 if (!presenceChecker.isFree(leftCoord)) {
 block.increaseTaxonCoordLeft(tID);
 block.addGapLeft(tID);
 } else {
 break;
 }
 }
 while (block.hasTaxon(tID)) {
 size_t rightCoord = block.getTaxonCoordsWithFlanks(tID).second;
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
 }
 assert(block.getTaxonIDsInBlock().size() == block.getNTaxInBlock());
 }*/

void trimSeededBlock(Seed& block, PresenceChecker& presenceChecker, const Options& options) {
	std::vector<size_t> taxIDs = block.getTaxonIDsInBlock();
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

size_t findBestKForMatchcount(size_t sIdx, const std::vector<size_t>& lcp, size_t nMatches) {
	if (nMatches == 0 || nMatches == 1) {
		return std::numeric_limits<size_t>::max();
	}
	size_t kmax = lcp[sIdx + 1];
	for (size_t i = 2; i < nMatches; ++i) {
		kmax = std::min(kmax, lcp[sIdx + i]);
	}
	return kmax;
}

size_t findMinK(size_t sIdx, const std::vector<size_t>& lcp, size_t nTax) {
	if (sIdx + nTax >= lcp.size()) {
		return 1;
	}
	// find the smallest value for k such that not more than nTax matches occur
	return findBestKForMatchcount(sIdx, lcp, nTax + 1) + 1; // best k for at least nTax + 1 matches, +1 more
}

void processExtendedBlockBuffer(std::vector<ExtendedBlock>& extendedBlockBuffer, const Options& options, SummaryStatistics& stats,
		BlockWriter& writer, const IndexedConcatenatedSequence& concat) {
#pragma omp parallel for
	for (size_t i = 0; i < extendedBlockBuffer.size(); ++i) {
		ExtendedBlock extendedBlock = extendedBlockBuffer[i];
#pragma omp critical
		stats.updateSummaryStatistics(extendedBlock, concat.nTax());
		if (!options.outpath.empty()) {
			writer.writeTemporaryBlockMSA(extendedBlock, concat.getConcatenatedSeq(), concat.nTax(), options);
		}
	}
}

std::vector<ExtendedBlock> processSeedInfoBuffer(std::vector<SeedInfo>& seedInfoBuffer, const IndexedConcatenatedSequence& concat,
		PresenceChecker& presenceChecker, const Options& options, ApproximateMatcher& approxMatcher) {
	std::vector<ExtendedBlock> extendedBlocks;

#pragma omp parallel for
	for (size_t i = 0; i < seedInfoBuffer.size(); ++i) {
		SeedInfo seedInfo = seedInfoBuffer[i];

		// make a SeededBlock out of the seedInfo
		Seed block(concat.nTax());
		size_t sIdx = seedInfo.saPos;
		size_t k = seedInfo.k;
		size_t matchCount = seedInfo.n;
		for (size_t i = sIdx; i < sIdx + matchCount; ++i) {
			block.addTaxon(
					posToTaxon(concat.getSuffixArray()[i], concat.getTaxonCoords(), concat.getConcatenatedSeq().size(),
							options.reverseComplement), concat.getSuffixArray()[i], concat.getSuffixArray()[i] + k - 1);
		}
		for (std::pair<size_t, size_t> extraOcc : seedInfo.extraOccs) {
			block.addTaxon(
					posToTaxon(extraOcc.first, concat.getTaxonCoords(), concat.getConcatenatedSeq().size(), options.reverseComplement),
					extraOcc.first, extraOcc.second);
		}
		// first check if the current seed candidate is hopeless
		trimSeededBlock(block, presenceChecker, options);
		if (block.getNTaxInBlock() < options.minTaxaPerBlock) {
			continue;
		}

		// TODO: If we want to extend the block already here, then we need a trimExtendedBlock function.

		// Final trimming and adding, this time in critical mode
#pragma omp critical
		{
			// Extend the seed
			trivialExtension(block, concat.getConcatenatedSeq(), presenceChecker, concat.nTax(), options);
			ExtendedBlock extendedBlock = extendBlock(block, concat.getConcatenatedSeq(), concat.nTax(), presenceChecker, options);

			// TODO: Correctly implement extended block trimming, this is needed because of the parallelization
			//trimExtendedBlock(extendedBlock, presenceChecker, options);
			presenceChecker.reserveExtendedBlock(extendedBlock);
			extendedBlocks.push_back(extendedBlock);
		}
	}
	return extendedBlocks;
}

std::vector<SeedInfo> extractSeedInfos(const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker,
		const Options& options) {
	std::vector<SeedInfo> res;
	size_t actSAPos = 0;
	double lastP = 0;

	ApproximateMatcher approxMatcher(true);

#pragma omp parallel for schedule(dynamic)
	for (size_t sIdx = 0; sIdx < concat.getSuffixArray().size() - options.minTaxaPerBlock; ++sIdx) {
		bool ok = true;
		for (size_t i = 1; i < options.minTaxaPerBlock; ++i) {
			if (concat.getLcpArray()[sIdx + i] < options.minK) {
				ok = false;
				break;
			}
		}
		if (!ok) {
			continue;
		}

		size_t startPos = concat.getSuffixArray()[sIdx];
		size_t k = options.minK;
		if ((startPos + k >= concat.getConcatenatedSeq().size() || !presenceChecker.isFree(startPos, startPos + k - 1))) {
			continue;
		}
		size_t matchCount = countMatches(sIdx, concat.getLcpArray(), k);

		size_t actMinK = 0;
		if (matchCount >= options.minTaxaPerBlock) {
			actMinK = std::min(options.maxK, findMinK(sIdx, concat.getLcpArray(), concat.nTax()));
			if (k < actMinK) {
				k = actMinK;
			}
		}
		while (matchCount >= options.minTaxaPerBlock) {
			//we can stop if at least one of the exact matches we've found occurs before the current match in the suffix array... if we don't do that, we don't find all matches!
			bool stopEarly = false;
			if (concat.getLcpArray()[sIdx] >= k) {
				stopEarly = true;
			}

			if (!stopEarly
					&& acceptSeed(sIdx, matchCount, k, concat.nTax(), concat.getSuffixArray(), presenceChecker, concat.getTaxonCoords(),
							concat.getConcatenatedSeq(), options)) {

				SeedInfo info(sIdx, k, matchCount);
				// Augment the seed with mismatches
				if (options.maxMismatches > 0) {
					std::unordered_set<size_t> taxIDs;
					for (size_t t = 0; t < matchCount; ++t) {
						taxIDs.insert(
								posToTaxon(concat.getSuffixArray()[sIdx + t], concat.getTaxonCoords(), concat.getConcatenatedSeq().size(),
										options.reverseComplement));
					}

					std::string pattern = concat.getConcatenatedSeq().substr(concat.getSuffixArray()[sIdx], k);
					std::vector<std::pair<size_t, size_t> > extraOccs = approxMatcher.findOccurrences(concat.getConcatenatedSeq(),
							concat.getSuffixArray(), presenceChecker, pattern, options.maxMismatches, 1, false);
					for (size_t i = 0; i < extraOccs.size(); ++i) {
						size_t taxID = posToTaxon(extraOccs[i].first, concat.getTaxonCoords(), concat.getConcatenatedSeq().size(),
								options.reverseComplement);
						if (taxID < concat.nTax() && taxIDs.find(taxID) == taxIDs.end()) {
							// check if the extra occurrence is fine
							if (extraOccs[i].second >= concat.getConcatenatedSeq().size()) {
								continue;
							}
							info.addExtraOcc(extraOccs[i]);
							taxIDs.insert(taxID);
						}
					}
					assert(info.n + info.extraOccs.size() <= concat.nTax());
				}

#pragma omp critical
				res.push_back(info);
				if (options.verboseDebug) {
#pragma omp critical
					std::cout << "Pushing back a seeded occurrence with " << info.n << " taxa and seed size " << info.k << "\n";
				}
				break;
			} else {
				if (stopEarly) {
					k = concat.getLcpArray()[sIdx];
				}

				if (k >= options.maxK || startPos + k + 1 >= concat.getConcatenatedSeq().size() || !presenceChecker.isFree(startPos + k)) { // no further extension of seed length, or newly added character would be already taken anyway
					break;
				}
				k++;
				matchCount = countMatches(sIdx, concat.getLcpArray(), k);
				// TODO: Find the largest k that would fit to the current count of matches
				size_t bestK = findBestKForMatchcount(sIdx, concat.getLcpArray(), matchCount);
				if (bestK < k) {
					throw std::runtime_error("the cat jumps over the fox");
				}
				k = bestK;
				if (k >= options.maxK) {
					break;
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

void extractExtendedBlocks(const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker, BlockWriter& writer,
		SummaryStatistics& stats, const Options& options, size_t minK, size_t nMin, size_t nMax) {
	PresenceChecker seedingPresenceChecker(presenceChecker);
	std::cout << "Extracting seeded blocks...\n";
	std::vector<SeedInfo> seededBlockInfos = extractSeedInfos(concat, seedingPresenceChecker, options);
	std::cout << "seeded block infos.size(): " << seededBlockInfos.size() << "\n";
	std::cout << "Processing seeds...\n";
	std::sort(seededBlockInfos.begin(), seededBlockInfos.end(), std::greater<SeedInfo>());

	ApproximateMatcher approxMatcher(options.mismatchesOnly);
	std::vector<SeedInfo> buffer;
	size_t lastN = seededBlockInfos[0].n;
	buffer.push_back(seededBlockInfos[0]);
	for (size_t i = 1; i < seededBlockInfos.size(); ++i) {
		if (seededBlockInfos[i].n == lastN) {
			buffer.push_back(seededBlockInfos[i]);
		} else {
			std::cout << "Number of seeded blocks with " << lastN << " exact matches: " << buffer.size() << "\n";
			std::vector<ExtendedBlock> extendedBlockBuffer = processSeedInfoBuffer(buffer, concat, presenceChecker, options, approxMatcher);
			std::cout << "Assembling " << extendedBlockBuffer.size() << " extended blocks with " << lastN << " exact taxa...\n";
			processExtendedBlockBuffer(extendedBlockBuffer, options, stats, writer, concat);
			buffer.clear();
			lastN = seededBlockInfos[i].n;
			buffer.push_back(seededBlockInfos[i]);
		}
	}
	// process the last buffer
	std::cout << "Number of seeded blocks with " << lastN << " exact matches: " << buffer.size() << "\n";
	std::vector<ExtendedBlock> extendedBlockBuffer = processSeedInfoBuffer(buffer, concat, presenceChecker, options, approxMatcher);
	std::cout << "Assembling " << extendedBlockBuffer.size() << " extended blocks with " << lastN << " exact taxa...\n";
	processExtendedBlockBuffer(extendedBlockBuffer, options, stats, writer, concat);
	// TODO: Remove me again, this is just out of curiosity
	//std::vector<Superseed> superseeds = buildSuperseeds(seededBlocks, concat.getConcatenatedSeq(), presenceChecker, concat.nTax(), options);
	//std::cout << "We have made " << superseeds.size() << " superseeds out of " << seededBlocks.size() << " seeded blocks.\n";
}
