/*
 * block_extraction.cpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#include "block_extraction.hpp"
#include "block_helper_functions.hpp"
#include "dna_functions.hpp"
#include "indexing/approx_matching.hpp"

#include <unordered_set>
#include <cmath>
#include <queue>
#include <algorithm>

#include "indexing/suffix_array_classic.hpp"

#include "build_superseeds.hpp"
#include "block_extension.hpp"

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

bool canStay(size_t pos, const std::vector<IndexedTaxonCoords>& taxonCoords, const std::vector<size_t>& wantedTaxa) {
	for (size_t i = 0; i < wantedTaxa.size(); ++i) {
		if (taxonCoords[wantedTaxa[i]].contains(pos)) {
			return true;
		}
	}
	return false;
}

// Shrink arrays... needed for quartet-based stuff.
std::pair<std::vector<size_t>, std::vector<size_t> > shrinkArrays(const IndexedConcatenatedSequence& concat,
		const std::vector<IndexedTaxonCoords>& taxonCoords, const std::vector<size_t>& wantedTaxa, const Options& options) {
	std::vector<size_t> resSA;
	std::vector<size_t> resLCP;

	bool recomputeNeeded = false;

	for (size_t i = 0; i < concat.getSuffixArray().size(); ++i) {
		if (canStay(concat.getSuffixArray()[i], taxonCoords, wantedTaxa)) {
			resSA.push_back(concat.getSuffixArray()[i]);
			size_t lcpVal = concat.getLcpArray()[i];
			if (recomputeNeeded) {
				// recompute lcpVal
				size_t lTop = std::min(concat.getConcatSize(), options.maxK);
				lcpVal = longestCommonPrefix(concat.getConcatenatedSeq(), resSA[resSA.size() - 2], resSA[resSA.size() - 1], lTop);
				recomputeNeeded = false;
			}
			resLCP.push_back(lcpVal);
		} else {
			recomputeNeeded = true;
		}
	}
	return std::make_pair(resSA, resLCP);
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

bool acceptSeed(size_t actSAPos, size_t matchCount, const std::vector<std::pair<size_t, size_t> >& extraOccs, size_t k, size_t nTax,
		const std::vector<size_t>& SA, PresenceChecker& presenceChecker, const std::vector<IndexedTaxonCoords>& taxonCoords,
		const std::string& T, const Options& options, size_t maxMatches) {
	if (matchCount + extraOccs.size() > nTax || matchCount + extraOccs.size() > maxMatches /*|| matchCount > options.maxTaxaPerBlock*/) { // easy test for paralogy
		return false;
	}
	if (matchCount + extraOccs.size() < options.minTaxaPerBlock) { // easy test for not enough taxa
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
	for (size_t i = 0; i < extraOccs.size(); ++i) {
		size_t taxID = posToTaxon(extraOccs[i].first, taxonCoords, concatSize, options.reverseComplement);
		if (takenTaxa.find(taxID) != takenTaxa.end()) { // multiple match in taxon
			return false;
		}
		takenTaxa.insert(taxID);
	}

	size_t flankOffset = 0;
	if (!options.dynamicFlanks && options.fixedFlanks) {
		flankOffset = options.flankWidth;
	}

	// check for presence/absence of the whole region, including flanks
	for (size_t i = actSAPos; i < actSAPos + matchCount; ++i) {
		if (flankOffset > SA[i] || SA[i] + k - 1 + flankOffset >= concatSize) {
			return false;
		}
		if (!presenceChecker.isFree(SA[i] - flankOffset, SA[i] + k - 1 + flankOffset)) {
			return false;
		}
	}

	for (size_t i = 0; i < extraOccs.size(); ++i) {
		if (flankOffset > extraOccs[i].first || extraOccs[i].second + flankOffset >= concatSize) {
			return false;
		}
		if (!presenceChecker.isFree(extraOccs[i].first - flankOffset, extraOccs[i].second + flankOffset)) {
			return false;
		}
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

size_t findLargestK(const std::vector<size_t>& lcp, const std::vector<size_t>& sa, const std::vector<IndexedTaxonCoords>& taxonCoords,
		size_t concatSize, bool rc, const PresenceChecker& presenceChecker) {
	size_t max = 0;
	for (size_t i = 0; i < lcp.size() - 4; ++i) {
		size_t act = std::min(std::min(std::min(lcp[i], lcp[i + 1]), lcp[i + 2]), lcp[i + 3]);
		if (act < max)
			continue;
		size_t nMatches = countMatches(i, lcp, act);

		// check for presence/absence of the whole region, including flanks
		for (size_t j = i; j < i + nMatches; ++j) {
			if (sa[j] + act - 1 + 0 >= concatSize) {
				continue;
			}
			if (!presenceChecker.isFree(sa[j], sa[j] + act - 1)) {
				continue;
			}
		}

		std::unordered_set<size_t> taxa;
		for (size_t j = 0; j < nMatches; ++j) {
			taxa.insert(posToTaxon(sa[i + j], taxonCoords, concatSize, rc));
		}
		if (taxa.size() == nMatches) {
			max = act;
		}
	}
	return std::min((size_t) 150, max);
}

// TODO: Re-add mismatches and indels in seeds
std::vector<Seed> extractSeededBlocksMismatchAugmentationOnly(const std::string& T, size_t nTax, const std::vector<size_t>& SA,
		const std::vector<size_t>& lcp, PresenceChecker& presenceChecker, const std::vector<IndexedTaxonCoords>& taxonCoords,
		const Options& options, size_t minK, size_t minMatches, size_t maxMatches) {
	std::vector<Seed> res;
	size_t actSAPos = 0;
	double lastP = 0;

	ApproximateMatcher approxMatcher(true);

#pragma omp parallel for schedule(dynamic)
	for (size_t sIdx = 0; sIdx < SA.size(); ++sIdx) {
		size_t startPos = SA[sIdx];
		size_t k = minK;
		if ((startPos + k >= T.size() || !presenceChecker.isFree(startPos, startPos + k - 1))) {
			continue;
		}
		size_t matchCount = countMatches(sIdx, lcp, k);
		std::vector<std::pair<size_t, size_t> > extraOccs;

		while (matchCount >= std::max(options.minTaxaPerBlock, minMatches)) {
			//we can stop if at least one of the exact matches we've found occurs before the current match in the suffix array... if we don't do that, we don't find all matches!
			bool stopEarly = false;
			if (lcp[sIdx] >= k) {
				stopEarly = true;
			}

			if (!stopEarly && acceptSeed(sIdx, matchCount, extraOccs, k, nTax, SA, presenceChecker, taxonCoords, T, options, maxMatches)) {
				if (options.largeSeeds) {
					// try to find largest seed size that still gets accepted
					size_t bestK = k;
					while (acceptSeed(sIdx, matchCount, extraOccs, k, nTax, SA, presenceChecker, taxonCoords, T, options, maxMatches)) {
						bestK = k;
						if (k == options.maxK || startPos + k + 1 >= T.size() || !presenceChecker.isFree(startPos + k)) { // no further extension of seed length, or newly added character would be already taken anyway
							break;
						}
						k++;
						matchCount = countMatches(sIdx, lcp, k);
					}
					k = bestK;
					matchCount = countMatches(sIdx, lcp, k);
				}

				Seed block(nTax);
				for (size_t i = sIdx; i < sIdx + matchCount; ++i) {
					block.addTaxon(posToTaxon(SA[i], taxonCoords, T.size(), options.reverseComplement), SA[i], SA[i] + k - 1);
				}

				if (options.maxMismatches > 0) {
					std::string pattern = T.substr(startPos, k);
					extraOccs = approxMatcher.findOccurrences(T, SA, presenceChecker, pattern, options.maxMismatches, 1, false);

					// TODO: Maybe only add those approximate matches that don't collide with the exact matches we already have?
					for (size_t i = 0; i < extraOccs.size(); ++i) {
						size_t taxID = posToTaxon(extraOccs[i].first, taxonCoords, T.size(), options.reverseComplement);
						if (!block.hasTaxon(taxID)) {
							size_t flankOffset = 0;
							if (!options.dynamicFlanks && options.fixedFlanks) {
								flankOffset = options.flankWidth;
							}
							// check if the extra occurrence is fine
							if (flankOffset > extraOccs[i].first || extraOccs[i].second + flankOffset >= T.size()) {
								continue;
							}
							if (!presenceChecker.isFree(extraOccs[i].first - flankOffset, extraOccs[i].second + flankOffset)) {
								continue;
							}
							block.addTaxon(posToTaxon(extraOccs[i].first, taxonCoords, T.size(), options.reverseComplement),
									extraOccs[i].first, extraOccs[i].second);
						}
					}
				}
				computeBestCaseMaxSizes(block, T, presenceChecker, nTax);
#pragma omp critical
				{
					res.push_back(block);
					if (options.verboseDebug) {
						std::cout << "Pushing back a seeded block with " << block.getNTaxInBlock() << " taxa and seed size "
								<< block.getSeedSize() << "\n";
					}
				}
				if (!options.quartetFlavor) {
					double progress = (double) 100 * sIdx / SA.size(); // TODO: Fix this, this looks kinda wrong in parallel mode
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
				break;
			} else {
				if (k == options.maxK || startPos + k + 1 >= T.size() || !presenceChecker.isFree(startPos + k)) { // no further extension of seed length, or newly added character would be already taken anyway
					break;
				}
				k++;
				matchCount = countMatches(sIdx, lcp, k);
			}
		}
	}
	return res;
}

// TODO: Re-add mismatches and indels in seeds
std::vector<Seed> extractSeededBlocks(const std::string& T, size_t nTax, const std::vector<size_t>& SA, const std::vector<size_t>& lcp,
		PresenceChecker& presenceChecker, const std::vector<IndexedTaxonCoords>& taxonCoords, const Options& options, size_t minK, size_t minMatches,
		size_t maxMatches) {
	if (options.mismatchAugmentationOnly) {
		return extractSeededBlocksMismatchAugmentationOnly(T, nTax, SA, lcp, presenceChecker, taxonCoords, options, minK, minMatches, maxMatches);
	}

	std::vector<Seed> res;
	size_t actSAPos = 0;
	double lastP = 0;

	ApproximateMatcher approxMatcher(true);

#pragma omp parallel for schedule(dynamic)
	for (size_t sIdx = 0; sIdx < SA.size(); ++sIdx) {
		size_t startPos = SA[sIdx];
		size_t k = options.minK;
		if ((startPos + k >= T.size() || !presenceChecker.isFree(startPos, startPos + k - 1))) {
			continue;
		}
		size_t matchCount = countMatches(sIdx, lcp, k);
		std::vector<std::pair<size_t, size_t> > extraOccs;

		while (matchCount + extraOccs.size() >= std::max(options.minTaxaPerBlock, minMatches)) {
			//we can stop if at least one of the exact matches we've found occurs before the current match in the suffix array... if we don't do that, we don't find all matches!
			bool stopEarly = false;
			if (lcp[sIdx] >= k) {
				stopEarly = true;
			}

			if (!stopEarly && acceptSeed(sIdx, matchCount, extraOccs, k, nTax, SA, presenceChecker, taxonCoords, T, options, maxMatches)) {
				if (options.largeSeeds) {
					// try to find largest seed size that still gets accepted
					size_t bestK = k;
					while (acceptSeed(sIdx, matchCount, extraOccs, k, nTax, SA, presenceChecker, taxonCoords, T, options, maxMatches)) {
						bestK = k;
						if (k == options.maxK || startPos + k + 1 >= T.size() || !presenceChecker.isFree(startPos + k)) { // no further extension of seed length, or newly added character would be already taken anyway
							break;
						}
						k++;
						matchCount = countMatches(sIdx, lcp, k);
					}
					k = bestK;
					matchCount = countMatches(sIdx, lcp, k);
				}

				if (options.maxMismatches > 0) {
					std::string pattern = T.substr(startPos, k);
					extraOccs = approxMatcher.findOccurrences(T, SA, presenceChecker, pattern, options.maxMismatches, 1, false);

					// TODO: Maybe only add those approximate matches that don't collide with the exact matches we already have?

					if (!acceptSeed(sIdx, matchCount, extraOccs, k, nTax, SA, presenceChecker, taxonCoords, T, options, maxMatches)) {
						if (k == options.maxK || startPos + k + 1 >= T.size() || !presenceChecker.isFree(startPos + k)) { // no further extension of seed length, or newly added character would be already taken anyway
							break;
						}
						k++;
						matchCount = countMatches(sIdx, lcp, k);
						continue;
					}
				}

				Seed block(nTax);
				for (size_t i = sIdx; i < sIdx + matchCount; ++i) {
					block.addTaxon(posToTaxon(SA[i], taxonCoords, T.size(), options.reverseComplement), SA[i], SA[i] + k - 1);
				}
				for (size_t i = 0; i < extraOccs.size(); ++i) {
					block.addTaxon(posToTaxon(extraOccs[i].first, taxonCoords, T.size(), options.reverseComplement), extraOccs[i].first,
							extraOccs[i].second);
				}
				computeBestCaseMaxSizes(block, T, presenceChecker, nTax);
#pragma omp critical
				{
					res.push_back(block);
					if (options.verboseDebug) {
						std::cout << "Pushing back a seeded block with " << block.getNTaxInBlock() << " taxa and seed size "
								<< block.getSeedSize() << "\n";
					}
				}
				if (!options.quartetFlavor) {
					double progress = (double) 100 * sIdx / SA.size(); // TODO: Fix this, this looks kinda wrong in parallel mode
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
				break;
			} else {
				if (k == options.maxK || startPos + k + 1 >= T.size() || !presenceChecker.isFree(startPos + k)) { // no further extension of seed length, or newly added character would be already taken anyway
					break;
				}
				k++;
				matchCount = countMatches(sIdx, lcp, k);
			}
		}
	}
	return res;
}

std::vector<Seed> filterSeededBlocks(std::vector<Seed>& seededBlocks, const std::string& T, size_t nTax, PresenceChecker& presenceChecker,
		const Options& options) {
	std::vector<Seed> res;
	for (size_t i = 0; i < seededBlocks.size(); ++i) {
		if (presenceChecker.isFine(seededBlocks[i])) {
			res.push_back(seededBlocks[i]);
			trivialExtension(seededBlocks[i], T, presenceChecker, nTax);
			presenceChecker.reserveSeededBlock(seededBlocks[i]);
		}
	}
	return res;
}

void extractExtendedBlocksPreselectSeeds(const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker, BlockWriter& writer,
		SummaryStatistics& stats, const Options& options, size_t minK, size_t nMin, size_t nMax) {
// this version assumes iterativeSeeding and preselectSeeds. It extends the blocks only after the seeds have been chosen.
	if (!options.quartetFlavor)
		std::cout << "Extracting seeded blocks MIAU...\n";
	std::vector<Seed> seededBlocks;
	for (size_t i = nMax; i >= nMin; i--) {
		std::vector<Seed> actSeededBlocks = extractSeededBlocks(concat.getConcatenatedSeq(), concat.nTax(), concat.getSuffixArray(),
				concat.getLcpArray(), presenceChecker, concat.getTaxonCoords(), options, minK, i, i);
		actSeededBlocks = filterSeededBlocks(actSeededBlocks, concat.getConcatenatedSeq(), concat.nTax(), presenceChecker, options);
		std::sort(actSeededBlocks.begin(), actSeededBlocks.end(), std::greater<Seed>());
		std::cout << "Found " << actSeededBlocks.size() << " new seeded blocks with " << i << " matches.\n";
		for (size_t j = 0; j < actSeededBlocks.size(); ++j) {
			seededBlocks.push_back(actSeededBlocks[j]);
		}
	}

// TODO: Remove me again, this is just out of curiosity
//std::vector<Superseed> superseeds = buildSuperseeds(seededBlocks, T, presenceChecker, nTax, options);

	std::cout << "seededBlocks.size(): " << seededBlocks.size() << "\n";
	if (!options.quartetFlavor)
		std::cout << "Assembling extended blocks...\n";
	double lastP = 0;
	std::sort(seededBlocks.begin(), seededBlocks.end(), std::greater<Seed>());
	for (size_t i = 0; i < seededBlocks.size(); ++i) {
		Seed seededBlock = seededBlocks[i];
//trivialExtension(seededBlock, concat.getConcatenatedSeq(), presenceChecker, concat.nTax());
		ExtendedBlock extendedBlock = extendBlock(seededBlock, concat.getConcatenatedSeq(), concat.nTax(), presenceChecker, options);
// check if the extended block can still be accepted.
		if (presenceChecker.isFineWithoutSeed(extendedBlock)) {
			presenceChecker.reserveExtendedBlock(extendedBlock);
			std::vector<std::string> msa = extendedBlock.msaWrapper.assembleMSA();
			if (options.verboseDebug) {
				std::cout << "Pushing back a block with alignment: \n";
				for (size_t i = 0; i < msa.size(); ++i) {
					std::cout << msa[i] << "\n";
				}
			}
			stats.updateSummaryStatistics(extendedBlock, concat.nTax());
			writer.writeTemporaryBlockMSA(extendedBlock, concat.nTax());
		}
		if (!options.quartetFlavor) {
			double progress = (double) 100 * i / seededBlocks.size();
			if (progress > lastP + 1) {
				std::cout << progress << " %\n";
				lastP = progress;
			}
		}
	}
}

void extractExtendedBlocks(const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker, BlockWriter& writer,
		SummaryStatistics& stats, const Options& options, size_t minK, size_t nMin, size_t nMax) {
	if (options.preselectSeeds) {
		return extractExtendedBlocksPreselectSeeds(concat, presenceChecker, writer, stats, options, minK, nMin, nMax);
	}

	if (!options.quartetFlavor)
		std::cout << "Extracting seeded blocks...\n";
	std::vector<Seed> seededBlocks;
	if (options.iterativeSeeding) {
		PresenceChecker seedingPresenceChecker(presenceChecker);
		for (size_t i = nMax; i >= nMin; i--) {
			std::vector<Seed> actSeededBlocks = extractSeededBlocks(concat.getConcatenatedSeq(), concat.nTax(), concat.getSuffixArray(),
					concat.getLcpArray(), seedingPresenceChecker, concat.getTaxonCoords(), options, minK, i, i);
			actSeededBlocks = filterSeededBlocks(actSeededBlocks, concat.getConcatenatedSeq(), concat.nTax(), seedingPresenceChecker,
					options);
			std::sort(actSeededBlocks.begin(), actSeededBlocks.end(), std::greater<Seed>());
			std::cout << "Found " << actSeededBlocks.size() << " new seeded blocks with " << i << " matches.\n";

			if (options.iterativeExtension) {
				double lastP = 0;
				for (size_t i = 0; i < actSeededBlocks.size(); ++i) {
					Seed seededBlock = actSeededBlocks[i];
					if (!presenceChecker.isFine(seededBlock))
						continue;
					//trivialExtension(seededBlock, T, presenceChecker, nTax);
					ExtendedBlock extendedBlock = extendBlock(seededBlock, concat.getConcatenatedSeq(), concat.nTax(), presenceChecker,
							options);
					// check if the extended block can still be accepted.
					if (presenceChecker.isFine(extendedBlock)) {
						presenceChecker.reserveExtendedBlock(extendedBlock);
						seedingPresenceChecker.reserveExtendedBlock(extendedBlock);
						std::vector<std::string> msa = extendedBlock.msaWrapper.assembleMSA();
						if (options.verboseDebug) {
							std::cout << "Pushing back a block with alignment: \n";
							for (size_t i = 0; i < msa.size(); ++i) {
								std::cout << msa[i] << "\n";
							}
						}
						stats.updateSummaryStatistics(extendedBlock, concat.nTax());
						writer.writeTemporaryBlockMSA(extendedBlock, concat.nTax());
					}
					if (!options.quartetFlavor) {
						double progress = (double) 100 * i / actSeededBlocks.size();
						if (progress > lastP + 1) {
							std::cout << progress << " %\n";
							lastP = progress;
						}
					}
				}
				std::cout << "Finished extension of the newly discovered seeded blocks with at least " << i << " matches.\n";
				std::cout << "Current extendedBlocks.size(): " << stats.getCurrentNBlocks() << "\n";
			} else {
				for (size_t j = 0; j < actSeededBlocks.size(); ++j) {
					seededBlocks.push_back(actSeededBlocks[j]);
				}
			}
		}
	} else {
		seededBlocks = extractSeededBlocks(concat.getConcatenatedSeq(), concat.nTax(), concat.getSuffixArray(), concat.getLcpArray(),
				presenceChecker, concat.getTaxonCoords(), options, minK, nMin, nMax);
	}

// TODO: Remove me again, this is just out of curiosity
//std::vector<Superseed> superseeds = buildSuperseeds(seededBlocks, T, presenceChecker, nTax, options);

	if (!options.iterativeExtension) {
		std::cout << "seededBlocks.size(): " << seededBlocks.size() << "\n";
		if (!options.quartetFlavor)
			std::cout << "Assembling extended blocks...\n";
		double lastP = 0;
		std::sort(seededBlocks.begin(), seededBlocks.end(), std::greater<Seed>());
		for (size_t i = 0; i < seededBlocks.size(); ++i) {
			Seed seededBlock = seededBlocks[i];
			if (!presenceChecker.isFine(seededBlock))
				continue;
			trivialExtension(seededBlock, concat.getConcatenatedSeq(), presenceChecker, concat.nTax());
			ExtendedBlock extendedBlock = extendBlock(seededBlock, concat.getConcatenatedSeq(), concat.nTax(), presenceChecker, options);
			// check if the extended block can still be accepted.
			if (presenceChecker.isFine(extendedBlock)) {
				presenceChecker.reserveExtendedBlock(extendedBlock);
				std::vector<std::string> msa = extendedBlock.msaWrapper.assembleMSA();
				if (options.verboseDebug) {
					std::cout << "Pushing back a block with alignment: \n";
					for (size_t i = 0; i < msa.size(); ++i) {
						std::cout << msa[i] << "\n";
					}
				}
				stats.updateSummaryStatistics(extendedBlock, concat.nTax());
				writer.writeTemporaryBlockMSA(extendedBlock, concat.nTax());
			}
			if (!options.quartetFlavor) {
				double progress = (double) 100 * i / seededBlocks.size();
				if (progress > lastP + 1) {
					std::cout << progress << " %\n";
					lastP = progress;
				}
			}
		}
	}

}

std::vector<ExtendedBlock> extractExtendedBlocks(const std::string& T, size_t nTax, const std::vector<size_t>& SA,
		const std::vector<size_t>& lcp, PresenceChecker& presenceChecker, const std::vector<IndexedTaxonCoords>& taxonCoords,
		const Options& options, size_t minK, size_t nMin, size_t nMax) {
	std::vector<ExtendedBlock> res;
	if (!options.quartetFlavor)
		std::cout << "Extracting seeded blocks...\n";

	std::vector<Seed> seededBlocks;
	if (options.iterativeSeeding) {
		PresenceChecker seedingPresenceChecker(presenceChecker);
		for (size_t i = nMax; i >= nMin; i--) {
			std::vector<Seed> actSeededBlocks = extractSeededBlocks(T, nTax, SA, lcp, seedingPresenceChecker, taxonCoords, options, minK, i, i);
			actSeededBlocks = filterSeededBlocks(actSeededBlocks, T, nTax, seedingPresenceChecker, options);
			std::sort(actSeededBlocks.begin(), actSeededBlocks.end(), std::greater<Seed>());
			std::cout << "Found " << actSeededBlocks.size() << " new seeded blocks with at least " << i << " matches.\n";

			if (options.iterativeExtension) {
				double lastP = 0;
				for (size_t i = 0; i < actSeededBlocks.size(); ++i) {
					Seed seededBlock = actSeededBlocks[i];
					if (!presenceChecker.isFine(seededBlock))
						continue;
					//trivialExtension(seededBlock, T, presenceChecker, nTax);
					ExtendedBlock extendedBlock = extendBlock(seededBlock, T, nTax, presenceChecker, options);
					// check if the extended block can still be accepted.
					if (presenceChecker.isFine(extendedBlock)) {
						presenceChecker.reserveExtendedBlock(extendedBlock);
						seedingPresenceChecker.reserveExtendedBlock(extendedBlock);
						std::vector<std::string> msa = extendedBlock.msaWrapper.assembleMSA();
						if (options.verboseDebug) {
							std::cout << "Pushing back a block with alignment: \n";
							for (size_t i = 0; i < msa.size(); ++i) {
								std::cout << msa[i] << "\n";
							}
						}
						res.push_back(extendedBlock);
					}
					if (!options.quartetFlavor) {
						double progress = (double) 100 * i / actSeededBlocks.size();
						if (progress > lastP + 1) {
							std::cout << progress << " %\n";
							lastP = progress;
						}
					}
				}
				std::cout << "Finished extension of the newly discovered seeded blocks with at least " << i << " matches.\n";
				std::cout << "Current extendedBlocks.size(): " << res.size() << "\n";
			} else {
				for (size_t j = 0; j < actSeededBlocks.size(); ++j) {
					seededBlocks.push_back(actSeededBlocks[j]);
				}
			}
		}
	} else {
		seededBlocks = extractSeededBlocks(T, nTax, SA, lcp, presenceChecker, taxonCoords, options, nMin, nTax,
				minK);
	}

// TODO: Remove me again, this is just out of curiosity
//std::vector<Superseed> superseeds = buildSuperseeds(seededBlocks, T, presenceChecker, nTax, options);

	if (!options.iterativeExtension) {
		std::cout << "seededBlocks.size(): " << seededBlocks.size() << "\n";
		if (!options.quartetFlavor)
			std::cout << "Assembling extended blocks...\n";
		double lastP = 0;
		std::sort(seededBlocks.begin(), seededBlocks.end(), std::greater<Seed>());
		for (size_t i = 0; i < seededBlocks.size(); ++i) {
			Seed seededBlock = seededBlocks[i];
			if (!presenceChecker.isFine(seededBlock))
				continue;
			trivialExtension(seededBlock, T, presenceChecker, nTax);
			ExtendedBlock extendedBlock = extendBlock(seededBlock, T, nTax, presenceChecker, options);
			// check if the extended block can still be accepted.
			if (presenceChecker.isFine(extendedBlock)) {
				presenceChecker.reserveExtendedBlock(extendedBlock);
				std::vector<std::string> msa = extendedBlock.msaWrapper.assembleMSA();
				if (options.verboseDebug) {
					std::cout << "Pushing back a block with alignment: \n";
					for (size_t i = 0; i < msa.size(); ++i) {
						std::cout << msa[i] << "\n";
					}
				}
				res.push_back(extendedBlock);
			}
			if (!options.quartetFlavor) {
				double progress = (double) 100 * i / seededBlocks.size();
				if (progress > lastP + 1) {
					std::cout << progress << " %\n";
					lastP = progress;
				}
			}
		}
	}
	return res;
}

void extractExtendedBlocksDecreasingNminMinK(const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker,
		BlockWriter& writer, SummaryStatistics& stats, const Options& options) {
	size_t nMin = concat.nTax();
	while (nMin >= options.minTaxaPerBlock) {
		size_t minK = findLargestK(concat.getLcpArray(), concat.getSuffixArray(), concat.getTaxonCoords(), concat.getConcatSize(),
				options.reverseComplement, presenceChecker);
		minK = std::min(minK, options.maxK);
		while (minK >= options.minK) {
			std::cout << "current minK: " << minK << "\n";
			extractExtendedBlocks(concat, presenceChecker, writer, stats, options, minK, nMin, nMin);
			minK--;
		}
		std::cout << "current number of extended blocks with nMin >= " << nMin << ": " << stats.getCurrentNBlocks() << "\n";
		nMin--;
	}
}
