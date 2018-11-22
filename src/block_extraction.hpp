/*
 * block_extraction.hpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <vector>
#include <cstdlib>
#include <string>

#include "presence_checker.hpp"
#include "options.hpp"
#include "extended_block.hpp"
#include "indexed_concat.hpp"
#include "seed.hpp"
#include "block_writer.hpp"
#include "summary_stats.hpp"

std::vector<Seed> extractSeededBlocks(const std::string& T, size_t nTax, const std::vector<size_t>& SA, const std::vector<size_t>& lcp,
		PresenceChecker& presenceChecker, const std::vector<IndexedTaxonCoords>& taxonCoords, const Options& options);
Seed nextSeededBlock(size_t& actSAPos, const std::string& T, size_t nTax, const std::vector<size_t>& SA, const std::vector<size_t>& lcp,
		PresenceChecker& presenceChecker, const std::vector<IndexedTaxonCoords>& taxonCoords, const Options& options);

ExtendedBlock extendBlock(const Seed& block, const std::string& T, size_t nTax, PresenceChecker& presenceChecker, const Options& options);
std::vector<ExtendedBlock> extractExtendedBlocks(const std::string& T, size_t nTax, const std::vector<size_t>& SA,
		const std::vector<size_t>& lcp, PresenceChecker& presenceChecker, const std::vector<IndexedTaxonCoords>& taxonCoords,
		const Options& options, size_t minK, size_t nMin, size_t nMax);

void extractExtendedBlocks(const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker, BlockWriter& writer,
		SummaryStatistics& stats, const Options& options, size_t minK, size_t nMin, size_t nMax);

void extractExtendedBlocksDecreasingNminMinK(const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker, BlockWriter& writer,
		SummaryStatistics& stats, const Options& options);
