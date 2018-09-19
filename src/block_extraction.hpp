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
#include "seeded_block.hpp"
#include "extended_block.hpp"
#include "aligned_block.hpp"
#include "indexed_concat.hpp"

std::pair<std::vector<size_t>, std::vector<size_t> > shrinkArrays(const IndexedConcatenatedSequence& concat, const std::vector<std::pair<size_t, size_t> >& taxonCoords,
		const std::vector<size_t>& wantedTaxa, const Options& options);

std::vector<SeededBlock> extractSeededBlocks(const std::string& T, size_t nTax, const std::vector<size_t>& SA,
		const std::vector<size_t>& lcp, PresenceChecker& presenceChecker, const std::vector<std::pair<size_t, size_t> >& taxonCoords,
		const Options& options);
SeededBlock nextSeededBlock(size_t& actSAPos, const std::string& T, size_t nTax, const std::vector<size_t>& SA,
		const std::vector<size_t>& lcp, PresenceChecker& presenceChecker, const std::vector<std::pair<size_t, size_t> >& taxonCoords,
		const Options& options);

ExtendedBlock extendBlock(const SeededBlock& block, const std::string& T, size_t nTax, PresenceChecker& presenceChecker,
		const Options& options);
std::vector<ExtendedBlock> extractExtendedBlocks(const std::string& T, size_t nTax, const std::vector<size_t>& SA,
		const std::vector<size_t>& lcp, PresenceChecker& presenceChecker, const std::vector<std::pair<size_t, size_t> >& taxonCoords,
		const Options& options);

std::vector<AlignedBlock> extractAlignedBlocks(const std::string& T, size_t nTax, const std::vector<size_t>& SA,
		const std::vector<size_t>& lcp, PresenceChecker& presenceChecker, const std::vector<std::pair<size_t, size_t> >& taxonCoords,
		const Options& options);
