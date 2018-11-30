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

void extractExtendedBlocks(const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker, BlockWriter& writer,
		SummaryStatistics& stats, const Options& options, size_t minK, size_t maxK, size_t flankWidth);
