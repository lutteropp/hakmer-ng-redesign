/*
 * block_extension.hpp
 *
 *  Created on: Nov 7, 2018
 *      Author: sarah
 */

#pragma once

#include <vector>
#include <string>
#include <utility>

#include "options.hpp"
#include "seed.hpp"
#include "extended_block.hpp"
#include "presence_checker.hpp"

std::pair<size_t, size_t> computeBestCaseMaxSizes(Seed& seededBlock, const std::string& T, PresenceChecker& presenceChecker, size_t nTax); // TODO: Do we still need this?

void trivialExtension(Seed& seededBlock, const std::string& T, PresenceChecker& presenceChecker, size_t nTax);

ExtendedBlock extendBlock(const Seed& seededBlock, const std::string& T, size_t nTax, PresenceChecker& presenceChecker,
		const Options& options);