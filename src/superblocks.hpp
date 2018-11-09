/*
 * superblocks.hpp
 *
 *  Created on: Nov 7, 2018
 *      Author: sarah
 */

#pragma once

#include "options.hpp"
#include "superseed.hpp"
#include "presence_checker.hpp"

std::vector<Superseed> buildSuperseeds(const std::vector<Seed>& seeds, const std::string& T, const PresenceChecker& presenceChecker, size_t nTax, const Options& options);
