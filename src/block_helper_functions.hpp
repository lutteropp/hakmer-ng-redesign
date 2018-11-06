/*
 * block_helper.hpp
 *
 *  Created on: Sep 4, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include "extended_block.hpp"
#include "seed.hpp"

std::string createMissingString(size_t len);
std::string extractTaxonSequence(ExtendedBlock& block, size_t taxID);
std::string extractTaxonSequence(const Seed& block, size_t taxID, const std::string& T);
