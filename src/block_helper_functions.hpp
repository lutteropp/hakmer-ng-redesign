/*
 * block_helper.hpp
 *
 *  Created on: Sep 4, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include "aligned_block.hpp"
#include "extended_block.hpp"
#include "seeded_block.hpp"

std::string extractTaxonSequence(const AlignedBlock& block, size_t taxID);
std::string extractTaxonSequence(const ExtendedBlock& block, size_t taxID, const std::string& T);
std::string extractTaxonSequence(const SeededBlock& block, size_t taxID, const std::string& T);
std::string createMissingString(size_t len);
