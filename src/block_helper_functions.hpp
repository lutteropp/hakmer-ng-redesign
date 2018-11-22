/*
 * block_helper.hpp
 *
 *  Created on: Sep 4, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <stddef.h>
#include <string>
#include "extended_block.hpp"
#include "seed.hpp"

std::string extractTaxonSequence(ExtendedBlock& block, size_t taxID);
std::string extractTaxonSequence(const Seed& block, size_t taxID, const std::string& T);
