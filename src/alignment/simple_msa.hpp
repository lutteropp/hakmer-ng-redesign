/*
 * simple_msa.hpp
 *
 *  Created on: Nov 14, 2018
 *      Author: sarah
 */

#pragma once

#include <stddef.h>
#include <string>
#include <vector>

#include "simple_coords.hpp"
#include "../extended_block.hpp"
#include "../options.hpp"

std::vector<std::string> computeMSA(const std::vector<std::string>& seqs);
std::vector<std::string> computeMSA(const std::vector<SimpleCoords>& seqCoords, const std::string& T, size_t nTax);
std::vector<std::string> computeMSA(const ExtendedBlock& block, const std::string& T, size_t nTax, const Options& options);
