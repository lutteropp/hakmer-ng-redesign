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

inline std::vector<std::string> computeMSA(const std::vector<std::string>& seqs);
inline std::vector<std::string> computeMSA(const std::vector<SimpleCoords>& seqCoords, const std::string& T, size_t nTax);
inline std::vector<std::string> computeMSA(const ExtendedBlock& block, const std::string& T, size_t nTax);
