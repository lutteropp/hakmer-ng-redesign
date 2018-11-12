/*
 * block_writer.hpp
 *
 *  Created on: Nov 12, 2018
 *      Author: sarah
 */

#pragma once

#include <fstream>
#include <vector>
#include <string>

#include "extended_block.hpp"

class BlockWriter {
public:
	BlockWriter(size_t nTax);
	void writeTemporaryBlockMSA(ExtendedBlock& block, size_t nTax);
	void assembleFinalSupermatrix(const std::vector<std::string>& taxonLabels, const std::string& outpath);
private:
	void deleteTemporaryFiles();
	std::vector<std::ofstream> tempMSAFiles;
};
