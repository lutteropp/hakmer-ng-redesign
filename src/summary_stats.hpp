/*
 * summary_stats.hpp
 *
 *  Created on: Nov 12, 2018
 *      Author: sarah
 */

#pragma once

#include <cstdlib>
#include <iostream>
#include "extended_block.hpp"
#include "options.hpp"

class SummaryStatistics {
public:
	void updateSummaryStatistics(const ExtendedBlock& block, size_t nTax) {
		size_t rowSize = block.getSeedSize() + block.getLeftFlankSize() + block.getRightFlankSize();
		seqDataUsed += rowSize * block.getNTaxInBlock();
		nMissingData += rowSize * (nTax - block.getNTaxInBlock());
		leftFlankSum += block.getLeftFlankSize();
		rightFlankSum += block.getRightFlankSize();
		seedSizeSum += block.getSeedSize();
		nTaxSum += block.getNTaxInBlock();
		nBlocks++;
	}
	void printSummaryStatistics(size_t nTax, size_t totalSeqData, const Options& options) {
		std::cout << "Number of extracted blocks: " << nBlocks << "\n";
		std::cout << "Percentage of reconstructed sequence data: " << ((double) seqDataUsed * 100) / totalSeqData << " %\n";
		std::cout << "Percentage of missing data: " << ((double) nMissingData * 100) / (nMissingData + seqDataUsed) << " %\n";
		std::cout << "Average left flank size: " << (double) leftFlankSum / nBlocks << "\n";
		std::cout << "Average right flank size: " << (double) rightFlankSum / nBlocks << "\n";
		std::cout << "Average seed size: " << (double) seedSizeSum / nBlocks << "\n";
		std::cout << "Average nTax in block: " << (double) nTaxSum / nBlocks << "\n";

		if (!options.infopath.empty()) {
			std::ofstream info(options.infopath);
			info << "Number of extracted blocks: " << nBlocks << "\n";
			info << "Percentage of reconstructed sequence data: " << ((double) seqDataUsed * 100) / totalSeqData << " %\n";
			info << "Percentage of missing data: " << ((double) nMissingData * 100) / (nMissingData + seqDataUsed) << " %\n";
			info << "Average left flank size: " << (double) leftFlankSum / nBlocks << "\n";
			info << "Average right flank size: " << (double) rightFlankSum / nBlocks << "\n";
			info << "Average seed size: " << (double) seedSizeSum / nBlocks << "\n";
			info << "Average nTax in block: " << (double) nTaxSum / nBlocks << "\n";
			info.close();
		}
	}
	size_t getCurrentNBlocks() const {
		return nBlocks;
	}
private:
	size_t seqDataUsed = 0;
	size_t nMissingData = 0;
	size_t rightFlankSum = 0;
	size_t leftFlankSum = 0;
	size_t seedSizeSum = 0;
	size_t nTaxSum = 0;
	size_t nBlocks = 0;
};
