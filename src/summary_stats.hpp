/*
 * summary_stats.hpp
 *
 *  Created on: Nov 12, 2018
 *      Author: sarah
 */

#pragma once

#include <stddef.h>
#include <fstream>
#include <iostream>
#include <string>

#include "extended_block.hpp"
#include "options.hpp"

class SummaryStatistics {
public:
	SummaryStatistics(size_t nTax) {
		taxonUsage.resize(nTax);
		for (size_t i = 0; i < nTax; ++i) {
			taxonUsage[i] = 0;
		}
	}

	void updateSummaryStatistics(const ExtendedBlock& block, size_t nTax) {
		size_t rowSize = block.getAverageSeedSize() + block.getAverageLeftFlankSize() + block.getAverageRightFlankSize();
		//seqDataUsed += rowSize * block.getNTaxInBlock();

		seqDataUsed += block.getTotalBasesUsed();

		nMissingData += rowSize * (nTax - block.getNTaxInBlock());
		leftFlankSum += block.getAverageLeftFlankSize();
		rightFlankSum += block.getAverageRightFlankSize();
		seedSizeSum += block.getAverageSeedSize();
		nTaxSum += block.getNTaxInBlock();

		for (size_t tID : block.getTaxonIDsInBlock()) {
			taxonUsage[tID] += block.getSizeWithFlanks(tID);
		}

		nBlocks++;
	}
	void printSummaryStatistics(const IndexedConcatenatedSequence& concat, const Options& options) {
		size_t nTax = concat.nTax();
		size_t totalSeqData = concat.getConcatSize();

		std::cout << "Number of extracted blocks: " << nBlocks << "\n";
		std::cout << "Percentage of total reconstructed sequence data: " << ((double) seqDataUsed * 100) / totalSeqData << " %\n";
		std::cout << "Percentage of missing data: " << ((double) nMissingData * 100) / (nMissingData + seqDataUsed) << " %\n";
		std::cout << "Average left flank size: " << (double) leftFlankSum / nBlocks << "\n";
		std::cout << "Average right flank size: " << (double) rightFlankSum / nBlocks << "\n";
		std::cout << "Average seed size: " << (double) seedSizeSum / nBlocks << "\n";
		std::cout << "Average nTax in block: " << (double) nTaxSum / nBlocks << "\n";

		// print percentage of reconstructed sequence data for each taxon
		std::cout << "\nPercentage of reconstructed sequence data per taxon:\n";
		for (size_t i = 0; i < nTax; ++i) {
			double p = ((double) taxonUsage[i] * 100) / concat.getTaxonCoords(i).getTotalLength();
			std::cout << concat.getTaxonLabels()[i] << ": " << p << " %\n";
		}
		std::cout << "\n";

		if (!options.infopath.empty()) {
			std::ofstream info(options.infopath);
			info << "Number of extracted blocks: " << nBlocks << "\n";
			info << "Percentage of total reconstructed sequence data: " << ((double) seqDataUsed * 100) / totalSeqData << " %\n";
			info << "Percentage of missing data: " << ((double) nMissingData * 100) / (nMissingData + seqDataUsed) << " %\n";
			info << "Average left flank size: " << (double) leftFlankSum / nBlocks << "\n";
			info << "Average right flank size: " << (double) rightFlankSum / nBlocks << "\n";
			info << "Average seed size: " << (double) seedSizeSum / nBlocks << "\n";
			info << "Average nTax in block: " << (double) nTaxSum / nBlocks << "\n";

			// print percentage of reconstructed sequence data for each taxon
			info << "\nPercentage of reconstructed sequence data per taxon:\n";
			for (size_t i = 0; i < nTax; ++i) {
				double p = ((double) taxonUsage[i] * 100) / concat.getTaxonCoords(i).getTotalLength();
				info << concat.getTaxonLabels()[i] << ": " << p << " %\n";
			}

			info.close();
		}
	}
	double getAmountSeqDataUsed(size_t totalSeqData) const {
		return ((double) seqDataUsed) / totalSeqData;
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
	std::vector<size_t> taxonUsage;
};
