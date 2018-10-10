/*
 * star_msa.hpp
 *
 *  Created on: Oct 9, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <vector>
#include <string>
#include "split_pairwise_alignment.hpp"

// 1.) Build all O(n^2) pairwise alignments.
// 2.) Choose center sequence S that is closest to all other sequences (this is, sum of edit distances to the other sequences is minimal).
// 3.) Incrementally merge the pairwise alignments into one MSA (Progressive alignment). Gaps in S from the pairwise alignment stay, insert gaps to the other sequences as needed.

class StarMSA {
public:
	StarMSA();
	void init(size_t nTax);
	std::vector<std::string> assembleMSA();
	double pairwiseDistance(size_t idx1, size_t idx2);
	double normalizedPairwiseDistance(size_t idx1, size_t idx2);
	void addCharsLeft(const std::vector<char>& chars);
	void addCharsRight(const std::vector<char>& chars);
	void setSeeds(const std::vector<std::string>& seeds);
	void setSeeds(const std::string& seed);
	void shrinkDownToLeftFlank(size_t newLeftFlankSize);
	void shrinkDownToRightFlank(size_t newRightFlankSize);
private:
	void addToMSA(size_t taxonToAdd, std::vector<std::string>& msa, size_t centerSequenceIdx);
	size_t nTax;
	TwoDimMatrix<SplitPairwiseAlignment> pairwiseAlignments;
	std::vector<std::string> msa;
	bool msaValid;
};
