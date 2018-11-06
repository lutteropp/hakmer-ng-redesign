/*
 * split_pairwise_alignment.hpp
 *
 *  Created on: Oct 9, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <string>

#include "pairwise_alignment.hpp"

class SplitPairwiseAlignment {
public:
	SplitPairwiseAlignment();
	void addCharsLeft(char a, char b);
	void addCharsRight(char a, char b);
	void setSeed(const std::string& seed1, const std::string& seed2);
	void setSeed(const std::string& seed);
	std::pair<std::string, std::string> extractAlignment();
	void printAlignment();
	size_t getAlignmentWidth();
	std::string getS1() const;
	std::string getS2() const;
	double pairwiseDistance();
	void shrinkDownToLeftFlank(size_t newLeftFlankSize);
	void shrinkDownToRightFlank(size_t newRightFlankSize);
	std::pair<std::string, std::string> extractRightFlankAlignment();
	std::pair<std::string, std::string> extractReversedLeftFlankAlignment();
private:
	bool singleSeed;
	std::string singleSeedSequence;
	PairwiseAlignment leftFlank;
	PairwiseAlignment seed;
	PairwiseAlignment rightFlank;
	std::pair<std::string, std::string> ali;
	bool aliValid;
};
