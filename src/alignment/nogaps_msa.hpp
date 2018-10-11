/*
 * nogaps_msa.hpp
 *
 *  Created on: Oct 9, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <string>
#include <vector>
#include "twodim_matrix.hpp"

class NoGapsMSA {
public:
	NoGapsMSA();
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
	size_t getAlignmentWidth() const;
private:
	bool isGapCharacter(char c);
	size_t nTax;
	TwoDimMatrix<size_t> pairwiseHammingDistances;
	std::vector<std::string> sequencesLeft, sequencesMiddle, sequencesRight;
	size_t width;
	std::vector<std::string> msa;
	bool msaValid;
};
