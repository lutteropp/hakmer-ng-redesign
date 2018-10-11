/*
 * msa_wrapper.hpp
 *
 *  Created on: Oct 11, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <vector>
#include <string>

#include "nogaps_msa.hpp"
#include "star_msa.hpp"

class MSAWrapper {
public:
	MSAWrapper(bool noGaps) :
			noGaps(noGaps) {
	}
	void init(size_t nTax) {
		if (noGaps) {
			noGapsMSA.init(nTax);
		} else {
			starMSA.init(nTax);
		}
	}
	std::vector<std::string> assembleMSA() {
		if (noGaps) {
			return noGapsMSA.assembleMSA();
		} else {
			return starMSA.assembleMSA();
		}
	}
	double pairwiseDistance(size_t idx1, size_t idx2) {
		if (noGaps) {
			return noGapsMSA.pairwiseDistance(idx1, idx2);
		} else {
			return starMSA.pairwiseDistance(idx1, idx2);
		}
	}
	double normalizedPairwiseDistance(size_t idx1, size_t idx2) {
		if (noGaps) {
			return noGapsMSA.normalizedPairwiseDistance(idx1, idx2);
		} else {
			return starMSA.normalizedPairwiseDistance(idx1, idx2);
		}
	}
	void addCharsLeft(const std::vector<char>& chars) {
		if (noGaps) {
			noGapsMSA.addCharsLeft(chars);
		} else {
			starMSA.addCharsLeft(chars);
		}
	}
	void addCharsRight(const std::vector<char>& chars) {
		if (noGaps) {
			noGapsMSA.addCharsRight(chars);
		} else {
			starMSA.addCharsRight(chars);
		}
	}
	void setSeeds(const std::vector<std::string>& seeds) {
		if (noGaps) {
			noGapsMSA.setSeeds(seeds);
		} else {
			starMSA.setSeeds(seeds);
		}
	}
	void setSeeds(const std::string& seed) {
		if (noGaps) {
			noGapsMSA.setSeeds(seed);
		} else {
			starMSA.setSeeds(seed);
		}
	}
	void shrinkDownToLeftFlank(size_t newLeftFlankSize) {
		if (noGaps) {
			noGapsMSA.shrinkDownToLeftFlank(newLeftFlankSize);
		} else {
			starMSA.shrinkDownToLeftFlank(newLeftFlankSize);
		}
	}
	void shrinkDownToRightFlank(size_t newRightFlankSize) {
		if (noGaps) {
			noGapsMSA.shrinkDownToRightFlank(newRightFlankSize);
		} else {
			starMSA.shrinkDownToRightFlank(newRightFlankSize);
		}
	}
	size_t getAlignmentWidth() {
		if (noGaps) {
			return noGapsMSA.getAlignmentWidth();
		} else {
			return starMSA.getAlignmentWidth();
		}
	}
private:
	bool noGaps;
	NoGapsMSA noGapsMSA;
	StarMSA starMSA;
};
