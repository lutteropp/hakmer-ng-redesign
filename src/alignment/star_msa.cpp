/*
 * star_msa.cpp
 *
 *  Created on: Oct 9, 2018
 *      Author: Sarah Lutteropp
 */

#include "star_msa.hpp"

#include <stdexcept>
#include <iostream>

StarMSA::StarMSA() :
		nTax(0), msaValid(false) {
}
void StarMSA::init(size_t nTax) {
	this->nTax = nTax;
	pairwiseAlignments.init(nTax, nTax);
	normalizedPairwiseDistances.init(nTax, nTax);
}
std::vector<std::string> StarMSA::assembleMSA() {
	if (msaValid) {
		return msa;
	}
	msa.resize(nTax);
	// first, find the sequence with the smallest distance sum to the others
	std::vector<int> distanceSums(nTax, 0);
	for (size_t i = 0; i < nTax; ++i) {
		for (size_t j = i + 1; j < nTax; ++j) {
			distanceSums[i] += pairwiseAlignments.entryAt(i, j).pairwiseDistance();
			distanceSums[j] += pairwiseAlignments.entryAt(i, j).pairwiseDistance();
		}
	}
	size_t smallestIdx = 0;
	int smallestDist = distanceSums[0];
	for (size_t i = 1; i < nTax; ++i) {
		if (distanceSums[i] < smallestDist) {
			smallestIdx = i;
			smallestDist = distanceSums[i];
		}
	}

	for (size_t i = 0; i < nTax; ++i) {
		if (i == smallestIdx)
			continue;
		addToMSA(i, msa, smallestIdx);
	}

	// perform a sanity check for the MSA
	size_t w = 0;
	for (size_t i = 0; i < msa.size(); ++i) {
		if (msa[i].empty())
			continue;
		if (w == 0) {
			w = msa[i].size();
		}
		if (msa[i].size() != w) {
			std::cout << "w: " << w << "\n";
			std::cout << "msa[i].size(): " << msa[i].size() << "\n";
			std::cout << "i: " << i << "\n";

			std::cout << "The entire MSA:\n";
			for (size_t j = 0; j < msa.size(); ++j) {
				std::cout << msa[j] << "\n";
			}
			std::cout << "\n";
			std::cout << "The pairwise MSAs:\n";
			for (size_t j = 0; j < pairwiseAlignments.size(); ++j) {
				for (size_t k = j + 1; k < pairwiseAlignments.size(); ++k) {
					pairwiseAlignments.entryAt(j, k).printAlignment();
				}
			}
			throw std::runtime_error("THIS SHOULD NOT HAPPEN! ALIGNMENT WIDTH IS NOT EQUAL");
		}
	}
	msaValid = true;

	// rescue the normalized pairwise distances
	for (size_t i = 0; i < nTax; ++i) {
		for (size_t j = 0; j < nTax; ++j) {
			normalizedPairwiseDistances.entryAt(i, j) = normalizedPairwiseDistance(i, j);
		}
	}
	// now the pairwise alignments are not needed anymore, most likely not even the pairwise distances... free them.
	pairwiseAlignments.shrinkDownTo(0, 0);

	return msa;
}
double StarMSA::pairwiseDistance(size_t idx1, size_t idx2) {
	if (idx1 == idx2)
		return 0.0;
	size_t firstIdx = std::min(idx1, idx2);
	size_t secondIdx = std::max(idx1, idx2);
	return pairwiseAlignments.entryAt(firstIdx, secondIdx).pairwiseDistance();
}
double StarMSA::normalizedPairwiseDistance(size_t idx1, size_t idx2) {
	if (msaValid) {
		return normalizedPairwiseDistances.entryAt(idx1, idx2);
	}

	if (idx1 == idx2)
		return 0.0;
	size_t firstIdx = std::min(idx1, idx2);
	size_t secondIdx = std::max(idx1, idx2);
	double editDist = pairwiseAlignments.entryAt(firstIdx, secondIdx).pairwiseDistance();
	size_t s1Size = pairwiseAlignments.entryAt(firstIdx, secondIdx).getS1().size();
	size_t s2Size = pairwiseAlignments.entryAt(firstIdx, secondIdx).getS2().size();
	double normEditDist = (2.0 * editDist) / (GAP_PENALTY * (s1Size + s2Size) + editDist);
	return normEditDist;
}

void StarMSA::addCharsLeft(const std::vector<char>& chars) {
	for (size_t i = 0; i < nTax - 1; ++i) {
		for (size_t j = i + 1; j < nTax; ++j) {
			pairwiseAlignments.entryAt(i, j).addCharsLeft(chars[i], chars[j]);
		}
	}
	msaValid = false;
}
void StarMSA::addCharsRight(const std::vector<char>& chars) {
	for (size_t i = 0; i < nTax - 1; ++i) {
		for (size_t j = i + 1; j < nTax; ++j) {
			pairwiseAlignments.entryAt(i, j).addCharsRight(chars[i], chars[j]);
		}
	}
	msaValid = false;
}
void StarMSA::setSeeds(const std::vector<std::string>& seeds) {
	for (size_t i = 0; i < nTax - 1; ++i) {
		for (size_t j = i + 1; j < nTax; ++j) {
			pairwiseAlignments.entryAt(i, j).setSeed(seeds[i], seeds[j]);
		}
	}
	msaValid = false;
}
void StarMSA::setSeeds(const std::string& seed) {
	for (size_t i = 0; i < nTax - 1; ++i) {
		for (size_t j = i + 1; j < nTax; ++j) {
			pairwiseAlignments.entryAt(i, j).setSeed(seed);
		}
	}
	msaValid = false;
}
void StarMSA::shrinkDownToLeftFlank(size_t newLeftFlankSize) {
	for (size_t i = 0; i < nTax - 1; ++i) {
		for (size_t j = i + 1; j < nTax; ++j) {
			pairwiseAlignments.entryAt(i, j).shrinkDownToLeftFlank(newLeftFlankSize);
		}
	}
	msaValid = false;
}
void StarMSA::shrinkDownToRightFlank(size_t newRightFlankSize) {
	for (size_t i = 0; i < nTax - 1; ++i) {
		for (size_t j = i + 1; j < nTax; ++j) {
			pairwiseAlignments.entryAt(i, j).shrinkDownToRightFlank(newRightFlankSize);
		}
	}
	msaValid = false;
}

size_t StarMSA::getAlignmentWidth() {
	if (!msaValid) {
		assembleMSA();
	}
	for (size_t i = 0; i < msa.size(); ++i) {
		if (!msa[i].empty()) {
			return msa[i].size();
		}
	}
	return 0;
}

void StarMSA::addToMSA(size_t taxonToAdd, std::vector<std::string>& msa, size_t centerSequenceIdx) {
	if (taxonToAdd == centerSequenceIdx) {
		throw std::runtime_error("This should not happen: taxonToAdd == centerSequenceIndex here");
	}

	size_t firstIdx = std::min(taxonToAdd, centerSequenceIdx);
	size_t secondIdx = std::max(taxonToAdd, centerSequenceIdx);
	std::pair<std::string, std::string> alignmentToAdd = pairwiseAlignments.entryAt(firstIdx, secondIdx).extractAlignment();

	std::string aliSeqCenter;
	std::string aliSeqNewTaxon;
	if (centerSequenceIdx == firstIdx) {
		aliSeqCenter = alignmentToAdd.first;
		aliSeqNewTaxon = alignmentToAdd.second;
	} else {
		aliSeqCenter = alignmentToAdd.second;
		aliSeqNewTaxon = alignmentToAdd.first;
	}

	if (aliSeqCenter.size() != aliSeqNewTaxon.size()) {
		throw std::runtime_error("This should not happen: The sequences in the pairwise alignment have different sizes");
	}

	if (msa[centerSequenceIdx].empty()) {
		msa[centerSequenceIdx] += aliSeqCenter;
		msa[taxonToAdd] += aliSeqNewTaxon;
		return;
	}

	size_t i_ali = 0;
	size_t i_msa = 0;

	while (i_ali < aliSeqCenter.size() && i_msa < msa[centerSequenceIdx].size()) {
		if (msa[centerSequenceIdx][i_msa] == aliSeqCenter[i_ali]) {
			msa[taxonToAdd] += aliSeqNewTaxon[i_ali];
			i_msa++;
			i_ali++;
		} else if (msa[centerSequenceIdx][i_msa] == '-') {
			msa[taxonToAdd] += "-";
			i_msa++;
		} else if (aliSeqCenter[i_ali] == '-') {
			// add gap to all sequences already added to the msa, at position i
			for (size_t j = 0; j < msa.size(); ++j) {
				if (msa[j].empty())
					continue;
				if (j == taxonToAdd) {
					continue;
				}
				std::string left = msa[j].substr(0, i_msa);
				std::string right = msa[j].substr(i_msa, std::string::npos);
				msa[j] = left + "-" + right;
			}
			msa[taxonToAdd] += aliSeqNewTaxon[i_ali];
			i_ali++;
			i_msa++; // because we inserted the gap
		} else {
			std::cout << "msa[centerSequenceIdx][i_msa]: " << msa[centerSequenceIdx][i_msa] << "\n";
			std::cout << "aliSeqCenter[i_ali]: " << aliSeqCenter[i_ali] << "\n";
			throw std::runtime_error("This should not happen");
		}
	}

	// add the trailing gaps to the newly added taxon
	while (i_msa < msa[centerSequenceIdx].size()) {
		msa[taxonToAdd] += "-";
		i_msa++;
	}
	// add the trailing gaps/ remaining taxonToAdd sequence to the MSA
	while (i_ali < aliSeqNewTaxon.size()) {
		for (size_t j = 0; j < msa.size(); ++j) {
			if (msa[j].empty()) {
				continue;
			}
			if (j == taxonToAdd) {
				continue;
			}
			msa[j] += '-';
		}
		msa[taxonToAdd] += aliSeqNewTaxon[i_ali];
		i_ali++;
	}
}
