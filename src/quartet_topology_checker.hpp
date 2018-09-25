#pragma once

#ifdef WITH_GENESIS

#include <vector>
#include <array>
#include <string>
#include <unordered_map>

#include "../../genesis/lib/genesis/genesis.hpp"
#include "tree_information.hpp"
#include "indexed_concat.hpp"
#include "indexing/suffix_array_classic.hpp"

class TopologyChecker {
public:
	TopologyChecker(const IndexedConcatenatedSequence& concat, const std::string& speciesTreePath, const std::string& geneTreesPath,
			const std::string& multiSPAMPath);
	double getQICScore(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx, const std::array<size_t, 3>& counts);
	bool sameTopologyAsReference(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx, const std::array<size_t, 3>& counts);
	std::array<size_t, 3> findGeneTopologyCounts(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx);
	std::array<size_t, 3> findMultiSPAMTopologyCounts(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx);

	size_t pairwiseDistSum(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx);
	double quartetDifficulty(size_t uIdx, size_t vIdx, size_t wIdx, size_t zIdx);
	bool longBranchAttractionIssue(size_t uIdx, size_t vIdx, size_t wIdx, size_t zIdx);
private:
	size_t findReferenceTopology(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx);
	size_t findGeneTopology(size_t uIdx, size_t vIdx, size_t wIdx, size_t zIdx, TreeInformation& info);
	size_t rootIdx;
	std::vector<size_t> fastaIdxToRefIdx;
	TreeInformation informationReferenceTree;
	std::vector<TreeInformation> informationGeneTrees;
	std::unordered_map<std::string, size_t> taxonToReferenceID;

	std::vector<std::unordered_map<size_t, size_t> > refIDToGeneID;

	std::string textMultiSPAM;
	SuffixArrayClassic saMultiSPAM;
};

#endif
