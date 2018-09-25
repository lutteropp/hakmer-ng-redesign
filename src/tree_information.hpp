#pragma once

#ifdef WITH_GENESIS

#include "../genesis/lib/genesis/genesis.hpp"
#include <vector>
#include <memory>
#include <algorithm>

/**
 * Helper class that computes the lowest common ancestor of two nodes in respect to any given root node,
 * as well as the distance between any two nodes in number of edges.
 */
class TreeInformation {
public:
	TreeInformation();
	void init(genesis::tree::Tree const& tree);
	/*
	 * @brief Returns the index of the lowest common ancestor of the nodes at index uIdx and vIdx with respect to root node at rootIdx.
	 */
	size_t lowestCommonAncestorIdx(size_t uIdx, size_t vIdx, size_t rootIdx);
	unsigned distanceInEdges(size_t uIdx, size_t vIdx);
	double distanceInBranchLengths(size_t uIdx, size_t vIdx);
	size_t getRootIdx();
	double connectingQuartetPathDistance(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx);
private:
	size_t myRootIndex; /**< index of the root node */
	genesis::utils::RangeMinimumQuery rrmq;
	size_t rmqQueryCorrectOrder(size_t i, size_t j);
	std::vector<size_t> eulerTourNodes; /**< the nodes of the tree visited in an Eulerian tour order */
	std::vector<size_t> firstOccurrenceInEulerTour; /**< the nodes' first occurrence in an Eulerian tour */
	std::vector<size_t> dist_to_root; /**< distance of a given node to the root in terms of number of edges */
	std::vector<double> dist_to_root_branches; /**< distance of a given node to the root in terms of sum of branch lengths */
};

#endif
