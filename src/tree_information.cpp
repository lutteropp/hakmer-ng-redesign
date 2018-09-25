#ifdef WITH_GENESIS

#include "tree_information.hpp"
#include "../genesis/lib/genesis/genesis.hpp"

/**
 * Compute the distance in edges between two nodes, using their lowest common ancestor.
 * @param uIdx index of the first node in the tree
 * @param vIdx index of the second node in the tree
 */
unsigned TreeInformation::distanceInEdges(size_t uIdx, size_t vIdx) {
	size_t lcaIdx = lowestCommonAncestorIdx(uIdx, vIdx, myRootIndex);
	return dist_to_root[uIdx] + dist_to_root[vIdx] - 2 * dist_to_root[lcaIdx];
}

/**
 * Compute the distance in sum of branch lengths between two nodes, using their lowest common ancestor.
 * @param uIdx index of the first node in the tree
 * @param vIdx index of the second node in the tree
 */
double TreeInformation::distanceInBranchLengths(size_t uIdx, size_t vIdx) {
	size_t lcaIdx = lowestCommonAncestorIdx(uIdx, vIdx, myRootIndex);
	return dist_to_root_branches[uIdx] + dist_to_root_branches[vIdx] - 2 * dist_to_root_branches[lcaIdx];
}

/**
 * Given a quartet with topology ab|cd, return the sum of branch lengths connecting lca(a,b) to lca(c,d).
 * @param aIdx index of node a in the tree
 * @param bIdx index of node b in the tree
 * @param cIdx index of node c in the tree
 * @param dIdx index of node d in the tree
 */
double TreeInformation::connectingQuartetPathDistance(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx) {
	size_t lcaABR = lowestCommonAncestorIdx(aIdx, bIdx, myRootIndex);
	size_t lcaCDR = lowestCommonAncestorIdx(cIdx, dIdx, myRootIndex);
	if (lcaABR == lcaCDR) {
		size_t lca = lowestCommonAncestorIdx(aIdx, bIdx, lcaCDR);
		return distanceInBranchLengths(lca, lcaCDR);
	} else {
		size_t lca = lowestCommonAncestorIdx(cIdx, dIdx, lcaABR);
		return distanceInBranchLengths(lca, lcaABR);
	}
}

/**
 * Query the RangeMinimumQuery to find the index of the smallest entry between positions i and j in the euler tour levels.
 * This returns the ID of the lowest common ancestor of the nodes with IDs i and j, in respect to the tree's root node.
 * @param i the starting position
 * @param j the ending position
 */
size_t TreeInformation::rmqQueryCorrectOrder(size_t i, size_t j) {
	if (i <= j)
		return rrmq.query(i, j);
	else
		return rrmq.query(j, i);
}

/**
 * Return the ID of the root node in the tree.
 */
size_t TreeInformation::getRootIdx() {
	return myRootIndex;
}

/**
 * Return the ID of the lowest common ancestor of the nodes at uIdx and vIdx, using rootIdx as the root node ID.
 * @param uIdx ID of the node u
 * @param vIdx ID of the node v
 * @param rootIdx the ID of the node to be used as the root node
 */
size_t TreeInformation::lowestCommonAncestorIdx(size_t uIdx, size_t vIdx, size_t rootIdx) {
	size_t uEulerIdx = firstOccurrenceInEulerTour[uIdx];
	size_t vEulerIdx = firstOccurrenceInEulerTour[vIdx];
	size_t rootEulerIdx = firstOccurrenceInEulerTour[rootIdx];

	if (rootIdx == myRootIndex) {
		return eulerTourNodes[rmqQueryCorrectOrder(uEulerIdx, vEulerIdx)];
	} else { // take the "odd man out", see http://stackoverflow.com/questions/25371865/find-multiple-lcas-in-unrooted-tree
		size_t candidateOne = eulerTourNodes[rmqQueryCorrectOrder(uEulerIdx, vEulerIdx)];
		size_t candidateTwo = eulerTourNodes[rmqQueryCorrectOrder(uEulerIdx, rootEulerIdx)];
		size_t candidateThree = eulerTourNodes[rmqQueryCorrectOrder(vEulerIdx, rootEulerIdx)];
		if (candidateOne == candidateTwo) {
			return candidateThree;
		} else if (candidateOne == candidateThree) {
			return candidateTwo;
		} else {
			return candidateOne;
		}
	}
}

/**
 * @param tree the tree to build the TreeInformation for.
 */
void TreeInformation::init(genesis::tree::Tree const &tree) {
	dist_to_root = node_path_length_vector(tree);
	dist_to_root_branches = node_branch_length_distance_vector(tree);
	// Arrays for LCA computation
	std::vector<int> eulerTourLevels;

	firstOccurrenceInEulerTour.resize(tree.node_count());
	std::fill(firstOccurrenceInEulerTour.begin(), firstOccurrenceInEulerTour.end(), std::numeric_limits<size_t>::max());

	for (auto it : genesis::tree::eulertour(tree)) {
		eulerTourNodes.push_back(it.node().index());
		eulerTourLevels.push_back(dist_to_root[it.node().index()]);
		if (firstOccurrenceInEulerTour[it.node().index()] == std::numeric_limits<size_t>::max()) {
			firstOccurrenceInEulerTour[it.node().index()] = eulerTourNodes.size() - 1;
		}
	}
	myRootIndex = tree.root_node().index();

	rrmq = genesis::utils::RangeMinimumQuery(eulerTourLevels);
}

/**
 * @param tree the tree to build the TreeInformation for.
 */
TreeInformation::TreeInformation() {
	myRootIndex = 0;
}

#endif
