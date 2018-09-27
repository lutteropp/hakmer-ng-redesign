#ifdef WITH_GENESIS

#include "quartet_topology_checker.hpp"
#include <fstream>

TopologyChecker::TopologyChecker() {
}

void TopologyChecker::init(const IndexedConcatenatedSequence& concat, const std::string& speciesTreePath,
		const std::string& geneTreesPath, const std::string& multiSPAMPath) {
	// read the species tree; including Sebastians quickfix
	genesis::tree::Tree speciesTree = genesis::tree::DefaultTreeNewickReader().from_file(speciesTreePath);
	if (speciesTree.root_node().degree() == 1) {
		std::string newick = genesis::tree::DefaultTreeNewickWriter().to_string(speciesTree);

		int c = 0;
		size_t a = 0;
		size_t b = newick.size();
		size_t d = newick.size();
		for (size_t i = 0; i < newick.size() and b == newick.size(); ++i) {
			if (newick[i] == '(')
			c++;
			if (newick[i] == ')')
			c--;
			if (c == 2 and a == 0)
			a = i;
			if (c == 1 and a > 0) {
				b = i;
				d = b;
				if (newick[i + 1] == ':') {
					d++;
					while (newick[d] != ',' and newick[d] != ')')
					d++;
				}
			}
		}
		std::string newick2 = newick.substr(0, a) + newick.substr(a + 1, b - a - 1)
		+ newick.substr(d, newick.size() - d);
		speciesTree = genesis::tree::DefaultTreeNewickReader().from_string(newick2);
	}
	rootIdx = speciesTree.root_node().index();

	fastaIdxToRefIdx.resize(concat.nTax());
	for (size_t i = 0; i < speciesTree.node_count(); ++i) {
		if (speciesTree.node_at(i).is_leaf()) {
			std::string leafName = speciesTree.node_at(i).data<genesis::tree::DefaultNodeData>().name;
			taxonToReferenceID[leafName] = i;
			bool found = false;
			const std::vector<std::string> taxLabels = concat.getTaxonLabels();
			for (size_t j = 0; j < taxLabels.size(); ++j) {
				if (taxLabels[j] == leafName) {
					fastaIdxToRefIdx[j] = i;
					found = true;
					break;
				}
			}
			if (!found) {
				throw std::runtime_error("Could not find any taxon named " + leafName + " in the FASTA file");
			}
		}
	}

	informationReferenceTree.init(speciesTree);

	if (!geneTreesPath.empty()) {
		genesis::utils::InputStream instream(
				genesis::utils::make_unique<genesis::utils::FileInputSource>(geneTreesPath));
		auto itTree = genesis::tree::NewickInputIterator(instream, genesis::tree::DefaultTreeNewickReader());
		size_t i = 0;
		while (itTree) { // iterate over the set of evaluation trees
			genesis::tree::Tree const& tree = *itTree;
			std::unordered_map < size_t, size_t > refIDToGeneIDCurr;
			for (size_t i = 0; i < tree.node_count(); ++i) {
				if (tree.node_at(i).is_leaf()) {
					std::string leafName = tree.node_at(i).data<genesis::tree::DefaultNodeData>().name;
					if (taxonToReferenceID.find(leafName) != taxonToReferenceID.end()) {
						refIDToGeneIDCurr[taxonToReferenceID[leafName]] = i;
					}
				}
			}
			refIDToGeneID.push_back(refIDToGeneIDCurr);
			informationGeneTrees.push_back(TreeInformation());
			informationGeneTrees[informationGeneTrees.size() - 1].init(tree);
			++itTree;
		}
	}

	if (!multiSPAMPath.empty()) {
		std::ifstream infile(multiSPAMPath);
		std::string text {istreambuf_iterator<char>(infile), istreambuf_iterator<char>()};
		textMultiSPAM = text;
		saMultiSPAM.buildSuffixArray(textMultiSPAM);
	}
}

size_t TopologyChecker::pairwiseDistSum(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx) {
	size_t res = 0;
	size_t uIdx = fastaIdxToRefIdx[aIdx];
	size_t vIdx = fastaIdxToRefIdx[bIdx];
	size_t wIdx = fastaIdxToRefIdx[cIdx];
	size_t zIdx = fastaIdxToRefIdx[dIdx];
	res += informationReferenceTree.distanceInEdges(uIdx, vIdx);
	res += informationReferenceTree.distanceInEdges(uIdx, wIdx);
	res += informationReferenceTree.distanceInEdges(uIdx, zIdx);
	res += informationReferenceTree.distanceInEdges(vIdx, wIdx);
	res += informationReferenceTree.distanceInEdges(vIdx, zIdx);
	res += informationReferenceTree.distanceInEdges(wIdx, zIdx);
	return res;
}

bool TopologyChecker::longBranchAttractionIssue(size_t uIdx, size_t vIdx, size_t wIdx, size_t zIdx) {
	size_t aIdx = fastaIdxToRefIdx[uIdx];
	size_t bIdx = fastaIdxToRefIdx[vIdx];
	size_t cIdx = fastaIdxToRefIdx[wIdx];
	size_t dIdx = fastaIdxToRefIdx[zIdx];
	QuartetTopology refTopo = findReferenceTopology(uIdx, vIdx, wIdx, zIdx);
	size_t a, b, c, d;
	if (refTopo == QuartetTopology::AB_CD) {
		a = aIdx;
		b = bIdx;
		c = cIdx;
		d = dIdx;
	} else if (refTopo == QuartetTopology::AC_BD) {
		a = aIdx;
		b = cIdx;
		c = bIdx;
		d = dIdx;
	} else if (refTopo == QuartetTopology::AD_BC) {
		a = aIdx;
		b = dIdx;
		c = bIdx;
		d = cIdx;
	} else {
		throw std::runtime_error("Quartet has star topology");
	}

	double ab = informationReferenceTree.distanceInBranchLengths(a, b);
	double ac = informationReferenceTree.distanceInBranchLengths(a, c);
	double ad = informationReferenceTree.distanceInBranchLengths(a, d);
	double bc = informationReferenceTree.distanceInBranchLengths(b, c);
	double bd = informationReferenceTree.distanceInBranchLengths(b, d);
	double cd = informationReferenceTree.distanceInBranchLengths(c, d);

	std::vector<double> dists = {ab, ac, ad, bc, bd, cd};
	double minDist = ab;
	for (size_t i = 0; i < dists.size(); ++i) {
		minDist = std::min(minDist, dists[i]);
	}

	if (minDist == ab || minDist == cd) {
		return false;
	} else {
		return true;
	}
}

double TopologyChecker::quartetDifficulty(size_t uIdx, size_t vIdx, size_t wIdx, size_t zIdx) {
	size_t aIdx = fastaIdxToRefIdx[uIdx];
	size_t bIdx = fastaIdxToRefIdx[vIdx];
	size_t cIdx = fastaIdxToRefIdx[wIdx];
	size_t dIdx = fastaIdxToRefIdx[zIdx];

	QuartetTopology refTopo = findReferenceTopology(uIdx, vIdx, wIdx, zIdx);

	size_t a, b, c, d;
	if (refTopo == QuartetTopology::AB_CD) { // topology ab|cd
		a = aIdx;
		b = bIdx;
		c = cIdx;
		d = dIdx;
	} else if (refTopo == QuartetTopology::AC_BD) { // topology ac|bd
		a = aIdx;
		b = cIdx;
		c = bIdx;
		d = dIdx;
	} else if (refTopo == QuartetTopology::AD_BC) { // topology ad|bc
		a = aIdx;
		b = dIdx;
		c = bIdx;
		d = cIdx;
	} else {
		throw std::runtime_error("Cannot compute difficulty because the quartet has star topology");
	}

	// the larger the return value, the easier the quartet topology inference

	double minDist = std::min(informationReferenceTree.distanceInBranchLengths(a, b),
			informationReferenceTree.distanceInBranchLengths(c, d));

	return informationReferenceTree.connectingQuartetPathDistance(a, b, c, d) / minDist;
}

QuartetTopology TopologyChecker::findReferenceTopology(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx) {
	size_t uIdx = fastaIdxToRefIdx[aIdx];
	size_t vIdx = fastaIdxToRefIdx[bIdx];
	size_t wIdx = fastaIdxToRefIdx[cIdx];
	size_t zIdx = fastaIdxToRefIdx[dIdx];

	size_t lca_uv = informationReferenceTree.lowestCommonAncestorIdx(uIdx, vIdx, rootIdx);
	size_t lca_uw = informationReferenceTree.lowestCommonAncestorIdx(uIdx, wIdx, rootIdx);
	size_t lca_uz = informationReferenceTree.lowestCommonAncestorIdx(uIdx, zIdx, rootIdx);
	size_t lca_vw = informationReferenceTree.lowestCommonAncestorIdx(vIdx, wIdx, rootIdx);
	size_t lca_vz = informationReferenceTree.lowestCommonAncestorIdx(vIdx, zIdx, rootIdx);
	size_t lca_wz = informationReferenceTree.lowestCommonAncestorIdx(wIdx, zIdx, rootIdx);

	if (informationReferenceTree.distanceInEdges(lca_uv, lca_wz)
			> informationReferenceTree.distanceInEdges(lca_uw, lca_vz)
			&& informationReferenceTree.distanceInEdges(lca_uv, lca_wz)
			> informationReferenceTree.distanceInEdges(lca_uz, lca_vw)) {
		return QuartetTopology::AB_CD;
	} else if (informationReferenceTree.distanceInEdges(lca_uw, lca_vz)
			> informationReferenceTree.distanceInEdges(lca_uv, lca_wz)
			&& informationReferenceTree.distanceInEdges(lca_uw, lca_vz)
			> informationReferenceTree.distanceInEdges(lca_uz, lca_vw)) {
		return QuartetTopology::AC_BD;
	} else if (informationReferenceTree.distanceInEdges(lca_uz, lca_vw)
			> informationReferenceTree.distanceInEdges(lca_uv, lca_wz)
			&& informationReferenceTree.distanceInEdges(lca_uz, lca_vw)
			> informationReferenceTree.distanceInEdges(lca_uw, lca_vz)) {
		return QuartetTopology::AD_BC;
	} else {
		return QuartetTopology::STAR;
	}
}

QuartetTopology TopologyChecker::findGeneTopology(size_t uIdx, size_t vIdx, size_t wIdx, size_t zIdx, TreeInformation& info) {
	size_t lca_uv = info.lowestCommonAncestorIdx(uIdx, vIdx, rootIdx);
	size_t lca_uw = info.lowestCommonAncestorIdx(uIdx, wIdx, rootIdx);
	size_t lca_uz = info.lowestCommonAncestorIdx(uIdx, zIdx, rootIdx);
	size_t lca_vw = info.lowestCommonAncestorIdx(vIdx, wIdx, rootIdx);
	size_t lca_vz = info.lowestCommonAncestorIdx(vIdx, zIdx, rootIdx);
	size_t lca_wz = info.lowestCommonAncestorIdx(wIdx, zIdx, rootIdx);

	if (info.distanceInEdges(lca_uv, lca_wz) > info.distanceInEdges(lca_uw, lca_vz)
			&& info.distanceInEdges(lca_uv, lca_wz) > info.distanceInEdges(lca_uz, lca_vw)) {
		return QuartetTopology::AB_CD;
	} else if (info.distanceInEdges(lca_uw, lca_vz) > info.distanceInEdges(lca_uv, lca_wz)
			&& info.distanceInEdges(lca_uw, lca_vz) > info.distanceInEdges(lca_uz, lca_vw)) {
		return QuartetTopology::AC_BD;
	} else if (info.distanceInEdges(lca_uz, lca_vw) > info.distanceInEdges(lca_uv, lca_wz)
			&& info.distanceInEdges(lca_uz, lca_vw) > info.distanceInEdges(lca_uw, lca_vz)) {
		return QuartetTopology::AD_BC;
	} else {
		return QuartetTopology::STAR;
	}
}

bool TopologyChecker::sameTopologyAsReference(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx,
		QuartetTopology topo) {
	QuartetTopology refTopology = findReferenceTopology(aIdx, bIdx, cIdx, dIdx);
	return (topo == refTopology);
}

double TopologyChecker::getQICScore(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx,
		const std::array<size_t, 3>& counts) {
	if (counts[0] + counts[1] + counts[2] == 0) {
		return 0.0;
	}

	QuartetTopology refTopology = findReferenceTopology(aIdx, bIdx, cIdx, dIdx);
	size_t q1, q2, q3;
	if (refTopology == QuartetTopology::AB_CD) { // ab|cd
		q1 = 0;
		q2 = 1;
		q3 = 2;
	} else if (refTopology == QuartetTopology::AC_BD) { // ac|bd
		q1 = 1;
		q2 = 0;
		q3 = 2;
	} else if (refTopology == QuartetTopology::AD_BC) { // ad|bc
		q1 = 2;
		q2 = 1;
		q3 = 0;
	} else {
		throw std::runtime_error("The quartet has star-topology in the reference tree");
	}

	size_t sum = counts[q1] + counts[q2] + counts[q3];
	double p_q1 = (double) counts[q1] / sum;
	double p_q2 = (double) counts[q2] / sum;
	double p_q3 = (double) counts[q3] / sum;

	// qic = 1 + p_q1 * log(p_q1) + p_q2 * log(p_q2) + p_q3 * log(p_q3);
	double qic = 1;
	if (p_q1 != 0)
	qic += p_q1 * log(p_q1) / log(3);
	if (p_q2 != 0)
	qic += p_q2 * log(p_q2) / log(3);
	if (p_q3 != 0)
	qic += p_q3 * log(p_q3) / log(3);

	if (counts[q1] < counts[q2] || counts[q1] < counts[q3]) {
		return qic * -1;
	} else {
		return qic;
	}
}

std::array<size_t, 3> TopologyChecker::findGeneTopologyCounts(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx) {
	std::array < size_t, 3> counts = {0,0,0};

	for (size_t i = 0; i < informationGeneTrees.size(); ++i) {
		if (refIDToGeneID[i].find(fastaIdxToRefIdx[aIdx]) == refIDToGeneID[i].end() || refIDToGeneID[i].find(fastaIdxToRefIdx[bIdx]) == refIDToGeneID[i].end() || refIDToGeneID[i].find(fastaIdxToRefIdx[cIdx]) == refIDToGeneID[i].end() || refIDToGeneID[i].find(fastaIdxToRefIdx[dIdx]) == refIDToGeneID[i].end()) {
			continue;
		}

		size_t uIdx = refIDToGeneID[i][fastaIdxToRefIdx[aIdx]];
		size_t vIdx = refIDToGeneID[i][fastaIdxToRefIdx[bIdx]];
		size_t wIdx = refIDToGeneID[i][fastaIdxToRefIdx[cIdx]];
		size_t zIdx = refIDToGeneID[i][fastaIdxToRefIdx[dIdx]];

		QuartetTopology topo = findGeneTopology(uIdx, vIdx, wIdx, zIdx, informationGeneTrees[i]);
		if (topo == QuartetTopology::AB_CD) {
			counts[0]++;
		} else if (topo == QuartetTopology::AC_BD) {
			counts[1]++;
		} else if (topo == QuartetTopology::AD_BC) {
			counts[2]++;
		}
	}
	return counts;
}

std::array<size_t, 3> TopologyChecker::findMultiSPAMTopologyCounts(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx) {
	std::array < size_t, 3> counts = {0,0,0};
	std::string abcd = std::to_string(aIdx) + "," + std::to_string(bIdx) + "|" + std::to_string(cIdx) + "," + std::to_string(dIdx) + ":";
	std::string abdc = std::to_string(aIdx) + "," + std::to_string(bIdx) + "|" + std::to_string(dIdx) + "," + std::to_string(cIdx) + ":";
	std::string acbd = std::to_string(aIdx) + "," + std::to_string(cIdx) + "|" + std::to_string(bIdx) + "," + std::to_string(dIdx) + ":";
	std::string acdb = std::to_string(aIdx) + "," + std::to_string(cIdx) + "|" + std::to_string(dIdx) + "," + std::to_string(bIdx) + ":";
	std::string adbc = std::to_string(aIdx) + "," + std::to_string(dIdx) + "|" + std::to_string(bIdx) + "," + std::to_string(cIdx) + ":";
	std::string adcb = std::to_string(aIdx) + "," + std::to_string(dIdx) + "|" + std::to_string(cIdx) + "," + std::to_string(bIdx) + ":";

	std::string bacd = std::to_string(bIdx) + "," + std::to_string(aIdx) + "|" + std::to_string(cIdx) + "," + std::to_string(dIdx) + ":";
	std::string badc = std::to_string(bIdx) + "," + std::to_string(aIdx) + "|" + std::to_string(dIdx) + "," + std::to_string(cIdx) + ":";
	std::string bcad = std::to_string(bIdx) + "," + std::to_string(cIdx) + "|" + std::to_string(aIdx) + "," + std::to_string(dIdx) + ":";
	std::string bcda = std::to_string(bIdx) + "," + std::to_string(cIdx) + "|" + std::to_string(dIdx) + "," + std::to_string(aIdx) + ":";
	std::string bdac = std::to_string(bIdx) + "," + std::to_string(dIdx) + "|" + std::to_string(aIdx) + "," + std::to_string(cIdx) + ":";
	std::string bdca = std::to_string(bIdx) + "," + std::to_string(dIdx) + "|" + std::to_string(cIdx) + "," + std::to_string(aIdx) + ":";

	std::string cabd = std::to_string(cIdx) + "," + std::to_string(aIdx) + "|" + std::to_string(bIdx) + "," + std::to_string(dIdx) + ":";
	std::string cadb = std::to_string(cIdx) + "," + std::to_string(aIdx) + "|" + std::to_string(dIdx) + "," + std::to_string(bIdx) + ":";
	std::string cbad = std::to_string(cIdx) + "," + std::to_string(bIdx) + "|" + std::to_string(aIdx) + "," + std::to_string(dIdx) + ":";
	std::string cbda = std::to_string(cIdx) + "," + std::to_string(bIdx) + "|" + std::to_string(dIdx) + "," + std::to_string(aIdx) + ":";
	std::string cdab = std::to_string(cIdx) + "," + std::to_string(dIdx) + "|" + std::to_string(aIdx) + "," + std::to_string(bIdx) + ":";
	std::string cdba = std::to_string(cIdx) + "," + std::to_string(dIdx) + "|" + std::to_string(bIdx) + "," + std::to_string(aIdx) + ":";

	std::string dabc = std::to_string(dIdx) + "," + std::to_string(aIdx) + "|" + std::to_string(bIdx) + "," + std::to_string(cIdx) + ":";
	std::string dacb = std::to_string(dIdx) + "," + std::to_string(aIdx) + "|" + std::to_string(cIdx) + "," + std::to_string(bIdx) + ":";
	std::string dbac = std::to_string(dIdx) + "," + std::to_string(bIdx) + "|" + std::to_string(aIdx) + "," + std::to_string(cIdx) + ":";
	std::string dbca = std::to_string(dIdx) + "," + std::to_string(bIdx) + "|" + std::to_string(cIdx) + "," + std::to_string(aIdx) + ":";
	std::string dcab = std::to_string(dIdx) + "," + std::to_string(cIdx) + "|" + std::to_string(aIdx) + "," + std::to_string(bIdx) + ":";
	std::string dcba = std::to_string(dIdx) + "," + std::to_string(cIdx) + "|" + std::to_string(bIdx) + "," + std::to_string(aIdx) + ":";

	// ab|cd topology:
	counts[0] += saMultiSPAM.countExactMatches(abcd, textMultiSPAM);
	counts[0] += saMultiSPAM.countExactMatches(abdc, textMultiSPAM);
	counts[0] += saMultiSPAM.countExactMatches(bacd, textMultiSPAM);
	counts[0] += saMultiSPAM.countExactMatches(badc, textMultiSPAM);
	counts[0] += saMultiSPAM.countExactMatches(cdab, textMultiSPAM);
	counts[0] += saMultiSPAM.countExactMatches(cdba, textMultiSPAM);
	counts[0] += saMultiSPAM.countExactMatches(dcab, textMultiSPAM);
	counts[0] += saMultiSPAM.countExactMatches(dcba, textMultiSPAM);
	// ac|bd topology:
	counts[1] += saMultiSPAM.countExactMatches(acbd, textMultiSPAM);
	counts[1] += saMultiSPAM.countExactMatches(acdb, textMultiSPAM);
	counts[1] += saMultiSPAM.countExactMatches(cabd, textMultiSPAM);
	counts[1] += saMultiSPAM.countExactMatches(cadb, textMultiSPAM);
	counts[1] += saMultiSPAM.countExactMatches(bdac, textMultiSPAM);
	counts[1] += saMultiSPAM.countExactMatches(bdca, textMultiSPAM);
	counts[1] += saMultiSPAM.countExactMatches(dbac, textMultiSPAM);
	counts[1] += saMultiSPAM.countExactMatches(dbca, textMultiSPAM);
	// ad|bc topology:
	counts[2] += saMultiSPAM.countExactMatches(adbc, textMultiSPAM);
	counts[2] += saMultiSPAM.countExactMatches(adcb, textMultiSPAM);
	counts[2] += saMultiSPAM.countExactMatches(dabc, textMultiSPAM);
	counts[2] += saMultiSPAM.countExactMatches(dacb, textMultiSPAM);
	counts[2] += saMultiSPAM.countExactMatches(bcad, textMultiSPAM);
	counts[2] += saMultiSPAM.countExactMatches(bcda, textMultiSPAM);
	counts[2] += saMultiSPAM.countExactMatches(cbda, textMultiSPAM);
	counts[2] += saMultiSPAM.countExactMatches(cbad, textMultiSPAM);

	return counts;
}

#endif
