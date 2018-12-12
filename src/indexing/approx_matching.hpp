#pragma once

#include <stddef.h>
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <cassert>
#include <unordered_set>

#include "suffix_array_classic.hpp"
#include "../presence_checker.hpp"

/**
 * Do approximate string matching with 01*0 lossless seeds by C.Vroland et al., 2016
 */

struct pair_hash {
	inline std::size_t operator()(const std::pair<size_t, size_t> & v) const {
		return v.first * 31 + v.second;
	}
};

struct Occ {
	size_t i;
	size_t begin;
};

struct ErrorOcc {
	size_t pos;
	size_t errors;
};

class ApproximateMatcher {
public:
	ApproximateMatcher(bool mismatchesOnly);
	std::vector<std::pair<size_t, size_t> > findFewOccurrences(const std::string& seq, const std::vector<size_t>& SA, const std::vector<size_t>& lcp, PresenceChecker& checker,
			const std::string& pattern, size_t maxErrors, size_t minErrors, bool keepOverlaps, const IndexedConcatenatedSequence& concat, std::vector<size_t>& taxPresence);
private:
	std::string extractText(const std::string& seq, size_t firstPos, size_t lastPos);
	std::vector<Occ> leftUntilExact(size_t jIdx, size_t lastPos, const std::vector<std::string>& subpatterns, const std::string& seq);
	std::vector<ErrorOcc> leftExtend(size_t iIdx, const std::vector<std::string> &subpatterns, size_t pos, const std::string& seq);
	std::vector<std::pair<size_t, size_t> > rightExtend(size_t jIdx, const std::vector<std::string> &subpatterns, size_t posJ,
			size_t maxErrors, size_t p1Pos, size_t minErrors, const std::string& seq);
	bool mismatchesOnly;
};
