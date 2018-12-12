#include "approx_matching.hpp"
#include "../dna_functions.hpp"

size_t dpCoord(size_t s2Size, size_t i, size_t j) {
	return i * (s2Size + 1) + j;
}

int max3(int a, int b, int c) {
	return std::max(a, std::max(b, c));
}

std::vector<int> fillDPMatrixLocal(const std::string& s1, const std::string& s2) {
	const int match = 1;
	const int mismatch = -1;
	const int gap = -2;

	// create a DP table
	size_t nrows = s1.size() + 1;
	size_t ncols = s2.size() + 1;
	std::vector<int> matrix(nrows * ncols);
	// initialize it
	for (size_t i = 0; i < nrows; ++i) {
		matrix[dpCoord(s2.size(), i, 0)] = 0;
	}
	for (size_t j = 0; j < ncols; ++j) {
		matrix[dpCoord(s2.size(), 0, j)] = 0;
	}
	// fill it
	for (size_t i = 1; i < nrows; ++i) {
		for (size_t j = 1; j < ncols; ++j) {
			int score = (ambiguousMatch(s1[i], s2[i])) ? match : mismatch;
			matrix[dpCoord(s2.size(), i, j)] = max3(matrix[dpCoord(s2.size(), i - 1, j - 1)] + score,
					matrix[dpCoord(s2.size(), i - 1, j)] + gap, matrix[dpCoord(s2.size(), i, j - 1)] + gap);
			if (matrix[dpCoord(s2.size(), i, j)] < 0) {
				matrix[dpCoord(s2.size(), i, j)] = 0;
			}
		}
	}
	return matrix;
}

class ExtensionInfo {
public:
	size_t firstPosInText;
	size_t lastPosInText;
	size_t errors;
};

size_t traceBack(size_t i, size_t j, const std::string& pattern, const std::string& seq, const std::vector<int>& dpMatrix) {
	// check if the backtracking has ended:
	if (dpMatrix[dpCoord(seq.size(), i, j)] == 0) {
		return j;
	}

	// check for the largest neighboring entry
	int largestI = -1;
	int largestJ = -1;
	int largestVal = -1;

	// check diagonal neighbor (match/mismatch)
	if (i > 0 && j > 0) {
		int newI = i - 1;
		int newJ = j - 1;
		int newVal = dpMatrix[dpCoord(seq.size(), newI, newJ)];
		if (newVal > largestVal) {
			largestI = newI;
			largestJ = newJ;
			largestVal = newVal;
		}
	}

	// check left neighbor (indel)
	if (j > 0) {
		int newI = i;
		int newJ = j - 1;
		int newVal = dpMatrix[dpCoord(seq.size(), newI, newJ)];
		if (newVal > largestVal) {
			largestI = newI;
			largestJ = newJ;
			largestVal = newVal;
		}
	}

	// check top neighbor (indel)
	if (i > 0) {
		int newI = i - 1;
		int newJ = j;
		int newVal = dpMatrix[dpCoord(seq.size(), newI, newJ)];
		if (newVal > largestVal) {
			largestI = newI;
			largestJ = newJ;
			largestVal = newVal;
		}
	}

	return traceBack(largestI, largestJ, pattern, seq, dpMatrix);
}

inline std::vector<ExtensionInfo> leftExtensionSearch(const std::string& pattern, const std::string& seq, unsigned int maxErrors,
		unsigned int minErrors) {
	std::vector<ExtensionInfo> res;
	size_t rows = pattern.size() + 1;
	size_t cols = seq.size() + 1;
	std::vector<int> dpMatrix = fillDPMatrixLocal(pattern, seq);
	size_t lastRowIdx = rows - 1;

	// improvement: only return approximate matches with the best errors score
	int bestErrorsScore = maxErrors + 1;
	for (size_t j = 0; j < cols; ++j) {
		int errors = pattern.size() - dpMatrix[dpCoord(seq.size(), lastRowIdx, j)];
		if (errors < bestErrorsScore && errors >= (int) minErrors) {
			bestErrorsScore = errors;
		}
	}

	if (bestErrorsScore >= (int) minErrors && bestErrorsScore <= (int) maxErrors) {
		for (size_t j = 0; j < cols; ++j) {
			int errors = pattern.size() - dpMatrix[dpCoord(seq.size(), lastRowIdx, j)];
			if (errors == bestErrorsScore) {
				// found an occurrence, let's trace it back to the first zero we find
				ExtensionInfo info;
				info.lastPosInText = j + 1; // +1 because there is an extra epsilon character in the matrix
				info.errors = errors;
				info.firstPosInText = traceBack(lastRowIdx, j, pattern, seq, dpMatrix) + 1; // +1 because there is an extra epsilon character in the matrix
				res.push_back(info);
			}
		}
	}
	return res;
}

inline std::vector<ExtensionInfo> rightExtensionSearch(const std::string& pattern, const std::string& text, unsigned int maxErrors,
		unsigned int minErrors) {
	std::vector<ExtensionInfo> res;
	std::string patternReversed(pattern);
	std::reverse(patternReversed.begin(), patternReversed.end());
	std::string textReversed(pattern);
	std::reverse(textReversed.begin(), textReversed.end());

	std::vector<ExtensionInfo> leftExtendResult = leftExtensionSearch(patternReversed, textReversed, maxErrors, minErrors);

	for (size_t i = 0; i < leftExtendResult.size(); ++i) {
		unsigned int firstPosReversed = leftExtendResult[i].firstPosInText;
		unsigned int lastPosReversed = leftExtendResult[i].lastPosInText;

		ExtensionInfo info;
		info.firstPosInText = text.size() - lastPosReversed - 1;
		info.lastPosInText = text.size() - firstPosReversed - 1;
		info.errors = leftExtendResult[i].errors;
		res.push_back(info);
	}

	return res;
}

ApproximateMatcher::ApproximateMatcher(bool mismatchesOnly) :
		mismatchesOnly(mismatchesOnly) {
}

template<typename T>
void addAll(std::vector<T>& first, const std::vector<T>& second) {
	for (T t : second) {
		first.push_back(t);
	}
}

template<typename T>
void addAll(std::unordered_set<T, pair_hash>& first, const std::vector<T>& second) {
	for (T t : second) {
		first.insert(t);
	}
}

template<typename T>
void addAll(std::vector<T>& first, const std::unordered_set<T, pair_hash>& second) {
	for (T t : second) {
		first.push_back(t);
	}
}

std::string ApproximateMatcher::extractText(const std::string& seq, size_t firstPos, size_t lastPos) {
	if (lastPos >= seq.size()) {
		throw std::runtime_error("lastPos is too large");
	}
	return seq.substr(firstPos, lastPos - firstPos + 1);
}

std::vector<Occ> ApproximateMatcher::leftUntilExact(size_t jIdx, size_t lastPos, const std::vector<std::string>& subpatterns,
		const std::string& seq) {
	std::vector<Occ> res;
	if (jIdx <= 0) {
		return res;
	}
	std::string newPattern = subpatterns[jIdx - 1];
	if (lastPos < newPattern.size() + 1) {
		return res;
	}
	std::string text = extractText(seq, lastPos - newPattern.size() - 1, lastPos - 1);

	if (lastPos - newPattern.size() == 0) { // at beginning of whole text document
		if (ambiguousEqual(newPattern, text) && jIdx == 1) { // exact match
			Occ o;
			o.i = 0;
			o.begin = lastPos - newPattern.size();
			res.push_back(o);
		}
		return res;
	}
	if (text.size() < newPattern.size()) {
		return res; // extension is impossible.
	}

	// here comes the standard case, this is, we consider the text preceding P_jIdx; of size subpatterns[jIdx-1].size() + 1.
	assert(text.size() == newPattern.size() + 1);

	// CHECK FOR EXACT MATCH OR SINGLE SUBSTITUTION, ALSO ADD SINGLE INSERTION BEFORE MATCH
	std::string sameSizeTextRight = text.substr(1, std::string::npos);
	// check for exact match
	if (ambiguousEqual(sameSizeTextRight, newPattern)) { // found i
		Occ o;
		o.i = jIdx - 1;
		o.begin = lastPos - newPattern.size();
		res.push_back(o);
		return res;
	} else {
		// check for substitution, that is: P[0..i-1]*P[i+1..n]
		for (size_t i = 0; i < newPattern.size(); ++i) {
			bool valid = true;
			for (size_t j = 0; j < newPattern.size(); ++j) {
				if (j != i && !ambiguousMatch(newPattern[j], sameSizeTextRight[j])) {
					valid = false;
					break;
				}
			}
			if (valid) { // found single mismatch error
				std::vector<Occ> continued = leftUntilExact(jIdx - 1, lastPos - newPattern.size(), subpatterns, seq);
				for (Occ o : continued) {
					res.push_back(o);
				}
				//return res;
			}
		}
	}

	if (!mismatchesOnly) {
		// CHECK FOR SINGLE INSERTION AFTER PATTERN (P[0..n]*)
		std::string sameSizeTextLeft = text.substr(0, newPattern.size());
		if (ambiguousEqual(sameSizeTextLeft, newPattern)) { // found insertion of character after pattern, this is, P[0..n]*
			std::vector<Occ> continued = leftUntilExact(jIdx - 1, lastPos - newPattern.size() - 1, subpatterns, seq);
			for (Occ o : continued) {
				res.push_back(o);
			}
		}

		// CHECK FOR SINGLE INSERTION WITHIN PATTERN (P[0..i-1]*P[i..n])
		for (size_t i = 1; i < newPattern.size(); ++i) {
			if (ambiguousEqual(text.substr(0, i), newPattern.substr(0, i))
					&& ambiguousEqual(text.substr(i + 1, std::string::npos), newPattern.substr(i, std::string::npos))) {
				// found insertion error
				std::vector<Occ> continued = leftUntilExact(jIdx - 1, lastPos - newPattern.size() - 1, subpatterns, seq);
				for (Occ o : continued) {
					res.push_back(o);
				}
			}
		}

		// CHECK FOR SINGLE DELETION (P[0..i-1]P[i..n])
		std::string textOneShorter = text.substr(2, std::string::npos);
		assert(textOneShorter.size() == newPattern.size() - 1);
		for (size_t i = 0; i < newPattern.size(); ++i) {
			if (jIdx == 1) { // this means, newPattern is P1
				if (i == 0 || i == newPattern.size() - 1)
					continue; // avoid deletion of first character
			}
			if (ambiguousEqual(textOneShorter.substr(0, i), newPattern.substr(0, i))
					&& ambiguousEqual(textOneShorter.substr(i, std::string::npos), newPattern.substr(i + 1, std::string::npos))) {
				// found deletion error
				std::vector<Occ> continued = leftUntilExact(jIdx - 1, lastPos - newPattern.size() + 1, subpatterns, seq);
				for (Occ o : continued) {
					res.push_back(o);
				}
			}
		}
	}
	return res;
}

/*
 * left extend Pi...Pj to get to P1...Pj with at most i-1 errors in P1...P{i-1}.
 * @þaram pos The position of the first character of subpattern i in the text
 */
std::vector<ErrorOcc> ApproximateMatcher::leftExtend(size_t iIdx, const std::vector<std::string> &subpatterns, size_t pos,
		const std::string& seq) {
	std::vector<ErrorOcc> res;
	std::string searchedPattern = ""; // This string has to be found with at most iIdx errors (because iIdx starts from zero, not from one).
	for (size_t i = 0; i < iIdx; ++i) {
		searchedPattern = searchedPattern + subpatterns[i];
	}
	if (searchedPattern.size() > pos) { // bail, because the pattern cannot be found
		return res;
	}
	if (searchedPattern.empty()) { // P1 has already been found
		ErrorOcc occ;
		occ.pos = pos - searchedPattern.size();
		occ.errors = 0;
		res.push_back(occ);
		return res;
	}
	if (!mismatchesOnly) {
		size_t startPos = 0;
		if (searchedPattern.size() + iIdx < pos) {
			startPos = pos - searchedPattern.size() - iIdx; // -iIdx because of insertion errors
		}
		std::string text = extractText(seq, startPos, pos - 1); // -iIdx because of insertion errors

		unsigned int minErrors = 0;
		unsigned int maxErrors = iIdx;
		std::vector<ExtensionInfo> leftExtensions = leftExtensionSearch(searchedPattern, text, maxErrors, minErrors);
		for (size_t i = 0; i < leftExtensions.size(); ++i) {
			ErrorOcc occ;
			occ.pos = startPos + leftExtensions[i].firstPosInText;
			occ.errors = leftExtensions[i].errors;
			res.push_back(occ);
		}
	} else {
		std::string text = extractText(seq, pos - searchedPattern.size(), pos - 1);
		size_t mis = 0;
		assert(text.size() == searchedPattern.size());
		for (size_t i = 0; i < text.size(); ++i) {
			if (!ambiguousMatch(text[i], searchedPattern[i])) {
				mis++;
				if (mis > iIdx) {
					break;
				}
			}
		}
		if (mis <= iIdx) {
			ErrorOcc occ;
			occ.pos = pos - searchedPattern.size();
			occ.errors = mis;
			res.push_back(occ);
		}
	}
	return res;
}

/*
 * right extend P1...Pj to get to P1...Pn with at most maxErrors errors in P{j+1}...Pn.
 */
std::vector<std::pair<size_t, size_t> > ApproximateMatcher::rightExtend(size_t jIdx, const std::vector<std::string> &subpatterns,
		size_t posJ, size_t maxErrors, size_t p1Pos, size_t minErrors, const std::string& seq) {
	std::vector<std::pair<size_t, size_t> > res;
	std::string searchedPattern = ""; // This string has to be found with at most maxErrors errors.
	for (size_t i = jIdx + 1; i < subpatterns.size(); ++i) {
		searchedPattern = searchedPattern + subpatterns[i];
	}
	if (searchedPattern.empty() && minErrors == 0) { // pattern has already been found
		if (posJ + subpatterns[jIdx].size() - 1 < seq.size()) {
			res.push_back(std::make_pair(p1Pos, posJ + subpatterns[jIdx].size() - 1));
		}
		return res;
	}
	if (!mismatchesOnly) {
		// isn't this just some edit distance thing?
		size_t lastTextPos = posJ + subpatterns[jIdx].size() + searchedPattern.size() - 1 + maxErrors;
		if (lastTextPos >= seq.size() && lastTextPos - maxErrors < seq.size()) {
			lastTextPos = seq.size() - 1;
		}

		if (lastTextPos < seq.size()) {
			size_t firstPos = posJ + subpatterns[jIdx].size();
			std::string text = extractText(seq, firstPos, lastTextPos);

			std::vector<ExtensionInfo> rightExtensions = rightExtensionSearch(searchedPattern, text, maxErrors, minErrors);
			for (size_t i = 0; i < rightExtensions.size(); ++i) {
				if (firstPos + rightExtensions[i].lastPosInText < seq.size()) {
					res.push_back(std::make_pair(p1Pos, firstPos + rightExtensions[i].lastPosInText));
				}
			}
		}
		return res;
	} else {
		size_t lastTextPos = posJ + subpatterns[jIdx].size() + searchedPattern.size() - 1;
		if (lastTextPos < seq.size()) {
			std::string text = extractText(seq, posJ + subpatterns[jIdx].size(), lastTextPos);
			assert(text.size() == searchedPattern.size());
			size_t mis = 0;
			for (size_t i = 0; i < text.size(); ++i) {
				if (!ambiguousMatch(text[i], searchedPattern[i])) {
					mis++;
					if (mis > maxErrors) {
						break;
					}
				}
			}
			if (mis <= maxErrors && mis >= minErrors) {
				if (posJ + subpatterns[jIdx].size() + searchedPattern.size() - 1 < seq.size()) {
					res.push_back(std::make_pair(p1Pos, posJ + subpatterns[jIdx].size() + searchedPattern.size() - 1));
				}
			}
		}
		return res;
	}
}

/*
 * Check if the potential occurrence area is already taken by another k-mer block and thus would be discarded
 */
bool alreadyTaken(const PresenceChecker& checker, size_t startPos, size_t size) {
	bool res = false;
	for (size_t i = startPos; i < startPos + size; ++i) {
		if (!checker.isFree(i)) {
			res = true;
			break;
		}
	}
	return res;
}

// Better version, as of Nov 10, 2015, handles boundary cases
/**
 * Perform a binary search in the suffix array to find matches.
 * @param patternSeqStartPos Index of the first base of the pattern in question in the input sequence data.
 * @param m Length of the pattern in question in the input sequence data.
 * @param kmerExtended3Prime Result giving the first and last index of the pattern occurrence in the suffix array as well as the total count.
 * @param S The input sequence data.
 */
bool binarySearch3Prime(const std::string& pattern, std::pair<size_t, size_t>& res, const std::string& text,
		const std::vector<size_t>& SA) {
	size_t l, h, mid, pos;
	size_t m = pattern.size();
	l = 0;
	h = text.size() - 1;
	while (l < h) { // find lower bound
		mid = l + (h - l) / 2; // calculate it this way to avoid overflows instead of (l+h)/2
		pos = SA[mid];
		if (pos > text.size() - m - 1) { /*printf ("Seq array out of bounds in bsearch; bailing\n");*/
			return false;
		}
		if (pattern.compare(0, m, text, pos, m) <= 0 || ambiguousEqual(text.substr(pos, m), pattern.substr(0, m))) {
			h = mid;
		} else {
			l = mid + 1;
		}
	}
	pos = SA[l];
	if (pos > text.size() - m - 1) { /*printf ("Seq array out of bounds in bsearch; bailing\n");*/
		return false;
	}

	if (!ambiguousEqual(text.substr(pos, m), pattern.substr(0, m))) {
		return false;
	} else {
		res.first = l;
	}

	l = 0;
	h = text.size() - 1;
	while (l < h) { // find upper bound
		if (h - l == 1) {
			mid = h; // treats integer rounding in mid calc for boundary case when interval is just [a,a+1]
		} else {
			mid = l + (h - l) / 2; // calculate it this way to avoid overflows instead of (l+h)/2
		}
		pos = SA[mid];
		if (pos > text.size() - m - 1) { /*printf ("Seq array out of bounds in bsearch; bailing\n");*/
			return false;
		}
		if (pattern.compare(0, m, text, pos, m) >= 0 || ambiguousEqual(text.substr(pos, m), pattern.substr(0, m))) {
			l = mid;
		} else {
			h = mid - 1;
		}
	}
	pos = SA[l];
	if (pos > text.size() - m - 1) { /*printf ("Seq array out of bounds in bsearch; bailing\n");*/
		return false;
	}
	if (!ambiguousEqual(text.substr(pos, m), pattern.substr(0, m))) {
		return false;
	} else {
		res.second = l - 1;
	}
	return true;
}

size_t exactMatches(const std::string& pattern, const std::string& text, std::vector<size_t> &matches, const std::vector<size_t>& SA) {
	std::pair<size_t, size_t> bl;
	binarySearch3Prime(pattern, bl, text, SA);
	for (size_t i = bl.first; i <= bl.second; ++i) {
		matches.push_back(SA[i]);
	}
	return matches.size();
}

std::vector<std::pair<size_t, size_t> > ApproximateMatcher::findOccurrences(const std::string& seq, const std::vector<size_t>& SA,
		PresenceChecker& checker, const std::string& pattern, size_t maxErrors, size_t minErrors, bool keepOverlaps) {
	std::vector<std::pair<size_t, size_t> > result;
	std::unordered_set<std::pair<size_t, size_t>, pair_hash> res;
// split the pattern into maxMismatches + 2 parts
	std::vector<std::string> subpatterns;
	subpatterns.resize(maxErrors + 2);
	for (size_t i = 0; i < maxErrors + 2; ++i) {
		size_t beginPos = pattern.size() * i / (maxErrors + 2);
		size_t endPos = pattern.size() * (i + 1) / (maxErrors + 2);
		subpatterns[i] = pattern.substr(beginPos, endPos - beginPos);
	}

	for (size_t j = 1; j < maxErrors + 2; ++j) { // iterate over possible parts for P_j
		std::vector<size_t> positionsPj;
		exactMatches(subpatterns[j], seq, positionsPj, SA);
		// search parts preceding P_j with at most one error, until a part is found exactly (this is P_i) -> via backtracking!
		// this gives us occurrences of P_i...P_j with exactly j-i-1 errors
//#pragma omp parallel for schedule(dynamic)
		for (size_t posJIdx = 0; posJIdx < positionsPj.size(); ++posJIdx) {
			size_t posJ = positionsPj[posJIdx];
			// optimization: stop the search if the position is already taken by some other occurrence
			if (alreadyTaken(checker, posJ, subpatterns[j].size())) {
				continue;
			}
			std::vector<Occ> occsPi = leftUntilExact(j, posJ, subpatterns, seq);
			// do left-extension to find P_1...P_{i-1} with at most i-1 errors
			// this gives us the occurrences of P_1...P_j with e errors
			for (Occ occI : occsPi) {
				size_t errors = j - occI.i - 1;
				std::vector<ErrorOcc> occsP1Pj = leftExtend(occI.i, subpatterns, occI.begin, seq);
				// do right-extension: e errors in P_1...P_j -> <= k-e errors in P_{j+1}...P_{k+2}
				for (ErrorOcc e : occsP1Pj) {
					errors += e.errors;

					std::vector<std::pair<size_t, size_t> > occsFullPattern = rightExtend(j, subpatterns, posJ, maxErrors - e.errors, e.pos,
							minErrors - errors, seq);
//#pragma omp critical
					addAll(res, occsFullPattern);
				}
			}
		}
	}

	addAll(result, res);
	std::sort(result.begin(), result.end());

	return result;
}
