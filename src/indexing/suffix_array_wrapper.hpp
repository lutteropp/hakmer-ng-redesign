#pragma once

#include <stddef.h>
#include <string>

#include "../indexing/suffix_array_classic.hpp"
#include "../indexing/suffix_array_fm.hpp"
#include "../definitions.hpp"
#include "../options.hpp"

class SuffixArrayWrapper {
public:
	void setCompression(bool comp) {
		useCompression = comp;
	}
	bool getCompression() {
		return useCompression;
	}
	void buildSuffixArray(const std::string& seq, size_t nTotalSites, const Options& options) {
		if (useCompression) {
			saFM.buildSuffixArray(seq, options);
		} else {
			saClassic.buildSuffixArray(seq, nTotalSites, options);
		}
	}
	inline size_t operator[](size_t pos) const {
		if (useCompression) {
			return saFM[pos];
		} else {
			return saClassic[pos];
		}
	}
	SuffixArrayClassic& getSAClassic() {
		return saClassic;
	}
	SuffixArrayFM& getSAFM() {
		return saFM;
	}
private:
	SuffixArrayClassic saClassic;
	SuffixArrayFM saFM;
	bool useCompression;
};
