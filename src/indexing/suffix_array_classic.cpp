#include "../indexing/suffix_array_classic.hpp"
#include <stdexcept>

const std::vector<size_t>& SuffixArrayClassic::getSA() const {
	return SA;
}

const std::vector<size_t>& SuffixArrayClassic::getLCP() const {
	return lcp;
}
