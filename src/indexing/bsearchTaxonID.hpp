#pragma once

#include "suffix_array_wrapper.hpp"
#include "../definitions.hpp"
#include "../kmer_block.hpp"
#include "../seq_data.hpp"

//TODO: SL: remove this because it should belong to the suffix array only
bool binarySearch3Prime(cint patternSeqStartPos, unsigned int m, std::pair<cint, cint>& res, const SeqData& S, SuffixArrayWrapper& SA);
bool binarySearch3Prime(const std::string& pattern, KmerBlock& kmerExtended3Prime, const SeqData& S, SuffixArrayWrapper& SA);
