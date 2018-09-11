/*
 * io.hpp
 *
 *  Created on: Sep 10, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <string>
#include <vector>
#include "options.hpp"
#include "indexed_concat.hpp"
#include "quartet_lookup_table.hpp"

#include "aligned_block.hpp"
#include "extended_block.hpp"
#include "seeded_block.hpp"

// read FASTA sequences
// read files containing FASTA contigs
// ---> into an IndexedConcatenatedSequence. (Also taking into account the reverse-complement if activated)

IndexedConcatenatedSequence readConcat(const Options& options);
void printInputStatistics(const std::string& filepath, bool contigs);

// write FASTA supermatrix
void writeFASTASupermatrix(const std::vector<AlignedBlock>& blocks, const std::vector<std::string>& taxonLabels, const std::string& filepath);
void writeFASTASupermatrix(const std::vector<ExtendedBlock>& blocks, const std::vector<std::string>& taxonLabels, const std::string& filepath, const std::string& T);
void writeFASTASupermatrix(const std::vector<SeededBlock>& blocks, const std::vector<std::string>& taxonLabels, const std::string& filepath, const std::string& T);
// optional: write block coords
// optional: write blocks separately

// read quartet topologies
// write quartet topologies
QuartetLookupTable<size_t> readQuartets(const std::string& filepath);
void writeQuartets(const QuartetLookupTable<size_t>& table, const std::vector<std::string>& taxonLabels, const std::string& filepath);

// that's probably all...
