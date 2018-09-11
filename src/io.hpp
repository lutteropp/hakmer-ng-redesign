/*
 * io.hpp
 *
 *  Created on: Sep 10, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <string>
#include "options.hpp"
#include "indexed_concat.hpp"
#include "quartet_lookup_table.hpp"

// read FASTA sequences
// read files containing FASTA contigs
// ---> into an IndexedConcatenatedSequence. (Also taking into account the reverse-complement if activated)

IndexedConcatenatedSequence readConcat(const Options& options);
void printInputStatistics(const std::string& filepath, bool contigs);

// write FASTA supermatrix
// optional: write block coords
// optional: write blocks separately

// read quartet topologies
// write quartet topologies
QuartetLookupTable<size_t> readQuartets(const std::string& filepath);
void writeQuartets(const QuartetLookupTable<size_t>& table, const std::vector<std::string>& taxonLabels, const std::string& filepath);

// that's probably all...
