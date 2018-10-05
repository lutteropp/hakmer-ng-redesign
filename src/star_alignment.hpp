/*
 * star_alignment.hpp
 *
 *  Created on: Oct 5, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <vector>
#include <string>

// 1.) Build all O(n^2) pairwise alignments.
// 2.) Choose center sequence S that is closest to all other sequences (this is, sum of edit distances to the other sequences is minimal).
// 3.) Incrementally merge the pairwise alignments into one MSA (Progressive alignment). Gaps in S from the pairwise alignment stay, add gaps to the other sequences as needed.

// Additional difficulty: We need to distinguish between left added flanks, seeds, and right added flanks.
// Still an open question: How to select the alignment penalties (indels, substitituions)?
// Open thing to do: Improve runtime and space-requirements of the pairwise alignments...
// Maybe we could use edlib for very performant pairwise alignment (but it only supports Levenshtein distance)?
