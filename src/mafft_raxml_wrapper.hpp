#pragma once

#include <string>
#include <array>
#include <vector>

#include "options.hpp"

#include "quartet_topology.hpp"

std::array<std::string, 4> mafftAlign(const std::string& prefix, const std::array<std::string, 4>& sequences);
QuartetTopology inferTopology(const std::string& prefix, const std::vector<std::array<std::string, 4> >& alignedSequences, const Options& options);
QuartetTopology inferTopology(const std::string& prefix, const std::array<std::string, 4>& alignedConcatenatedSequences, const Options& options);

QuartetTopology inferTopologyPlaceholder(const std::string& prefix, const std::vector<std::array<std::string, 4> >& alignedSequences, const Options& options);
