/*
 * Fswm.hpp
 *
 *  Created on: Jan 8, 2019
 *      Author: sarah
 */

#pragma once

std::vector<std::vector<double> > estimatePairwiseDistances(const std::string& fileName);

inline double harmonicMean(const std::vector<double>& values) {
	double sum = 0;
	for (size_t i = 0; i < values.size(); ++i) {
		sum += (1.0 / values[i]);
	}
	return 1.0 / (sum / values.size());
}

inline double avgSubRate(const std::string& fileName) {
	std::vector<std::vector<double> > pairwiseDist = estimatePairwiseDistances(fileName);
	std::vector<double> vals;
	for (size_t i = 0; i < pairwiseDist.size(); ++i) {
		for (size_t j = i + 1; j < pairwiseDist.size(); ++j) {
			vals.push_back(pairwiseDist[i][j]);
		}
	}
	return harmonicMean(vals);
}
