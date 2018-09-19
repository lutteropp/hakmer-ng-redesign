/*
 * distance_estimator.hpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <string>
#include <stdexcept>

class HammingDistanceEstimator {
public:
	HammingDistanceEstimator() : a(""), b(""), dist(0) {};
	HammingDistanceEstimator(const std::string& a, const std::string& b) {
		if (a.size() != b.size()) {
			throw std::runtime_error("Strings a and b must be of the same size!");
		}
		this->a = a;
		this->b = b;
		dist = 0;
		for (size_t i = 0; i < a.size(); ++i) {
			if (a[i] != b[i]) {
				dist++;
			}
		}
	}
	double distance();
	void addChars(char c1, char c2);
private:
	std::string a;
	std::string b;
	double dist;
};

inline double jukesCantorCorrection(double dist);
