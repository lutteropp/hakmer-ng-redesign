/*
 * log.hpp
 *
 *  Created on: Dec 18, 2018
 *      Author: sarah
 */

#pragma once

#include <fstream>
#include <iostream>
#include <string>

class Logger {
public:
	Logger(const std::string& infopath) {
		if (!infopath.empty()) {
			infoStream.open(infopath);
		}
	}
	void log(const std::string& message) {
#pragma omp critical
		{
			std::cout << message << "\n";
			if (infoStream.good()) {
				infoStream << message << "\n";
			}
		}
	}
	void close() {
		if (infoStream.good()) {
			infoStream.close();
		}
	}
private:
	std::ofstream infoStream;
};
