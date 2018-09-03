/*
 * dyn_prog.cpp
 *
 *  Created on: Jun 7, 2018
 *      Author: Sarah Lutteropp
 */

#include "dyn_prog.hpp"

void DynProg::dpIteration(int i, int j) {
	int sub = conf->getSubstitutionPenalty(a[i - 1 - dr.getMinRowIdx()], b[j - 1 - dr.getMinColIdx()]);
	int z = std::min(std::min(dr[i - 1][j].l + conf->getDeletionPenalty(), dr[i][j - 1].u + conf->getInsertionPenalty()), sub);
	dr[i][j].u = z - dr[i - 1][j].l;
	dr[i][j].l = z - dr[i][j - 1].u;
}

const std::string& DynProg::getA() {
	return a;
}

const std::string& DynProg::getB() {
	return b;
}

DynProg::DynProg(const std::string& a, const std::string& b, const DistConfig& config) {
	this->a = a;
	this->b = b;

	//std::cout << "create dynprog\n";
	//std::cout << " a = " << a << "\n";
	//std::cout << " b = " << b << "\n";

	conf = &config;
	dr.setMinRowIdx(0);
	dr.setMinColIdx(0);
	dr.setMaxRowIdx(a.size());
	dr.setMaxColIdx(b.size());
	for (int i = dr.getMinRowIdx() + 1; i <= dr.getMaxRowIdx(); ++i) {
		dr[i][dr.getMinColIdx()].u = conf->getDeletionPenalty();
	}
	for (int j = dr.getMinColIdx() + 1; j <= dr.getMaxColIdx(); ++j) {
		dr[dr.getMinRowIdx()][j].l = conf->getInsertionPenalty();
	}
	for (int i = dr.getMinRowIdx() + 1; i <= dr.getMaxRowIdx(); ++i) {
		for (int j = dr.getMinColIdx() + 1; j <= dr.getMaxColIdx(); ++j) {
			dpIteration(i, j);
		}
	}
	distValid_ = false;
	editDist_ = -1;
	normEditDist_ = -1;

	//std::cout << "created dynprog\n";
}

double DynProg::editDistance() {
	if (distValid_) {
		return editDist_;
	}

	double e;
	// corner case: one of the strings is empty
	if (a.size() == 0 && b.size() == 0) {
		e = 0;
	} else if (a.size() == 0 && b.size() > 0) {
		e = b.size() * conf->getInsertionPenalty();
	} else if (a.size() > 0 && b.size() == 0) {
		e = a.size() * conf->getDeletionPenalty();
	} else {
		e = 0;
		if (dr.getNRows() <= dr.getNCols()) { // iterate over i in D[i,n]
			e = (dr.getNCols() - 1) * conf->getInsertionPenalty();
			for (int i = dr.getMinRowIdx() + 1; i <= dr.getMaxRowIdx(); ++i) {
				e += dr[i][dr.getMaxColIdx()].u;
			}
		} else { // iterate over j in D[m,j]
			e = (dr.getNRows() - 1) * conf->getDeletionPenalty();
			for (int j = dr.getMinColIdx() + 1; j <= dr.getMaxColIdx(); ++j) {
				e += dr[dr.getMaxRowIdx()][j].l;
			}
		}
	}
	editDist_ = e;
	normEditDist_ = normalizedEditDistance(e);
	distValid_ = true;
	return e;
}

void DynProg::addCharARight(char c) {
	a += c;
	dr.setMaxRowIdx(dr.getMaxRowIdx() + 1);
	dr[dr.getMaxRowIdx()][dr.getMinColIdx()].u = conf->getDeletionPenalty();
	for (int j = dr.getMinColIdx() + 1; j <= dr.getMaxColIdx(); ++j) {
		dpIteration(dr.getMaxRowIdx(), j);
	}
	distValid_ = false;
}

void DynProg::addCharBRight(char c) {
	b += c;
	dr.setMaxColIdx(dr.getMaxColIdx() + 1);
	dr[dr.getMinRowIdx()][dr.getMaxColIdx()].l = conf->getInsertionPenalty();
	for (int i = dr.getMinRowIdx() + 1; i <= dr.getMaxRowIdx(); ++i) {
		dpIteration(i, dr.getMaxColIdx());
	}
	distValid_ = false;
}

void DynProg::addCharsRight(char cForA, char cForB) {
		addCharARight(cForA);
		addCharBRight(cForB);
}

void DynProg::updateDrColwise() {
	std::deque<int> prevChanged;
	// new first column now
	for (int i = dr.getMinRowIdx() + 1; i <= dr.getMaxRowIdx(); ++i) {
		dr[i][dr.getMinColIdx()].u = conf->getDeletionPenalty();
		prevChanged.push_back(i);
	}
	// set the missing L value
	dr[dr.getMinRowIdx()][dr.getMinColIdx() + 1].l = conf->getInsertionPenalty();

	// go column by column
	for (int j = dr.getMinColIdx() + 1; j <= dr.getMaxColIdx(); ++j) {
		if (prevChanged.empty()) {
			break;
		}

		std::deque<int> currChanged;
		while (!prevChanged.empty()) {
			int i = prevChanged.front();
			prevChanged.pop_front();

			// recompute values
			int oldU = dr[i][j].u;
			int oldL = dr[i][j].l;

			int sub = conf->getSubstitutionPenalty(a[i - 1 - dr.getMinRowIdx()],b[j - 1 - dr.getMinColIdx()]);
			int z = std::min(std::min(dr[i - 1][j].l + conf->getDeletionPenalty(), dr[i][j - 1].u + conf->getInsertionPenalty()),
					sub);
			int newU = z - dr[i - 1][j].l;
			int newL = z - dr[i][j - 1].u;
			dr[i][j].u = newU;
			dr[i][j].l = newL;
			if (oldU != newU) {
				currChanged.push_back(i);
			}

			if (oldL != newL) {
				if (prevChanged.empty() || prevChanged.front() != i + 1) {
					if (i + 1 <= dr.getMaxRowIdx()) {
						prevChanged.push_front(i + 1);
					}
				}
			}
		}
		prevChanged.swap(currChanged);
	}
}

void DynProg::updateDrRowwise() {
	std::deque<int> prevChanged;
	// new first row now
	for (int j = 1 + dr.getMinColIdx(); j <= dr.getMaxColIdx(); ++j) {
		dr[dr.getMinRowIdx()][j].l = conf->getInsertionPenalty();
		prevChanged.push_back(j);
	}

	// set the missing U value
	dr[dr.getMinRowIdx() + 1][dr.getMinColIdx()].u = conf->getDeletionPenalty();

	//go row by row
	for (int i = 1 + dr.getMinRowIdx(); i <= dr.getMaxRowIdx(); ++i) {
		if (prevChanged.empty()) {
			break;
		}

		std::deque<int> currChanged;
		while (!prevChanged.empty()) {
			int j = prevChanged.front();
			prevChanged.pop_front();

			// recompute values
			int oldU = dr[i][j].u;
			int oldL = dr[i][j].l;

			int sub = conf->getSubstitutionPenalty(a[i - 1 - dr.getMinRowIdx()], b[j - 1 - dr.getMinColIdx()]);
			int z = std::min(std::min(dr[i - 1][j].l + conf->getDeletionPenalty(), dr[i][j - 1].u + conf->getInsertionPenalty()),
					sub);
			int newU = z - dr[i - 1][j].l;
			int newL = z - dr[i][j - 1].u;
			dr[i][j].u = newU;
			dr[i][j].l = newL;
			if (oldL != newL) {
				currChanged.push_back(j);
			}

			if (oldU != newU) {
				if (prevChanged.empty() || prevChanged.front() != j + 1) {
					if (j + 1 <= dr.getMaxColIdx()) {
						prevChanged.push_front(j + 1);
					}
				}
			}
		}
		prevChanged.swap(currChanged);
	}
}

void DynProg::addCharALeft(char c) {
	std::string aNew = "";
	aNew += c;
	a = aNew + a;
	dr.setMinRowIdx(dr.getMinRowIdx() - 1);
	updateDrRowwise();
	distValid_ = false;
}

void DynProg::addCharBLeft(char c) {
	std::string bNew = "";
	bNew += c;
	b = bNew + b;
	dr.setMinColIdx(dr.getMinColIdx() - 1);
	updateDrColwise();
	distValid_ = false;
}

void DynProg::addCharsLeft(char cForA, char cForB) {
		addCharALeft(cForA);
		addCharBLeft(cForB);
}

void DynProg::printUMatrix() {
	std::cout << "U matrix begin:\n";
	for (int i = dr.getMinRowIdx() + 1; i <= dr.getMaxRowIdx(); ++i) {
		for (int j = dr.getMinColIdx() + 1; j <= dr.getMaxColIdx(); ++j) {
			std::cout << dr[i][j].u << " ";
		}
		std::cout << "\n";
	}
	std::cout << "U matrix end.\n";
}

void DynProg::printLMatrix() {
	std::cout << "L matrix begin:\n";
	for (int i = dr.getMinRowIdx() + 1; i <= dr.getMaxRowIdx(); ++i) {
		for (int j = dr.getMinColIdx() + 1; j <= dr.getMaxColIdx(); ++j) {
			std::cout << dr[i][j].l << " ";
		}
		std::cout << "\n";
	}
	std::cout << "L matrix end.\n";
}

double DynProg::normalizedEditDistance(double dist) {
	int alpha = std::max(conf->getInsertionPenalty(), conf->getDeletionPenalty());
	return (2.0 * dist) / (alpha * (a.size() + b.size()) + dist);
}

double DynProg::normalizedEditDistance() {
	if (distValid_) {
		return normEditDist_;
	} else {
		editDist_ = editDistance();
		int alpha = std::max(conf->getInsertionPenalty(), conf->getDeletionPenalty());
		normEditDist_ = (2.0 * editDist_) / (alpha * (a.size() + b.size()) + editDist_);
		distValid_ = true;
		return normEditDist_;
	}
}
