/*
 * quartet_topology.hpp
 *
 *  Created on: Sep 18, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <array>

enum class QuartetTopology {
	AB_CD, AC_BD, AD_BC, STAR
};

QuartetTopology topologyFromDistances(double ab_cd, double ac_bd, double ad_bc);
QuartetTopology topologyFromDistances(const std::array<double, 6>& pairwiseDist);
