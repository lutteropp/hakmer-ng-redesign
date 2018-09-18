/*
 * quartet_topology.cpp
 *
 *  Created on: Sep 18, 2018
 *      Author: Sarah Lutteropp
 */

#include "quartet_topology.hpp"

QuartetTopology topologyFromDistances(const std::array<double, 6>& pairwiseDist) {
	double ab_cd = pairwiseDist[0] + pairwiseDist[5];
	double ac_bd = pairwiseDist[1] + pairwiseDist[4];
	double ad_bc = pairwiseDist[2] + pairwiseDist[3];
	if (ab_cd < ac_bd && ab_cd < ad_bc) {
		return QuartetTopology::AB_CD;
	} else if (ac_bd < ab_cd && ac_bd < ad_bc) {
		return QuartetTopology::AC_BD;
	} else if (ad_bc < ab_cd && ad_bc < ac_bd) {
		return QuartetTopology::AD_BC;
	} else {
		return QuartetTopology::STAR;
	}
}
