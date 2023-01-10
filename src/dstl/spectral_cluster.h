/*
 * spectral_cluster.h
 *
 *  Created on: Jun 17, 2020
 *      Author: xuyang
 */

#ifndef DSTL_SPECTRAL_CLUSTER_H_
#define DSTL_SPECTRAL_CLUSTER_H_
#include <vector>
#include <cmath>

namespace NSPdstl{
class SpectralCluster {
public:
	typedef std::pair<std::vector<int>,std::vector<double>> INFO;
	//int0-serial, int1-nearest, d0-total, d1-nearest_distance
	SpectralCluster() {
		;
	}
	SpectralCluster(const std::vector<std::vector<double>>&m, double s): matx_(m), sigma_(s) {
		;
	}
	std::vector<INFO> getinfo();
	std::vector<std::vector<int>> cluster(double cutoff);
private:
	std::vector<std::vector<double>> matx_;
	double sigma_{1.0};
	std::vector<INFO> info_;
	std::vector<std::vector<double>> weighting();
};
}


#endif /* DSTL_SPECTRAL_CLUSTER_H_ */
