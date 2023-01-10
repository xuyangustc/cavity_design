/*
 * runkmeans.h
 *
 *  Created on: 2017年6月28日
 *      Author: hyliu
 */

#ifndef RUNKMEANS_H_
#define RUNKMEANS_H_
#include "kmeans/kmlocal/KMlocal.h"

#include <vector>
#include <memory>

namespace NSPkmeans {
class RunKMeans {
public:
	struct DataParameters{
		int datadim{0};
		int maxpoints{0};
		int npoints{0};
	};
	struct RunParameters {
		int k{0};
		KMterm term;
		int algorithm{HYBRID};
	};
	struct Cluster {
		std::vector<double> center;
		double distortion;
		std::vector<std::pair<int,double>> members;
	};
	enum ALGORITHM {LLOYDS,SWAP,EZ_HYBRID,HYBRID};
	void init_data(int dim, int maxpoints);
	void addpoint(const std::vector<double> & points);
	void finishdata();
	void setrun(int k, int algorithm);
	void setrun(int k, int algorithm, const KMterm & term);
	double run();
	const std::vector<Cluster> & clusters() const {return clusters_;}
	double getdistortion() const {
		return centers_->getDist();
	}
	~RunKMeans(){;}
private:
	std::shared_ptr<KMdata>  kmdata_;
	std::shared_ptr<KMfilterCenters> centers_;
	std::vector<Cluster> clusters_;
	DataParameters dataparameters_;
	RunParameters runparameters_;
};
}


#endif /* RUNKMEANS_H_ */
