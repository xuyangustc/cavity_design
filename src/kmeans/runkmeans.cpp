/*
 * runkmeans.cpp
 *
 *  Created on: 2017年6月28日
 *      Author: hyliu
 */

#include "kmeans/runkmeans.h"

using namespace NSPkmeans;

void RunKMeans::init_data(int dim, int maxpoints) {
	dataparameters_.datadim=dim;
	dataparameters_.maxpoints=maxpoints;
	dataparameters_.npoints=0;
	kmdata_=std::shared_ptr<KMdata> (new KMdata(dim,maxpoints));
}
void RunKMeans::addpoint(const std::vector<double> & point){
	assert(point.size()==dataparameters_.datadim);
	assert(dataparameters_.npoints< dataparameters_.maxpoints);
	for(int d=0;d<dataparameters_.datadim; ++d){
		(*kmdata_)[dataparameters_.npoints][d]=point[d];
	}
	++dataparameters_.npoints;
}

void RunKMeans::setrun(int k, int algorithm) {
 setrun(k,algorithm,KMterm(100, 0, 0, 0,		// run for 100 stages
		0.10,			// min consec RDL
		0.10,			// min accum RDL
		3,			// max run stages
		0.50,			// init. prob. of acceptance
		10,			// temp. run length
		0.95));
}
void RunKMeans::finishdata(){
	kmdata_->setNPts(dataparameters_.npoints);
    kmdata_->buildKcTree();
}
void RunKMeans::setrun(int k, int algorithm,const KMterm &term){
    runparameters_.k=k;
    runparameters_.algorithm=algorithm;
	runparameters_.term=term;
	centers_= std::shared_ptr<KMfilterCenters>(
			new KMfilterCenters(runparameters_.k, *kmdata_));
}

double RunKMeans::run(){
	std::shared_ptr<KMlocal> kmlocal;
	switch(runparameters_.algorithm){
	case LLOYD:
		kmlocal=std::shared_ptr<KMlocal> (new KMlocalLloyds(*centers_,runparameters_.term));
		break;
	case SWAP:
		kmlocal=std::shared_ptr<KMlocal> (new KMlocalSwap(*centers_,runparameters_.term));
		break;
	case EZ_HYBRID:
		kmlocal=std::shared_ptr<KMlocal> (new KMlocalEZ_Hybrid(*centers_,runparameters_.term));
		break;
	case HYBRID:
		kmlocal=std::shared_ptr<KMlocal> (new KMlocalHybrid(*centers_,runparameters_.term));
		break;
	default:
		;
	}
	*centers_=kmlocal->execute();
	KMctrIdxArray closeCtr = new KMctrIdx[dataparameters_.npoints];
	double* sqDist = new double[dataparameters_.npoints];
	centers_->getAssignments(closeCtr, sqDist);
	double *dists=centers_->getDists();
	clusters_.clear();
	for(int i=0;i<runparameters_.k;++i) {
		clusters_.push_back(Cluster());
		for(int d=0; d<dataparameters_.datadim;++d)
			clusters_.back().center.push_back((*centers_)[i][d]);
		clusters_.back().distortion=dists[i];
	}
	for(int i=0; i<dataparameters_.npoints; ++i) {
		clusters_[closeCtr[i]].members.push_back(std::make_pair(i,sqDist[i]));
	}
	delete[] closeCtr;
	delete[] sqDist;
	return  centers_->getDist(false)/(double)dataparameters_.npoints;
}
