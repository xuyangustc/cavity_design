/*
 * spectral_cluster.cpp
 *
 *  Created on: Jun 17, 2020
 *      Author: xuyang
 */
#include <algorithm>
#include <iostream>
#include <set>
#include "dstl/spectral_cluster.h"
using namespace NSPdstl;

std::vector<std::vector<double>> SpectralCluster::weighting() {
	std::vector<std::vector<double>> m=matx_;
#pragma omp parallel for schedule(dynamic, 1)
	for(int i=0;i<m.size();i++) {
		for(int j=0;j<m[i].size();j++) {
			double p = m[i][j]/sigma_;
			m[i][j] = exp(-p*p);
		}
	}
	return m;
}

std::vector<SpectralCluster::INFO> SpectralCluster::getinfo() {
	if(!info_.empty()) return info_;
	std::vector<std::vector<double>> ms=weighting();
	info_.resize(ms.size());
#pragma omp parallel for schedule(dynamic, 1)
	for(int i=0;i<ms.size();i++) {
		double s=0.0;
		for(int j=0;j<ms[i].size();j++) {
			s += ms[i][j];
		}
#pragma omp critical
		{
			info_[i] = {{i},{0,s}};
		}
	}
	std::sort(info_.begin(),info_.end(),[](INFO p1,INFO p2)->bool{return p1.second[1]>p2.second[1];});
	info_[0].first.push_back(0);
//#pragma omp parallel for schedule(dynamic, 1)
	for(int i=0;i<info_.size();i++) {
		if(i==0) continue;
		int ii=info_[i].first[0];
		int n;
		double min=100000000.0;
		for(int j=0;j<i;j++) {
			int jj = info_[j].first[0];
			if(matx_[ii][jj]<min) {
				min = matx_[ii][jj];
				n = jj;
			}
		}
		info_[i].first.push_back(n);
		info_[i].second[0] = min;
	}
	//std::cout <<'\t' <<info_.size() <<std::endl;
	return info_;
}

std::vector<std::vector<int>> SpectralCluster::cluster(double cutoff) {
	if(info_.empty()) getinfo();
	std::set<int> cens{0};
	std::vector<std::vector<int>> cls{{info_[0].first[0]}};
	for(int i=1;i<info_.size();i++) {
		int ii=info_[i].first[0]; //index of matx_
		bool in{false};
		for(int j=0;j<cls.size();j++) {
			int jj=cls[j][0]; //index of matx_
			if(matx_[ii][jj]>cutoff) continue;
			in=true;
			break;
		}
		if(!in) {
			cls.push_back({ii});
			cens.insert(ii);
		}
	}
	//std::cout <<cls.size() <<std::endl;
	for(int i=0;i<matx_.size();i++) {
		if(cens.find(i)!=cens.end()) continue;
		int n=-1;
		double min=10000000.0;
		for(int j=0;j<cls.size();j++) {
			int jj=cls[j][0];
			if(matx_[i][jj]<min) {
				n=j;
				min=matx_[i][jj];
			}
		}
		if(min>cutoff) std::cout <<n <<'\t' <<i <<'\t' <<min <<std::endl;
		cls[n].push_back(i);
	}
	return cls;
}
