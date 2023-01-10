/*
 * searchtree.cpp
 *
 *  Created on: 2020��2��13��
 *      Author: xuyang
 */

#include <map>

#include "dstl/searchtree.h"
using namespace SearchTree;

double TREE::GetPoint(double query) {
	static std::map<double,double,Comp_Double> mp;
#pragma omp critical
	{
		if(mp.empty()) {
			for(double d=Half;d<TreeMaxVal;d+=BinSize) {
				mp.insert({d,d});
			}
			for(double d=-Half;d>-TreeMaxVal;d-=BinSize) {
				mp.insert({d,d});
			}
		}
	}
	if(query>TreeMaxVal || query<-TreeMaxVal) {
		std::cout <<"Value is not in Range!!!" <<std::endl;
		exit(1);
	}
	if(mp.find(query)==mp.end()) {
		std::cout <<"Query Not Found     " <<query <<std::endl;
		if(query<0.0) query += 0.00000001;
		else query -= 0.00000001;
	}
	return mp.at(query);
}

std::vector<double> TREE::MakeKey(const std::vector<double> &query) {
	std::vector<double> key;
	for(double d:query) {
		key.push_back(GetPoint(d));
	}
	return key;
}

std::vector<std::vector<double>> TREE::GetBins(
		const std::vector<double> &query, const std::vector<double> &ranges1) {
	int nd = query.size() / ranges1.size();
	std::vector<double> ranges;
	for(double d:ranges1) {
		for(int i=0;i<nd;i++) ranges.push_back(d);
	}
	//for(double d:query) std::cout <<d <<'\t';std::cout <<std::endl;
	//for(double d:ranges) std::cout <<d <<'\t';std::cout <<std::endl;
	std::vector<std::vector<double>> vals;
	for(int i=0;i<query.size();i++) {
		double low = query[i]-ranges[i];
		double high = query[i]+ranges[i];
		low = GetPoint(low);
		high = GetPoint(high);
		std::vector<double> vd;
		for(double d=low;d<high+0.000001;d+=BinSize) {
			vd.push_back(d);
			//std::cout <<d <<' ';
		}
		//std::cout <<std::endl;
		vals.push_back(vd);
	}
	//std::cout <<vals.size() <<std::endl;

	std::vector<std::vector<double>> bins(1);
	for(int i=0;i<vals.size();i++) {
		bins[0].push_back(vals[i][0]);
	}
	for(int i=0;i<vals.size();i++) {
		int nz = bins.size();
		for(int j=1;j<vals[i].size();j++) {
			for(int k=0;k<nz;k++) {
				std::vector<double> vd = bins[k];
				vd[i] = vals[i][j];
				bins.push_back(vd);
			}
		}
	}
	//std::cout <<bins.size() <<std::endl;
	//exit(1);
	return bins;
}

std::map<int,std::vector<std::vector<double>>> TREE::GetBins_FixedNeighbor(
		const std::vector<double> &query, int neighbor) {
	std::vector<double> q0=MakeKey(query);
	std::vector<std::vector<std::pair<int,double>>> vals;
	for(int i=0;i<query.size();i++) {
		std::vector<std::pair<int,double>> vd;
		for(int j=-neighbor;j<=neighbor;j++) {
			vd.push_back({fabs(j),q0[i]+BinSize*(double)j});
			//std::cout <<vd.back().first <<std::endl;
		}
		vals.push_back(vd);
	}

	std::vector<std::vector<std::pair<int,double>>> bins(1);
	for(int i=0;i<vals.size();i++) {
		bins[0].push_back(vals[i][0]);
	}
	for(int i=0;i<vals.size();i++) {
		int nz = bins.size();
		for(int j=1;j<vals[i].size();j++) {
			for(int k=0;k<nz;k++) {
				std::vector<std::pair<int,double>> vd = bins[k];
				vd[i] = vals[i][j];
				bins.push_back(vd);
			}
		}
	}
	std::map<int,std::vector<std::vector<double>>> rts;
	for(auto &b:bins) {
		std::vector<double> vd;
		int i=0;
		for(auto &p:b) {
			vd.push_back(p.second);
			if(p.first>i) i=p.first;
		}
		if(rts.find(i)==rts.end()) rts.insert({i,{vd}});
		else rts.at(i).push_back(vd);
	}
	return rts;
}
/*
double TREE::f10(double val, double rad) {
	if(val>rad) return 0.0;
	double r1 = 0.1;
	if(val<r1) return 1.0;
	double r2 = (rad + r1) / 2.0;
	double r21 = r2 - r1;
	double a = 0.5 / r21 / r21;
	if(val<r2) {
		double delta = val-r1;
		return 1.0-a*delta*delta;
	}
	double delta = rad - val;
	return a*delta*delta;
}*/

double TREE::f10(double val, double rad) {
	if(val>rad) return 0.0;
	double r1 = 0.0;
	if(val<r1) return 1.0;
	double r2 = (rad + r1) / 2.0;
	double r21 = r2 - r1;
	double a = 0.5 / r21 / r21;
	if(val<r2) {
		double delta = val-r1;
		return 1.0-a*delta*delta;
	}
	double delta = rad - val;
	return a*delta*delta;
}

/*
 * v.size()==query.size()
 * query.size()==nd*ranges.size()
 * if(nd!=1) calculate distance between 2 points
 */
double TREE::score01(const std::vector<double> &v,
		const std::vector<double> &query, const std::vector<double>&ranges) {
	int nd = query.size() / ranges.size();
	double s1=1.0;
	std::vector<double> deltas;
	if(nd==1) {
		for(int i=0;i<v.size();i++) {
			deltas.push_back(fabs(v[i]-query[i]));
			if(deltas.back()>ranges[i]) return 0.0;
		}
		for(int i=0;i<deltas.size();i++) {
			s1 *= f10(deltas[i],ranges[i]);
		}
	} else {
		for(int i=0;i<ranges.size();i++) {
			for(int j=0;j<nd;j++) {
				double k=i*nd+j;
				deltas.push_back(fabs(v[k]-query[k]));
				if(deltas.back()>ranges[i]) return 0.0;
			}
		}
		for(int i=0;i<ranges.size();i++) {
			double dt=0.0;
			for(int j=0;j<nd;j++) {
				double k=i*nd+j;
				dt += deltas[k]*deltas[k];
			}
			dt = sqrt(dt);
			if(dt>ranges[i]) {
				return 0.0;
			} else s1 *= f10(dt,ranges[i]);
		}
	}
	return s1;
}

double TREE::score(const std::vector<double> &query, const std::vector<double>&ranges,
		std::vector<std::vector<double>>&stages) {
	double sc=0.0;
	//for(double d:query) std::cout <<d <<' ';
	//std::cout <<std::endl;
	//for(double d:ranges) std::cout <<d <<' ';
	//std::cout <<std::endl;
	//int temp=0;
	for(const auto &s0:stages) {
		const auto iter=tree.find(s0);
		if(iter==tree.end()) continue;
		//std::cout <<temp++ <<std::endl;
		for(const auto &v:iter->second) {
			double s1 = score01(v,query,ranges);
			//for(double d:v) std::cout <<d <<' ';
			//std::cout <<s1 <<std::endl;
			sc += s1;
		}
	}
	//std::cout <<sc <<std::endl;
	//exit(1);
	return sc;
}

double TREE::score(const std::vector<double> &query, const std::vector<double>&ranges) {
	std::vector<std::vector<double>> stages = GetBins(query,ranges);
	return score(query, ranges, stages);
}

std::vector<double> TREE::score(const std::vector<double> &query, const std::vector<std::vector<double>>&ranges,
		std::vector<std::vector<double>>&stages) {
	std::vector<double> scs(ranges.size(),0.0);
	for(const auto &s0:stages) {
		const auto iter=tree.find(s0);
		if(iter==tree.end()) continue;
		for(const auto &v:iter->second) {
			for(int i=0;i<ranges.size();i++) {
				scs[i] += score01(v,query,ranges[i]);
			}
		}
	}
	return scs;
}

std::vector<double> TREE::score(const std::vector<double> &query,
		const std::vector<std::vector<double>>&ranges, const std::vector<double>&range_bin) {
	std::vector<std::vector<double>> stages = GetBins(query,range_bin);
	return score(query, ranges, stages);
}


std::map<std::string,double> TREE::score2(const std::vector<double> &query,
		const std::vector<double>&ranges, const std::vector<std::vector<double>>&stages,
		const std::map<std::string,double>&ws) {
	std::map<std::string,double> scs;
	for(const auto &s0:stages) {
		const auto iter=tree2.find(s0);
		if(iter==tree2.end()) continue;
		for(const auto &p:iter->second) {
			if(scs.find(p.first)==scs.end()) scs.insert({p.first,0.0});
			auto ir=scs.find(p.first);
			for(const auto &v:p.second) {
				ir->second += score01(v,query,ranges);
			}
		}
	}
	for(auto &w:ws) {
		auto iter=scs.find(w.first);
		if(iter==scs.end()) continue;
		iter->second *= w.second;
	}
	return scs;
}

std::map<std::string,double> TREE::score2(const std::vector<double> &query,
		const std::vector<double>&ranges,const std::map<std::string,double>&ws) {
	std::vector<std::vector<double>> stages = GetBins(query,ranges);
	return score2(query, ranges, stages, ws);
}

double TREE::score21(const std::vector<double> &query, const std::vector<double>&ranges,
		const std::vector<std::vector<double>>&stages, std::string lb) {
	double sc=0.0;
	for(const auto &s0:stages) {
		const auto iter=tree2.find(s0);
		if(iter==tree2.end()) continue;
		auto ir=iter->second.find(lb);
		if(ir==iter->second.end()) continue;
		for(const auto &v:ir->second) {
			sc += score01(v,query,ranges);
		}
	}
	return sc;
}

double TREE::score21(const std::vector<double> &query, const std::vector<double>&ranges, std::string lb) {
	std::vector<std::vector<double>> stages = GetBins(query,ranges);
	return score21(query, ranges, stages, lb);
}

void TREE::buildtree(const std::vector<double>&dt) {
	std::vector<double> key = TREE::MakeKey(dt);
#pragma omp critical
	{
		auto iter=tree.find(key);
		if(iter==tree.end()) tree.insert({key,{dt}});
		else iter->second.push_back(dt);
	}
}











