/*
 * searchtree.h
 *
 *  Created on: 2020��2��13��
 *      Author: xuyang
 */

#ifndef DSTL_SEARCHTREE_H_
#define DSTL_SEARCHTREE_H_

#define BinSize 1.0
#define Half 0.5
//#define TreeMaxVal 30.0
#define TreeMaxVal 50.0

#include <unordered_map>
#include <vector>
#include <functional>
#include <cstdio>
#include <iostream>
#include "math.h"
#include <omp.h>

namespace SearchTree {
class Comp_Double {
public:
	bool operator() (double  v1, double v2) {
		if(v1<v2-Half) return true;
		return false;
	}
};
class Comp_Vector {
public:
	static bool Comp_Double_f(double v1, double v2, double diff) {
		if(v1<v2-diff) return true;
		return false;
	}
	bool operator() (const std::vector<double> &v1, const  std::vector<double> &v2) {
		for(int i=0;i<v1.size();i++) {
			if(Comp_Double_f(v1[i],v2[i],Half)) {
				return true;
			} else if(Comp_Double_f(v2[i],v1[i],Half)) {
				return false;
			}
		}
		return false;
	}
};
class TREE {
public:
	static double GetPoint(double query);
	static std::vector<double> MakeKey(const std::vector<double> &query);
	static std::vector<std::vector<double>> GetBins(
			const std::vector<double> &query, const std::vector<double> &ranges);
	static std::map<int,std::vector<std::vector<double>>> GetBins_FixedNeighbor(
			const std::vector<double> &query, int neighbor);
	static double f10(double val, double rad);
	static double score01(const std::vector<double> &v,
			const std::vector<double> &query, const std::vector<double>&ranges);

	/*static bool Equal_Double(double v1, double v2, double diff) {
		if(fabs(v1-v2)<diff) return true;
		return false;
	}
	struct KEY_COMP {
		bool operator()(const std::vector<double> &k1, const std::vector<double> &k2) const {
			for(int i=0;i<k1.size();i++) {
				if(!Equal_Double(k1[i],k2[i],Half)) {
					return false;
				}
			}
			return true;
		}
	};
	struct HashFunc {
		std::size_t operator() (const std::vector<double> &k) const {
			std::size_t t;
			for(int i=0;i<k.size();i++) {
				if(i==0) t=std::hash<double>()(k[i]);
				else t=t^std::hash<double>()(k[i]);
			}
			return t;
		}
	};

	std::unordered_map<std::vector<double>,std::vector<std::vector<double>>,HashFunc,KEY_COMP> tree;*/
	/*struct KEY {
		std::vector<double> key;
		KEY(const std::vector<double>&k):key(k) {
			;
		}
		bool operator==(const KEY &k) const {
			for(int i=0;i<key.size();i++) {
				if(!Equal_Double(key[i],k.key[i],Half)) {
					return false;
				}
			}
			return true;
		}
	};
	struct HashFunc {
		std::size_t operator() (const KEY &k) const {
			std::size_t t;
			for(int i=0;i<k.key.size();i++) {
				double d=GetPoint(k.key[i]);
				//d += 0.00001;
				//d *= 10.0;
				//int j=(int)d;
				//t=t^std::hash<int>()(j);
				//t=t^std::hash<double>()(d);
				//std::size_t t1=std::hash<int>()(j);
				//std::cout <<j <<' ' <<t1 <<std::endl;
				if(i==0) t=std::hash<double>()(d);
				else t=t^std::hash<double>()(d);
				//std::cout <<t ;
			}
			//std::cout <<std::endl;
			return t;
		}
	};
	std::unordered_map<KEY,std::vector<std::vector<double>>,HashFunc> tree;*/

	std::map<std::vector<double>,std::vector<std::vector<double>>,Comp_Vector> tree;
	void buildtree(const std::vector<double>&dt);
	double score(const std::vector<double> &query, const std::vector<double>&ranges);
	double score(const std::vector<double> &query, const std::vector<double>&ranges,
			std::vector<std::vector<double>>&stages);

	//return value correspond different range
	std::vector<double> score(const std::vector<double> &query, const std::vector<std::vector<double>>&ranges,
			std::vector<std::vector<double>>&stages);
	//range_bin must be the biggest range in ranges
	std::vector<double> score(const std::vector<double> &query,
			const std::vector<std::vector<double>>&ranges, const std::vector<double>&range_bin);

	std::map<std::vector<double>,std::map<std::string,std::vector<std::vector<double>>>,Comp_Vector> tree2;
	std::map<std::string,double> score2(const std::vector<double> &query, const std::vector<double>&ranges,
			const std::vector<std::vector<double>>&stages,
			const std::map<std::string,double>&ws=std::map<std::string,double>());
	std::map<std::string,double> score2(const std::vector<double> &query, const std::vector<double>&ranges,
			const std::map<std::string,double>&ws=std::map<std::string,double>());
	double score21(const std::vector<double> &query, const std::vector<double>&ranges,
			const std::vector<std::vector<double>>&stages, std::string lb);
	double score21(const std::vector<double> &query, const std::vector<double>&ranges, std::string lb);
};
}


#endif /* DSTL_SEARCHTREE_H_ */
