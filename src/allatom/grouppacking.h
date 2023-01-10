/*
 * grouppacking.h
 *
 *  Created on: Jun 17, 2019
 *      Author: xuyang
 */

#ifndef ALLATOM_GROUPPACKING_H_
#define ALLATOM_GROUPPACKING_H_

#include "allatom/pdbreader_xy.h"
#include "allatom/basic_info.h"
#include "geometry/localframe.h"
#include <memory>
#include <omp.h>
using namespace NSPgeometry;

namespace NSPallatom {


//typedef std::vector<std::string> GroupName; // 0-amino_acid, 1-atom

struct ChemGroup {
	std::string resname;
	//std::vector<std::string> cgname;

	int chainseq;
	int resseq;

	std::map<std::string,XYZ> crds;
	std::vector<std::pair<std::string,XYZ>> crds_v;
	std::set<int> disatom;
	void finddisatom();








	static XYZ RandomRotateWithFixedCenter(std::vector<XYZ>&cs);

	bool aa20cn{false};

	//std::vector<std::string> cgname;

	XYZ cacrd;
	std::vector<int> resseqs; // corresponding with crds_v
	std::vector<int> atomseqs; // use for forcefield

	ChemGroup(int c, int r):chainseq(c), resseq(r) {
		;
	}

	//for Atom3
	bool inmainchain() {
		std::set<std::string> mcas{"N","C","O","CA"};
		for(auto &c:crds_v) {
			if(mcas.find(c.first)==mcas.end()) return false;
		}
		return false;
	}

	const std::vector<std::string> makelabel() const {
		std::vector<std::string> lb{resname};
		for(auto &p:crds_v) lb.push_back(p.first);
		return lb;
	}

	XYZ center;
	void makecenter() {
		center=XYZ(0.0,0.0,0.0);
		for(auto &c:crds_v) {
			center = center + c.second;
		}
		double sz=(double)(crds_v.size());
		center = center / sz;
	}









	std::vector<std::pair<std::string,XYZ>> crd_crd;



	void findcenter();

	void rotate(const Rotation &rt);
	void translate(const XYZ &md);

	std::vector<std::map<std::string,XYZ>> decoys;
	void MakeDecoysWithFixedCenter(int num);
	XYZ findcenter(int n);

	std::vector<std::string> crd_atoms;
	void findcrdatom();
	std::vector<LocalFrame> makelocalframe(const std::map<std::string,XYZ>&cs);
	std::vector<LocalFrame> makelocalframe();
	std::vector<LocalFrame> makelocalframe(int n);
};

bool samechemgroup(const ChemGroup&c1, const ChemGroup&c2);

struct GroupPair {
	static int nbondin2resin1chain(int r1, std::string res1, std::string atom1,
			int r2, std::string res2, std::string atom2);
	static bool farinteraction(const ChemGroup&cg1, const ChemGroup&cg2,int nbond);
	static bool clash(const std::vector<std::string>&n1, const std::vector<XYZ>&cs1,
			const std::vector<std::string>&n2, const std::vector<XYZ>&cs2);
	std::shared_ptr<ChemGroup> g1 { nullptr };
	std::shared_ptr<ChemGroup> g2 { nullptr };
	double dis{-1.0};
	int nbond{-1};
	GroupPair(std::shared_ptr<ChemGroup> c1, std::shared_ptr<ChemGroup> c2):g1(c1), g2(c2) {
		//dis = sqrt((g1->center-g2->center).squarednorm());
	}

	bool samegroup() {
		if(g1->resname!=g2->resname) return false;
		if(g1->crds_v.size()!=g2->crds_v.size()) return false;
		for(int i=0;i<g1->crds_v.size();i++) {
			if(g1->crds_v[i].first!=g2->crds_v[i].first) return false;
		}
		return true;
	}

	double mindis{10000.0};
	void findmindis() {
		mindis=10000.0*10000.0;
		for(int i:g1->disatom) {
			for(int j:g2->disatom) {
				double d2=(g1->crds_v[i].second-g2->crds_v[j].second).squarednorm();
				if(d2<mindis) mindis=d2;
			}
		}
		mindis = sqrt(mindis);
	}
	void changeorder() {
		std::shared_ptr<ChemGroup> temp = g1;
		g1 = g2;
		g2 = temp;
	}
	void rerank();
	void rerank_samegroup(int i, int j);

	const std::vector<std::vector<std::string>> makelabel() const {
		return {g1->makelabel(), g2->makelabel()};
	}

	std::vector<double> atom3par5();
	std::vector<double> atom3par6();

	double minatomdis();



	/*void rerank() {
		if(g1->resname > g2->resname) {
			changeorder();
		} else if(g1->resname == g2->resname) {
			std::vector<std::string> v1, v2;
			for(auto &p:g1->crds) v1.push_back(p.first);
			for(auto &p:g2->crds) v2.push_back(p.first);
			for(int i=0;i<v1.size()&&i<v2.size();i++) {
				if(v1[i]>v2[i]) {
					changeorder();
				} else if(v1[i]<v2[i]) break;
			}
		}
	}*/



	//void rotate_g1(const Rotation &rt);
	//void rotate_g2(const Rotation &rt);
	//void translate_g2(double d);

	bool clash();

	std::vector<XYZ> crd4localframe();
	std::vector<XYZ> crd4localframe(int n);
};



class Atom3 {
public:
	Atom3(const std::vector<std::vector<Residue>>&ress,
			double dis, int itv):chains_(ress), MaxDistance_(dis), MinAtomInterval_(itv) {
		;
	}
	void findatom3();
	void cover2NM();
	void findatom3pair(bool removeclash);
	std::vector<std::shared_ptr<ChemGroup>> cgs() {return cgs_;}
	std::vector<GroupPair> gps() {return gps_;}
	std::vector<GroupPair> gpsinatomdis(double dis);
private:
	std::vector<std::vector<Residue>> chains_;
	double MaxDistance_{12.0};
	int MinAtomInterval_{6};
	std::vector<std::shared_ptr<ChemGroup>> cgs_;
	std::vector<GroupPair> gps_;
	std::vector<GroupPair> gps_inbond;
};




class GroupPacking {
public:
	GroupPacking(const std::vector<std::vector<Residue>>&ress, double dc=10.0):dcutoff_(dc), chains_(ress) {
		;
	}
	void findgroup();
	void findgrouppair(bool removemcmc, bool removeclash, bool samedouble);
	//std::vector<GroupPair>& gps() {return gps_;}
	std::vector<std::shared_ptr<ChemGroup>> cgs() {return cgs_;}
	std::vector<GroupPair>& scsc() {
		return gps_scsc;
	}
	std::vector<GroupPair>& scsc_1() {
		return gps_scsc_1;
	}
	std::vector<GroupPair>& scmc_1() {
		return gps_scmc_1;
	}
	std::vector<GroupPair>& mcsc_1() {
		return gps_mcsc_1;
	}
	std::vector<GroupPair>& mcsc() {
		return gps_mcsc;
	}
	std::vector<GroupPair>& gps() {
		return gps_all;
	}
private:
	std::vector<std::vector<Residue>> chains_;
	double dcutoff_;
	int nbondmin_{6};
	int seqitv_{2};
	std::vector<GroupPair> gps_scsc; // no order
	std::vector<GroupPair> gps_scsc_1;
	std::vector<GroupPair> gps_scmc_1;
	std::vector<GroupPair> gps_mcsc_1;
	std::vector<GroupPair> gps_mcsc; // no order
	std::vector<GroupPair> gps_mcmc; // no order
	std::vector<GroupPair> gps_all;
	std::vector<std::shared_ptr<ChemGroup>> cgs_;

	void addgroup(Residue &res, std::vector<std::string> &ans, std::string resname, int cid, int rid);
	bool gpfilter(GroupPair &gp, bool removeclash=true);
};








//f means 3 par is fixed: BinSize, Ndim, TreeMaxVal
class TREEf {
public:
	double BinSize{1.0};
	double BinSize2{6.0};
	int Ndim{6};
	double TreeMaxVal{30.0};
	std::vector<std::vector<double>> dts;
	static bool Comp_Double_d(double v1, double v2, double diff) {
		if(v1<v2-diff) return true;
		return false;
	}
	class Comp_Double {
	public:
		double half = 0.5;
		bool operator() (double  v1, double v2) {
			//return Comp_Double_d(v1,v2,BinSize/2.0);
			return Comp_Double_d(v1,v2,half);
		}
	};
	class Comp_Vector {
	public:
		bool operator() (const std::vector<double> &v1, const  std::vector<double> &v2) {
			//double half = BinSize / 2.0;
			double half = 0.5;
			double nd=6;
			for(int i=0;i<6;i++) {
				if(Comp_Double_d(v1[i],v2[i],half)) {
					return true;
				} else if(Comp_Double_d(v2[i],v1[i],half)) {
					return false;
				}
			}
			return false;
		}
	};
	class Comp_Vector2 {
	public:
		bool operator() (const std::vector<double> &v1, const  std::vector<double> &v2) {
			//double half = BinSize / 2.0;
			double half = 3.0;
			double nd=6;
			for(int i=0;i<6;i++) {
				if(Comp_Double_d(v1[i],v2[i],half)) {
					return true;
				} else if(Comp_Double_d(v2[i],v1[i],half)) {
					return false;
				}
			}
			return false;
		}
	};
	typedef std::map<std::vector<double>,std::set<long>,Comp_Vector> BRANCH;
	typedef std::map<std::vector<double>,BRANCH,Comp_Vector2> ROOT;
	ROOT tree;
	//std::map<std::vector<double>,std::vector<long>,Comp_Vector> tree;

	double GetPoint(double query);
	double GetPoint2(double query);
	std::vector<double> MakeKey(const std::vector<double> &query);
	std::vector<double> MakeKey2(const std::vector<double> &query);
	void BuildTree();
	void BuildTree_omp();

	std::vector<std::vector<double>> GetBins(
			const std::vector<double> &query, const std::vector<double> &ranges);
	std::set<long> findneighbor(const std::vector<double> &query, const std::vector<double>&ranges);

	std::vector<double> weights;
	double f10_old(double val, double rad);
	double f10(double val, double rad);
	double score1(std::set<long> &si, const std::vector<double> &qy, double rad);
	double score(const std::vector<double> &query, const std::vector<double>&ranges, double rad);
};



}


#endif /* ALLATOM_GROUPPACKING_H_ */
