/*
 *      Author: xuyang
 */
#include "allatom/secondary_structure.h"
#include "geometry/calculators.h"
#include "geometry/rotation.h"
#include "backbone/backbonebuilder.h"
#include "dstl/randomengine.h"
#include "dataio/parameters.h"
#include "sd/genchain.h"
#include "allatom/filter.h"
#include "allatom/loop.h"
#include "geometry/quatfit.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
using namespace NSPproteinrep;
using namespace NSPgeometry;
using namespace NSPallatom;
using namespace NSPsd;

std::vector<std::string> ReadLine(std::string parfile) {
	std::vector<std::string> ls;
	std::ifstream ifs(parfile);
	std::string readline;
	while(std::getline(ifs,readline)) {
		if(readline.empty() || readline[0]=='#') continue;
		ls.push_back(readline);
	}
	ifs.close();
	return ls;
}
std::vector<std::vector<Residue>> ReadPDB(std::string pdbfile, bool isnew) {
	PdbReader_xy pr;
	pr.init(pdbfile);
	if(isnew) return pr.chains_new();
	return pr.chains();
}
double RMSD(const std::vector<XYZ>&c0, const std::vector<XYZ>&c1) {
	NSPgeometry::QuatFit qf;
	return sqrt(qf.setup(c0,c1));
}

void Translate(std::vector<Residue> &ch, XYZ &trans) {
	for(Residue &r:ch) {
		for(auto &p:r.rds()) {
			p.second.crd = p.second.crd + trans;
		}
	}
}

void Rotate(std::vector<Residue> &ch, XYZ &fxp, XYZ &axis, double angle) {
	QuaternionCrd qc(axis,angle);
	Rotation rt(qc,fxp);
	for(auto &r:ch) {
		for(auto &c:r.rds()) {
			rt.apply(&(c.second.crd));
		}
	}
}

void Move2Ori(std::vector<Residue>&ch) {
	XYZ cen(0.0,0.0,0.0);
	for(Residue &r:ch) cen = cen + r.rds().at("CA").crd;
	cen = cen / (double)(ch.size());
	cen = -cen;
	Translate(ch,cen);
}

std::vector<Residue> MakeRandomHelix(std::vector<std::string>&seq) {
	std::vector<BackBoneSite> ch = BackBoneBuilder::buildhelixat(seq.size(), XYZ(0.0,0.0,0.0), XYZ(0.0,0.0,1.0), true);
	std::vector<Residue> rs;
	for(int i=0;i<ch.size();i++) rs.push_back(Residue(ch[i],seq[i]));
	Move2Ori(rs);
	return rs;
}

std::vector<Residue> MakeRandomStrand(std::vector<std::string>&seq) {
	std::vector<BackBoneSite> ch = BackBoneBuilder::buildstrandat(seq.size(), XYZ(0.0,0.0,0.0), XYZ(0.0,0.0,1.0), true);
	std::vector<Residue> rs;
	for(int i=0;i<ch.size();i++) rs.push_back(Residue(ch[i],seq[i]));
	Move2Ori(rs);
	return rs;
}

std::vector<Residue> MakeChain(std::vector<std::vector<double>>&dihs, std::vector<std::string>&seq) {
	std::vector<BackBoneSite> chain(dihs.size());
	genbackbonesite(nullptr, false, dihs[0][0], dihs[0][1], &chain[0]);
	for(int i=1;i<dihs.size();i++) {
		genbackbonesite(&chain[i-1], false, dihs[i][0], dihs[i][1], &chain[i]);
	}
	std::vector<Residue> rs;
	for(int i=0;i<chain.size();i++) rs.push_back(Residue(chain[i],seq[i]));
	Move2Ori(rs);
	return rs;
}

struct SS {
	std::string lb{"H"};
	std::string aa{"GLY"};
	bool phipsiisfixed{false};
	int lth_min{0}, lth_max{0};
	double rad{10.0};
	double trans_x{-1.0}, trans_y{-1.0}, trans_z{-1.0};
	double rot{180.0};
	bool doreverse{false};
};
std::vector<std::vector<Residue>> SSAssemble(std::vector<SS> &ss) {
	std::vector<std::vector<Residue>> chs;
	std::vector<double> helix{-57.0,-47.0};
	std::vector<double> parrstrand{-119.0,113.0};
	std::vector<double> antistrand{-139.0,135.0};
	for(int i=0;i<ss.size();i++) {
		std::vector<Residue> ch;
		std::vector<int> lths;
		//for(int l=-ss[i].change;l<=ss[i].change;l++) if(ss[i].lth+l>0) lths.push_back(ss[i].lth+l);
		for(int l=ss[i].lth_min;l<=ss[i].lth_max;l++) lths.push_back(l);
		int nlth = NSPdstl::RandomEngine<>::getinstance().intrng(0,lths.size()-1)();
		int lth = lths[nlth];
		std::vector<std::string> seq(lth,ss[i].aa);
		if(ss[i].phipsiisfixed) {
			std::vector<std::vector<double>> dihs;
			if(ss[i].lb[0]=='H') {
				for(int j=0;j<lth;j++) {
					dihs.push_back(helix);
				}
			} else if(ss[i].lb[0]=='E') {
				for(int j=0;j<lth;j++) {
					int pa = NSPdstl::RandomEngine<>::getinstance().intrng(0,1)();
					if(pa==0) dihs.push_back(parrstrand);
					else dihs.push_back(antistrand);
				}
			} else {
				std::cout <<"Can Not Recognize Label: " <<ss[i].lb <<std::endl;
				exit(1);
			}
			ch = MakeChain(dihs,seq);
		} else {
			if(ss[i].lb[0]=='H') {
				ch = MakeRandomHelix(seq);
			} else if(ss[i].lb[0]=='E') {
				ch = MakeRandomStrand(seq);
			} else {
				std::cout <<"Can Not Recognize Label: " <<ss[i].lb <<std::endl;
				exit(1);
			}
		}
		if(ch.empty()) {
			std::cout <<"Wrong!" <<std::endl;
			exit(1);
		}
		if(ss[i].doreverse) {
			XYZ axis_temp(1.0,0.0,0.0);
			XYZ fxp_temp(0.0,0.0,0.0);
			double ang_temp=180.0;
			Rotate(ch,fxp_temp,axis_temp,ang_temp);
		}
		XYZ cen(0.0,0.0,0.0);
		double maxangle=ss[i].rot;
		auto &rng=NSPdstl::RandomEngine<>::getinstance();
		XYZ axis=XYZ(rng.realrng(0.0,1.0),1.0);
		double angle=rng.realrng(0.0,maxangle*2.0)()-maxangle;
		Rotate(ch,cen,axis,angle);
		double maxtrans=ss[i].rad;
		XYZ trans=XYZ(rng.realrng(0.0,1.0),maxtrans,ss[i].trans_x,ss[i].trans_y,ss[i].trans_z);
		Translate(ch,trans);
		chs.push_back(ch);
	}
	return chs;
}
void SSAssemble(std::string parfile) {
	std::vector<std::string> ps = ReadLine(parfile);
	int seed=std::stoi(ps[0]);
	NSPdstl::RandomEngine<>::getinstance().reseed(seed);
	int nresult=std::stoi(ps[1]);
	std::vector<SS> vss;
	for(int i=2;i<ps.size();i+=2) {
		std::stringstream str(ps[i]);
		int nout;
		SS ss;
		bool dorestrict;
		str >>nout >>ss.lb >>ss.aa >>ss.phipsiisfixed >>ss.lth_min >>ss.lth_max >>ss.rad >>dorestrict;
		if(dorestrict) {
			str.clear();
			str <<ps[i+1];
			str >>ss.trans_x  >>ss.trans_y >>ss.trans_z >>ss.rot >>ss.doreverse;
		}
		for(int j=0;j<nout;j++) vss.push_back(ss);
	}
	for(int i=0;i<nresult;i++) {
		std::vector<std::vector<Residue>> chs=SSAssemble(vss);
		NSPallatom::residue2pdb("assemble_"+std::to_string(i)+".pdb",chs);
	}
}










DECOY pdb2decoy(std::string pdbfile) {
	PdbReader_xy pr;
	pr.init(pdbfile);
	return pr.chains_new();
}

void printpdb0(std::string outfile, const std::vector<std::vector<Residue>>& outpep) {
	std::ofstream ofs(outfile);
	int aid=0;
	for(int i=0;i<outpep.size();i++) {
		char cid='A'+i;
		for(int j=0;j<outpep[i].size();j++) {
			std::vector<PdbRecord> prs=outpep[i][j].getpdbrecords();
			for(PdbRecord &p:prs) {
				p.residueid = j;
				p.atomid = aid++;
				p.chainid = cid;
				ofs <<p.toString() <<std::endl;
			}
		}
	}
	ofs.close();
}

std::vector<Residue> make_chain(std::vector<std::vector<double>>&dihs, std::vector<std::string>&seq) {
	std::vector<BackBoneSite> chain(dihs.size());
	genbackbonesite(nullptr, false, dihs[0][0], dihs[0][1], &chain[0]);
	for(int i=1;i<dihs.size();i++) {
		genbackbonesite(&chain[i-1], false, dihs[i][0], dihs[i][1], &chain[i]);
	}
	std::vector<Residue> rs;
	for(int i=0;i<chain.size();i++) rs.push_back(Residue(chain[i],seq[i]));
	return rs;
}

void makess(int label, int lth, std::string outfile, std::string aa) {
	std::pair<double,double> helix{-57.0,-47.0};
	std::pair<double,double> parrstrand{-119.0,113.0};
	std::pair<double,double> antistrand{-139.0,135.0};
	std::map<int,std::pair<double,double>> dihs;
	dihs.insert({0,helix});
	dihs.insert({1,parrstrand});
	dihs.insert({2,antistrand});
	//XYZ start(st,0.0,0.0), dir(0.0,0.0,1.0);
	double phi=dihs.at(label).first;
	double psi=dihs.at(label).second;
	std::vector<std::vector<double>> dh(lth,{phi,psi});
	std::vector<std::string> seq(lth,aa);
	std::vector<Residue> rs=make_chain(dh,seq);
	std::vector<std::vector<Residue>> rss{rs};
	printpdb0(outfile,rss);
}










void make_random_chain(int lth, std::string outfile, std::string aa, int seed) {
	NSPdstl::RandomEngine<>::getinstance().reseed(seed);
	auto &rng=NSPdstl::RandomEngine<>::getinstance().realrng(0.0,360.0);
	std::vector<std::vector<double>> dh;
	std::vector<std::string> seq;
	for(int i=0;i<lth;i++) {
		double phi=rng()-180.0;
		double psi=rng()-180.0;
		dh.push_back({phi,psi});
		seq.push_back(aa);
	}
	std::vector<Residue> rs=make_chain(dh,seq);
	std::vector<std::vector<Residue>> rss{rs};
	printpdb0(outfile,rss);
}
std::vector<Residue> make_random_chain(int lth,
		std::vector<std::pair<int,int>> &hreg,
		std::vector<std::pair<int,int>> &sreg, std::set<int> &cissites) {
	BackBoneSite stbbs;
	genbackbonesite(nullptr, false, 0.0, 0.0, &stbbs);
	std::vector<BackBoneSite> bs=BackBoneBuilder::buildforwardbackbone(
			lth, stbbs, hreg, sreg, cissites);
	std::vector<Residue> ch;
	for(int i=0;i<bs.size();i++) {
		ch.push_back(Residue(bs[i],"GLY"));
	}
	return ch;
}
void make_random_chain(std::string infile, std::string outpdb, int seed) {
	std::vector<std::pair<int,int>> hreg, sreg;
	int lth=0;
	std::ifstream ifs(infile);
	std::string readline;
	while(std::getline(ifs,readline)) {
		if(readline.empty()) continue;
		if(readline[0]=='#') continue;
		std::stringstream ss(readline);
		std::string lb;
		int l;
		ss >>lb >>l;
		if(lb[0]=='H') hreg.push_back({lth,l});
		else if(lb[0]=='S') sreg.push_back({lth,l});
		lth +=l;
	}
	ifs.close();
	std::set<int> cissites;
	NSPdstl::RandomEngine<>::getinstance().reseed(seed);
	std::vector<std::vector<Residue>> chs{make_random_chain(lth, hreg, sreg, cissites)};
	residue2pdb(outpdb,chs);
}
bool FindClash(const std::vector<Residue> &ch, const std::vector<XYZ> &crd_clash) {
	double dis2=6.0*6.0;
	for(const Residue &r:ch) {
		for(const auto &p:r.rds()) {
			for(const XYZ &c:crd_clash) {
				double d2=(p.second.crd-c).squarednorm();
				if(d2<dis2) return true;
			}
		}
	}
	dis2=2.0*2.0;
	for(int i=0;i<ch.size();i++) {
		for(int j=i+5;j<ch.size();j++) {
			for(auto &pi:ch[i].rds()) {
				for(auto &pj:ch[j].rds()) {
					double d2=(pi.second.crd-pj.second.crd).squarednorm();
					if(d2<dis2) return true;
				}
			}
		}
	}
	return false;
}
void make_random_chain_assemble_Helix(std::string infofile) {
	std::vector<std::string> ls =ReadLine(infofile);
	int seed=std::stoi(ls[0]);
	NSPdstl::RandomEngine<>::getinstance().reseed(seed);
	int nout=std::stoi(ls[1]);
	int nh_min=std::stoi(ls[2]), nh_max=std::stoi(ls[3]);
	int minh=std::stoi(ls[4]), maxh=std::stoi(ls[5]);
	int minl=std::stoi(ls[6]), maxl=std::stoi(ls[7]);
	std::string outprefix=ls[8];

	XYZ crd_center;
	std::vector<XYZ> crd_clash;
	for(int i=9;i<ls.size();i++) {
		std::stringstream ss(ls[i]);
		XYZ c;
		ss >>c.x_ >>c.y_ >>c.z_;
		if(i==9) crd_center=c;
		else crd_clash.push_back(c);
	}
	auto &rng=NSPdstl::RandomEngine<>::getinstance();
	for(int i=0;i<nout;i++) {
		int nh=rng.intrng(nh_min,nh_max)();
		int lth=0;
		std::vector<std::pair<int,int>> hreg, sreg;
		std::set<int> cissites;
		hreg.push_back({lth,rng.intrng(minh,maxh)()});
		lth+=hreg.back().second;
		for(int j=0;j<nh-1;j++) {
			lth+=rng.intrng(minl,maxl)();
			hreg.push_back({lth,rng.intrng(minh,maxh)()});
			lth+=hreg.back().second;
		}
		std::cout <<i;
		for(auto &p:hreg) std::cout <<' ' <<p.second;
		std::cout <<' ' <<lth <<std::endl;
		std::vector<Residue> ch;
		do {
			ch=make_random_chain(lth, hreg, sreg, cissites);
			XYZ cen(0.0,0.0,0.0);
			for(Residue &r:ch) cen =cen+r.rds().at("CA").crd;
			cen =cen/(double)(ch.size());
			XYZ trans=crd_center-cen;
			for(Residue &r:ch) {
				for(auto &p:r.rds()) p.second.crd =p.second.crd+trans;
			}
		} while(FindClash(ch,crd_clash));
		std::vector<std::vector<Residue>> chs{ch};
		residue2pdb(outprefix+std::to_string(i),chs);
	}
}

int main(int argc, char ** argv) {
	make_random_chain_assemble_Helix(std::string(argv[1]));
}

