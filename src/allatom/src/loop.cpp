/*
 * loop.cpp
 *
 *  Created on: Jun 13, 2020
 *      Author: xuyang
 */
#include "allatom/loop.h"

using namespace NSPproteinrep;
using namespace NSPgeometry;
using namespace NSPsd;
using namespace NSPallatom;

//replace part of chains with lps
//chains.size() == lps.size(), which present number of chain
//chsten[0] : chst, resst, chen, resen
//loop_added.size() == chains.size()
/*std::vector<std::vector<Residue>> LOOP::NewChains(const std::vector<std::vector<Residue>>&chains,
			const std::vector<std::vector<Residue>>&lps, const std::vector<std::vector<int>>&chsten,
			std::vector<std::vector<int>>&loop_added) {
	assert(lps.size()==chsten.size());
	std::set<std::vector<int>> removed;
	for(int i=0;i<chsten.size();i++) {
		int c0=chsten[i][0];
		int st=chsten[i][1];
		int c1=chsten[i][2];
		int en=chsten[i][3];
		if(c0==c1) {
			for(int j=st+1;j<en;j++) {
				removed.insert({c0,j});
			}
		} else {
			for(int j=st+1;j<chains[c0].size();j++) {
				removed.insert({c0,j});
			}
			for(int j=0;j<en;j++) {
				removed.insert({c1,j});
			}
		}
	}
	std::map<std::vector<int>,Residue> ress;
	for(int i=0;i<chains.size();i++) {
		for(int j=0;j<chains[i].size();j++) {
			if(removed.find({i,j})==removed.end()) ress.insert({{i,j},chains[i][j]});
		}
	}
	//std::cout <<chains.size() <<' ' <<chains[0].size() <<' ' <<chains[1].size() <<std::endl;
	//std::cout <<lps.size() <<' ' <<lps[0].size() <<' ' <<lps[1].size() <<std::endl;
	//std::cout <<chsten.size() <<std::endl;
	//std::cout <<chsten[0].size() <<' ' <<chsten[0][0] <<' ' <<chsten[0][1] <<' ' <<chsten[0][2] <<' ' <<chsten[0][3] <<std::endl;
	//std::cout <<chsten[1].size() <<' ' <<chsten[1][0] <<' ' <<chsten[1][1] <<' ' <<chsten[1][2] <<' ' <<chsten[1][3] <<std::endl;
	//exit(1);
	std::vector<std::vector<Residue>> chs;
	while(!ress.empty()) {
		std::vector<Residue> ch;
		int cst=0;
		int rst=0;
		bool find{false};
		for(int i=0;i<chains.size();i++) {
			for(int j=0;j<chains[i].size();j++) {
				if(ress.find({i,j})==ress.end()) continue;
				//ch.push_back(chains[i][j]);
				//ress.erase({i,j});
				find=true;
				cst = i;
				rst = j;
				//std::cout <<"start: " <<i <<' ' <<j <<std::endl;
				break;
			}
			if(find) break;
		}
		//int ntime=0;
		while(true) {
			for(int i=rst;i<chains[cst].size();i++) {
				if(ress.find({cst,i})==ress.end()) break;
				ch.push_back(chains[cst][i]);
				ress.erase({cst,i});
				rst=i;
			}
			//std::cout <<cst <<" " <<rst <<" " <<ress.size() <<" " <<ch.size() <<std::endl;
			//std::cout <<chsten.size() <<" " <<cst <<" " <<rst <<std::endl;
			bool findloop{false};
			for(int i=0;i<chsten.size();i++) {
				//std::cout <<chsten.size() <<" " <<chsten[i][0] <<" " <<cst <<" " <<chsten[i][1] <<" " <<rst <<std::endl;
				if(chsten[i][0]==cst && chsten[i][1]==rst) {
					loop_added.push_back({chs.size(),ch.size(),ch.size()+lps[i].size()-1});
					for(const Residue &r:lps[i]) {
						ch.push_back(r);
					}
					cst=chsten[i][2];
					rst=chsten[i][3];
					//std::cout <<chsten[i][0] <<' ' <<chsten[i][1] <<' ' <<chsten[i][2] <<' ' <<chsten[i][3] <<std::endl;
					findloop=true;
					break;
				}
			}
			if(!findloop) break;
		}
		//break;
		//std::cout <<ch.size() <<std::endl;
		chs.push_back(ch);
	}
	//std::cout <<chs.size() <<std::endl;
	//NSPallatom::residue2pdb("out.pdb",chs);
	//exit(1);
	return chs;
}*/

std::vector<std::vector<Residue>> LOOP::NewChains(const std::vector<std::vector<Residue>>&chains,
			const std::vector<std::vector<Residue>>&lps, const std::vector<std::vector<int>>&chsten,
			std::vector<std::vector<int>>&loop_added) {
	assert(lps.size()==chsten.size());

	std::map<std::vector<int>,std::pair<std::vector<int>,std::vector<Residue>>> starts;
	std::map<int,std::vector<int>> stencomp;
	for(int i=0;i<chains.size();i++) stencomp.insert({i,{10000,10001}});
	for(int i=0;i<chsten.size();i++) {
		int cst=chsten[i][0];
		int rst=chsten[i][1];
		int cen=chsten[i][2];
		int ren=chsten[i][3];
		starts.insert({{cst,rst},{{cen,ren},lps[i]}});
		if(rst<stencomp.at(cst)[0]) {
			stencomp.at(cst)[0] = rst;
		}
		if(ren<stencomp.at(cen)[1]) {
			stencomp.at(cen)[1] = ren;
		}
	}

	//for one chain, if the min of start is smaller than min of end, then the first residue of this chain is a new start
	std::vector<int> new_chain_start;
	for(auto &p:stencomp) {
		if(p.second[0]<p.second[1]) new_chain_start.push_back(p.first);
	}
	if(new_chain_start.empty()) {
		std::cout <<"Can Not Find Starting Residue in the New Chain!!!" <<std::endl;
		exit(1);
	}

	std::vector<std::vector<Residue>> chs;
	for(int i=0;i<new_chain_start.size();i++) {
		std::vector<Residue> ch;
		int cst=new_chain_start[i];
		int rst=0;
		int cen=-1;
		int ren=-1;
		for(int j=0;j<chains[cst].size();j++) {
			ch.push_back(chains[cst][j]);
			if(starts.find({cst,j})!=starts.end()) {
				rst = j;
				cen = starts.at({cst,j}).first[0];
				ren = starts.at({cst,j}).first[1];
				break;
			}
		}
		if(cen==-1) {
			chs.push_back(ch);
			continue;
		}
		while(true) {
			std::vector<Residue> vr=starts.at({cst,rst}).second;
			loop_added.push_back({chs.size(),ch.size(),ch.size()+vr.size()-1});
			for(Residue &r:vr) ch.push_back(r);
			bool findstart{false};
			for(int j=ren;j<chains[cen].size();j++) {
				ch.push_back(chains[cen][j]);
				if(starts.find({cen,j})!=starts.end()) {
					cst = cen;
					rst = j;
					cen = starts.at({cst,j}).first[0];
					ren = starts.at({cst,j}).first[1];
					findstart=true;
					break;
				}
			}
			if(!findstart) break;
		}
		chs.push_back(ch);
	}
	return chs;
}

std::vector<std::vector<Residue>> LOOP::MakeNLoop(const std::vector<std::vector<Residue>>&chains,
		int nloop, int c0, int st, int c1, int en, int lth, int maxtry) {
	std::vector<std::pair<int,int>> hrs;
	std::vector<std::pair<int,int>> srs;
	return MakeNLinker(chains,hrs,srs,nloop,c0,st,c1,en,lth,maxtry);
}

std::vector<std::vector<Residue>> LOOP::MakeNLinker(const std::vector<std::vector<Residue>>&chains,
		const std::vector<std::pair<int,int>> & helixregions, //first is start position, second is length
		const std::vector<std::pair<int,int>> & strandregions,//first is start position, second is length
		int nloop, int c0, int st, int c1, int en, int lth, int maxtry) {
	BackBoneSite bs0=getbbs(chains,c0,st);
	BackBoneSite bs1=getbbs(chains,c1,en);
	std::set<int> cis;
	int ntime=0;
	std::vector<std::vector<BackBoneSite>> bbss;
	do {
		ntime++;
		auto bb=BackBoneBuilder::buildlinkers(lth, bs0, bs1, helixregions, strandregions, cis);
		for(auto &b:bb) bbss.push_back(*b);
		if(bbss.size()>=nloop) break;
	} while(ntime<maxtry);
	std::vector<std::vector<Residue>> lps;
	for(int i=0;i<nloop && i<bbss.size();i++) {
		std::vector<Residue> lp1;
		for(auto &bs:bbss[i]) lp1.push_back(Residue(bs,"GLY"));
		lps.push_back(lp1);
	}
	return lps;
}

BackBoneSite LOOP::getbbs(const std::vector<std::vector<Residue>>&chains, int c, int r) {
	BackBoneSite bs=chains[c][r].getbackbonesite();
	bs.data_[BackBoneSite::OMIGA] =180.0;
	if(r==0) {
		BackBoneSite bs0=chains[c][1].getbackbonesite();
		bs.data_[BackBoneSite::PHI] = bs0.phi(bs);
	} else {
		BackBoneSite bs0=chains[c][r-1].getbackbonesite();
		bs.data_[BackBoneSite::PHI] = bs.phi(bs0);
	}
	if(r==chains[c].size()-1) {
		bs.resetpsi();
	} else {
		BackBoneSite bs1=chains[c][r+1].getbackbonesite();
		bs.data_[BackBoneSite::PSI] = bs.psi(bs1);
	}
	return bs;
}

std::map<int,std::pair<double,SearchTree::TREE>> LoopLengthAndLocation::Build1Tree(std::string fn) {
	std::map<int,std::pair<double,SearchTree::TREE>> tes;
	std::map<int,int> tot;
	std::ifstream ifs(fn);
	std::string readline;
	while(std::getline(ifs,readline)) {
		std::stringstream ss(readline);
		int j;
		std::vector<double> v(6);
		ss >>j >>v[0] >>v[0] >>v[1] >>v[2] >>v[3] >>v[4] >>v[5];
		if(j<1 || j>10) continue;
		if(tot.find(j)==tot.end()) tot.insert({j,1});
		else tot.at(j)++;
		std::vector<double> key = SearchTree::TREE::MakeKey(v);
		if(tes.find(j)==tes.end()) tes.insert({j,{0.0,SearchTree::TREE()}});
		auto &te = tes.at(j).second.tree;
		auto iter=te.find(key);
		if(iter==te.end()) te.insert({key,{v}});
		else iter->second.push_back(v);
	}
	ifs.close();
	for(auto &p:tes) {
		p.second.first = 10000.0 / (double)(tot.at(p.first));
	}
	return tes;
}
std::map<std::string,std::map<int,std::pair<double,SearchTree::TREE>>> LoopLengthAndLocation::BuildTree() {
	static std::map<std::string,std::map<int,std::pair<double,SearchTree::TREE>>> trees;
	if(trees.empty()) {
		trees.insert({"EE",Build1Tree(NSPdataio::datafilename("looplengthandlocation/ee"))});
		trees.insert({"EH",Build1Tree(NSPdataio::datafilename("looplengthandlocation/eh"))});
		trees.insert({"HE",Build1Tree(NSPdataio::datafilename("looplengthandlocation/he"))});
		trees.insert({"HH",Build1Tree(NSPdataio::datafilename("looplengthandlocation/hh"))});
	}
	return trees;
}
std::map<std::string,std::map<int,double>> LoopLengthAndLocation::Criterion(double maxval) {
	static std::map<std::string,std::map<int,double>> criterion;
	if(criterion.size()!=0) return criterion;
	std::ifstream ifs(NSPdataio::datafilename("looplengthandlocation/cutoff20"));
	std::string readline;
	while(std::getline(ifs,readline)) {
		std::stringstream ss(readline);
		std::string key;
		int n;
		ss >>key >>n;
		ss.clear();
		std::map<int,double> mid;
		for(int i=0;i<n;i++) {
			std::getline(ifs,readline);
			std::stringstream ss1(readline);
			ss1 <<readline;
			int j;
			double c;
			ss1 >>j >>c;
			ss1.clear();
			if(c>maxval) c=maxval;
			mid.insert({j,c});
		}
		criterion.insert({key,mid});
	}
	ifs.close();
	return criterion;
}
std::map<std::string,std::map<int,std::vector<double>>> LoopLengthAndLocation::DisCut() {
	static std::map<std::string,std::map<int,std::vector<double>>> discut;
	if(discut.size()!=0) return discut;
	std::ifstream ifs(NSPdataio::datafilename("looplengthandlocation/distance_cutoff"));
	std::string readline;
	while(std::getline(ifs,readline)) {
		std::stringstream ss(readline);
		//std::cout <<readline <<std::endl;
		std::string key;
		int n;
		ss >>key >>n;
		ss.clear();
		std::map<int,std::vector<double>> mid;
		for(int i=0;i<n;i++) {
			std::getline(ifs,readline);
			//std::cout <<readline <<std::endl;
			std::stringstream ss1(readline);
			int j;
			std::vector<double> c(2);
			ss1 >>j >>c[0] >>c[1];
			ss1.clear();
			//std::cout <<j <<' ' <<c[0] <<' ' <<c[1] <<std::endl;
			mid.insert({j,c});
		}
		discut.insert({key,mid});
	}
	ifs.close();
	return discut;
}
std::vector<double> LoopLengthAndLocation::Ranges() {
	static std::vector<double> ranges;
	if(ranges.size()!=0) return ranges;
	for(int i=0;i<6;i++) ranges.push_back(4.0);
	return ranges;
}





LocalFrame LoopLengthAndLocation::nHelix(int c, int r, XYZ &c0) {
	c0 = chains[c][r].rds().at("CA").crd;
	XYZ c1 = chains[c][r-2].rds().at("CA").crd;
	XYZ c2 = chains[c][r-4].rds().at("CA").crd;
	return make_localframe(c0,c1,c2);
}
LocalFrame LoopLengthAndLocation::cHelix(int c, int r, XYZ &c0) {
	c0 = chains[c][r].rds().at("CA").crd;
	XYZ c1 = chains[c][r+2].rds().at("CA").crd;
	XYZ c2 = chains[c][r+4].rds().at("CA").crd;
	return make_localframe(c0,c1,c2);
}
LocalFrame LoopLengthAndLocation::nStrand(int c, int r, XYZ &c0) {
	c0 = chains[c][r].rds().at("CA").crd;
	XYZ c1 = chains[c][r].cb();
	XYZ c2 = chains[c][r-1].rds().at("CA").crd;
	return make_localframe(c0,c1,c2);
}
LocalFrame LoopLengthAndLocation::cStrand(int c, int r, XYZ &c0) {
	c0 = chains[c][r].rds().at("CA").crd;
	XYZ c1 = chains[c][r].cb();
	XYZ c2 = chains[c][r+1].rds().at("CA").crd;
	return make_localframe(c0,c1,c2);
}
std::vector<double> LoopLengthAndLocation::AddElement(XYZ&c0, LocalFrame&lf0, XYZ&c1, LocalFrame&lf1) {
	XYZ c00 = lf1.global2localcrd(c0);
	XYZ c11 = lf0.global2localcrd(c1);
	std::vector<double> v{c00.x_, c00.y_, c00.z_, c11.x_, c11.y_, c11.z_};
	return  v;
}

std::map<int,double> LoopLengthAndLocation::score(int c0, int r0, int c1, int r1, int minlth, int maxlth) {
	std::map<int,double> scs;
	std::string key{secstr[c0][r0],secstr[c1][r1]};
	if(secstr[c0][r0]=='H') {
		for(int i=0;i<5;i++) {
			if(secstr[c0][r0-i]!='H') return scs;
		}
	} else if(secstr[c0][r0]=='E') {
		for(int i=0;i<3;i++) {
			if(secstr[c0][r0-i]!='E') return scs;
		}
	} else return scs;
	if(secstr[c1][r1]=='H') {
		for(int i=0;i<5;i++) {
			if(secstr[c1][r1+i]!='H') return scs;
		}
	} else if(secstr[c1][r1]=='E') {
		for(int i=0;i<3;i++) {
			if(secstr[c1][r1+i]!='E') return scs;
		}
	} else return scs;

	//std::cout <<1 <<std::endl;
	XYZ crd0;
	LocalFrame lf0;
	if(secstr[c0][r0]=='H') {
		lf0 = nHelix(c0,r0,crd0);
	} else lf0 = nStrand(c0,r0,crd0);
	XYZ crd1;
	LocalFrame lf1;
	if(secstr[c1][r1]=='H') {
		lf1 = cHelix(c1,r1,crd1);
	} else lf1 = cStrand(c1,r1,crd1);

	double dis = sqrt((crd0-crd1).squarednorm());
	//std::cout <<dis <<std::endl;
	auto discut = DisCut();
	//for(auto &p1:discut) {
	//	std::cout <<p1.first <<' ' <<p1.second.size() <<std::endl;
	//	for(auto &p:p1.second) {
	//		std::cout <<p.first <<' ' <<p.second.size() <<' ' <<p.second[0] <<' ' <<p.second[1] <<std::endl;
	//	}
	//}
	std::map<int,std::vector<double>> dc = discut.at(key);
	std::vector<int> lth_allow;
	for(int i=minlth;i<=maxlth;i++) {
		if(dc.find(i)==dc.end()) continue;
		if(dis<dc.at(i)[0]) continue;
		if(dis>dc.at(i)[1]) continue;
		lth_allow.push_back(i);
	}
	//std::cout <<1 <<std::endl;
	std::vector<double> qs = AddElement(crd0,lf0,crd1,lf1);
	//for(double d:qs) std::cout <<' ' <<d;
	//std::cout <<std::endl;
	auto trees=BuildTree();
	//for(auto &tr:trees) {
	//	std::cout <<tr.first <<' ' <<tr.second.size() <<std::endl;
	//}
	auto ranges=Ranges();
	//for(double d:ranges) std::cout <<' ' <<d;
	//std::cout <<std::endl;
	for(int i:lth_allow) {
		std::pair<double,SearchTree::TREE> &p=trees.at(key).at(i);
		//std::cout <<1 <<std::endl;
		scs.insert({i,p.first*p.second.score(qs,ranges)});
		//std::cout <<2 <<std::endl;

	}
	return scs;
}

std::map<int,double> LoopLengthAndLocation::criterionfilter(int c0, int r0, int c1, int r1, std::map<int,double>&base) {
	std::string key{secstr[c0][r0],secstr[c1][r1]};
	std::map<int,double> newbase;
	auto criterion=Criterion();
	//for(auto &p1:criterion) {
	//	std::cout <<p1.first <<' ' <<p1.second.size() <<std::endl;
	//	for(auto &p:p1.second) {
	//		std::cout <<'\t' <<p.first <<' ' <<p.second <<std::endl;
	//	}
	//}
	//exit(1);
	for(auto &p:base) {
		//std::cout <<key <<' ' <<p.first <<std::endl;
		if(p.second<criterion.at(key).at(p.first)) continue;
		newbase.insert(p);
	}
	return newbase;
}
std::vector<std::vector<Residue>> LoopLengthAndLocation::looplengthtest(int c0, int r0, int c1, int r1, int lth) {
	int ntry=100;
	if(lth==2 || lth==1) ntry=10;
	std::vector<std::vector<Residue>> lps;
	int ntime=0;
	while(lps.empty() && ntime++<ntry) {
		lps = LOOP::MakeNLoop(chains, 1, c0, r0, c1, r1, lth);
	}
	return lps;
}
std::vector<std::pair<std::vector<int>,
			std::pair<double,std::vector<std::vector<Residue>>>>> LoopLengthAndLocation::linkfilter(
			int c0, int r0, int c1, int r1, bool test1more, std::map<int,double>&base) {
	std::vector<std::pair<std::vector<int>,std::pair<double,std::vector<std::vector<Residue>>>>> rts;
	std::vector<std::vector<Residue>> lps;
	for(auto &mp:base) {
		int lth = mp.first;
		lps = looplengthtest(c0,r0,c1,r1,lth);
		if(!lps.empty()) rts.push_back({{c0,r0,c1,r1,lth},{mp.second,lps}});
		lps.clear();
	}
	if(!test1more) return rts;
	for(auto &mp:base) {
		int lth = mp.first;
		lps = looplengthtest(c0,r0-1,c1,r1,lth+1);
		if(!lps.empty()) rts.push_back({{c0,r0-1,c1,r1,lth+1},{mp.second,lps}});
		lps.clear();
		lps = looplengthtest(c0,r0,c1,r1+1,lth+1);
		if(!lps.empty()) rts.push_back({{c0,r0,c1,r1+1,lth+1},{mp.second,lps}});
		lps.clear();
		lps = looplengthtest(c0,r0-1,c1,r1+1,lth+2);
		if(!lps.empty()) rts.push_back({{c0,r0-1,c1,r1+1,lth+2},{mp.second,lps}});
		lps.clear();
	}
	return rts;
}
std::map<std::vector<int>,std::pair<double,std::vector<std::vector<Residue>>>> LoopLengthAndLocation::lengthandlocation(
		std::vector<int>&sts, std::vector<int>&ens, int minlth, int maxlth, bool test1more) {
	std::map<std::vector<int>,std::pair<double,std::vector<std::vector<Residue>>>> rts;
	for(int i=sts[1];i<=sts[2];i++) {
		for(int j=ens[1];j<=ens[2];j++) {
			//std::cout <<1 <<std::endl;
			std::map<int,double> mid = score(sts[0],i,ens[0],j,minlth,maxlth);
			//std::cout <<2 <<std::endl;
			std::map<int,double> vi = criterionfilter(sts[0],i,ens[0],j,mid);
			//std::cout <<3 <<std::endl;
			auto lps = linkfilter(sts[0],i,ens[0],j,test1more,vi);
			//std::cout <<4 <<std::endl;
			for(auto &p:lps) {
				auto iter = rts.find(p.first);
				if(iter==rts.end()) rts.insert(p);
				else {
					for(auto &v:p.second.second) iter->second.second.push_back(v);
					if(p.second.first<iter->second.first) iter->second.first=p.second.first;
				}
			}
			//std::cout <<5 <<std::endl;
		}
	}
	return rts;
}










