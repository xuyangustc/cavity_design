/*
 * chemicalgroup.cpp
 *
 *  Created on: 2018年4月20日
 *      Author: hyliu
 */
#include "integrateddesign/chemicalgroup.h"
#include "dataio/datapaths.h"
#include "dataio/inputlines.h"
#include "geometry/calculators.h"
#include "sd/backboneff.h"
#include "sd/sidechainff.h"
#include "dstl/randomengine.h"
#include "geometry/rotation.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <cmath>
#define HBStaDis 5.0
using namespace NSPproteinrep;

std::map<std::string,std::string> SideChainChemGrp {
	{"TYR","phenol"},
	{"TRP","indole"},
	{"HIS","imidazole"},
	{"ARG","guanidyl"},
	{"LYS","amidogen"},
	{"SER","hydroxy"},
	{"THR","hydroxy"},
	{"ASP","carboxyl"},
	{"GLU","carboxyl"},
	{"ASN","amide"},
	{"GLN","amide"}
};

std::map<std::string,std::string> ResExampleForChemGrp {
	{"amide","ASN"},
	{"amidogen","LYS"},
	{"phenol","TYR"},
	{"carboxyl","ASP"},
	{"guanidyl","ARG"},
	{"hydroxy","THR"},
	{"imidazole","HIS"},
	{"indole","TRP"}
};

std::map<std::string,int> HBDonerAccepter {
	{"hydroxy_COH",0},
	{"carboxyl_COOH",-1},
	{"amide_CON_O",-1},
	{"amide_CON_N",1},
	{"amidogen_CNH3",1},
	{"guanidyl_CN3_NH",1},
	{"guanidyl_CN3_NH2",1},
	{"imidazole_CNCNC",0},
	{"phenol_C6OH",0},
	{"indole_CCNC6",1},
	{"peptide_CCONC_O",-1},
	{"peptide_CCONC_N",1}
};

std::map<std::string,std::string> UnitInChemGrp {
	{"hydroxy_COH","hydroxy"},
	{"carboxyl_COOH","carboxyl"},
	{"amide_CON_O","amide"},
	{"amide_CON_N","amide"},
	{"amidogen_CNH3","amidogen"},
	{"guanidyl_CN3_NH","guanidyl"},
	{"guanidyl_CN3_NH2","guanidyl"},
	{"imidazole_CNCNC","imidazole"},
	{"phenol_C6OH","phenol"},
	{"indole_CCNC6","indole"},
	{"peptide_CCONC_O","peptide_bond_O"},
	{"peptide_CCONC_N","peptide_bond_N"}
};

int InteractionUnit::getunittype(const std::string &unitname,
		NSPdataio::InputLines *il, int *lidx){
	static std::map<std::string,InteractionUnit>intunits;
	if(il){
		std::vector<std::string> words=(*il)[(*lidx)++];
		int nunits=std::stoi(words[0]);
		for(int n=0;n<nunits;++n){
			words =(*il)[(*lidx)++];
			int ut=std::stoi(words[1]);
			intunits.insert(std::make_pair(words[0],InteractionUnit(words[0],ut)));
		}
		return nunits;
	} else {
		return intunits.at(unitname).unittype;
	}
}
const ChemGrp & ChemGrp::getchemgrp(const std::string &name){
	static std::map<std::string,ChemGrp> chemgrps;
	static bool read{false};
	if(!read){
#ifdef _OPENMP
#pragma omp critical (readchemgrps)
		{
			if(!read){
#endif
		readchemgrps(chemgrps);
		read=true;
#ifdef _OPENMP
			}
		}
#endif
	}
	return chemgrps.at(name);
}
void ChemGrp::readchemgrps(std::map<std::string,ChemGrp> &chemgrps,
		const std::string &filename){
		std::string file=NSPdataio::datafilename(filename);
		NSPdataio::InputLines inputlines;
		inputlines.init(file,'#');
		int lidx=0;
		InteractionUnit::getunittype("", &inputlines,&lidx);
		std::string grpname;
		while (lidx<inputlines.size()){
			std::vector<std::string> &words=inputlines[lidx++];
			if(words[0]=="END"){
				continue;
			} else if(words[0]=="START" ){
				grpname=words[1];
				chemgrps.insert(std::make_pair(grpname,ChemGrp()));
			} else {
				ChemGrp &grp=chemgrps.at(grpname);
				if(words[0]=="UNIT") {
					UnitInGrp u;
					u.unitname=words[1];
					for(int a=0;a<3;++a){
						u.atoms.push_back(std::stoi(words[a+2]));
					}
					grp.units.push_back(u);
				} else {
					std::string modulename=words[0];
					std::vector<std::string> atoms;
					for(int a=1;a<words.size();++a){
						atoms.push_back(words[a]);
					}
					grp.modulesandatoms.insert(std::make_pair(modulename,atoms));
				}
			}
		}
}

std::vector<XYZ> buildgrpinmc(double b1, double b2, double ag, std::vector<int> &rk) {
	std::vector<XYZ> crds(3,XYZ(0.0,0.0,0.0));
	crds[1]=InternaltoXYZ(crds[0],b1);
	crds[2]=InternaltoXYZ(crds[1],crds[0],b2,ag);
	std::vector<XYZ> crds1;
	for(int i:rk) {
		crds1.push_back(crds[i]);
	}
	return crds1;
}

std::vector<XYZ> buildgrpinsc(std::string &res,
		std::vector<std::string> &atoms, std::vector<XYZ> csmc) {
	std::map<std::string,XYZ> atomcrd;
	auto &vsc=NSPsd::VSCType::getVSCType(res);
	std::vector<XYZ> r(4+vsc.nscatoms);
	r[0]=csmc[0];
	r[1]=csmc[1];
	r[2+vsc.nscatoms]=csmc[2];
	r[3+vsc.nscatoms]=csmc[3];
	atomcrd.insert({"N",r[0]});
	atomcrd.insert({"CA",r[1]});
	atomcrd.insert({"C",r[2+vsc.nscatoms]});
	atomcrd.insert({"O",r[3+vsc.nscatoms]});
	auto internalcrds=vsc.internalcrds;
	double deg=3.14159265/180.0;
	for(int i=0;i<vsc.nscatoms;++i){
		auto &ic=internalcrds[i];
		r[i+2]=InternaltoXYZ(r[ic[0].first],r[ic[1].first],
				r[ic[2].first],ic[0].second, ic[1].second*deg,ic[2].second*deg);
		atomcrd.insert({vsc.atomnames[i],r[i+2]});
	}
	std::vector<XYZ> crds;
	for(auto &a:atoms) {
		if(atomcrd.find(a)==atomcrd.end()) {
			std::cout <<"Chemical Groups are not corresponding to topology of residue!" <<std::endl;
			exit(1);
		}
		crds.push_back(atomcrd.at(a));
	}
	return crds;
}

std::vector<XYZ> ChemGrp::buildcrd(std::string grpname) {
	static std::map<std::string,std::vector<XYZ>> crds;
	if(crds.empty()) {
		std::string gpn;
		std::vector<int> pb;
		std::vector<XYZ> cs;
		double deg=3.14159265/180.0;
		double boc=backboneff.b0_co;
		double bcca=backboneff.b0_cac;
		double aocca=backboneff.t0_caco*deg;
		double bcn=backboneff.b0_cn;
		double bnca=backboneff.b0_nca;
		double acnca=backboneff.t0_cnca*deg;
		double ancac=backboneff.t0_ncac*deg;
		gpn="peptide_bond_N";
		pb={1,0,2};
		cs=buildgrpinmc(bcn,bnca,acnca,pb);
		crds.insert({gpn,cs});
		gpn="peptide_bond_O";
		pb={0,1,2};
		cs=buildgrpinmc(boc,bcca,aocca,pb);
		crds.insert({gpn,cs});
		std::vector<XYZ> csmc1=cs;//o,c,ca,n
		XYZ c_n=InternaltoXYZ(cs[2],cs[1],bnca,ancac);
		csmc1.push_back(c_n);
		std::vector<XYZ> csmc{csmc1[3],csmc1[2],csmc1[1],csmc1[0]};
		for(auto &gs:ResExampleForChemGrp) {
			std::vector<std::string> ans=getchemgrp(gs.first).modulesandatoms.at(gs.second);
			cs=buildgrpinsc(gs.second,ans,csmc);
			crds.insert({gs.first,cs});
		}
	}
	if(crds.find(grpname)==crds.end()) {
		std::cout <<"Chemical group which is not defined!" <<std::endl;
		exit(1);
	}
	return crds.at(grpname);
}

bool UnitPairInGrpPair::findunitpair(GrpInProtein &gip1, GrpInProtein &gip2, int &n1, int &n2) {
	auto &us1=gip1.grp->units;
	auto &us2=gip2.grp->units;
	bool sat{false};
	double dx=10000.0;
	for(int i1=0;i1<us1.size();i1++) {
		for(int i2=0;i2<us2.size();i2++) {
			int da1=InteractionUnit::getunittype(us1[i1].unitname);
			int da2=InteractionUnit::getunittype(us2[i2].unitname);
			if (da1==1 && da2==1) continue;
			if (da1==-1 && da2==-1) continue;
			int ix1=us1[i1].atoms[0];
			int ix2=us2[i2].atoms[0];
			XYZ &c1=gip1.crd[ix1];
			XYZ &c2=gip2.crd[ix2];
			double dx1=sqrt((c1-c2).squarednorm());
			if(dx1>HBStaDis) continue;
			sat=true;
			if(dx1<dx) {
				dx=dx1;
				n1=i1;
				n2=i2;
			}
		}
	}
	return sat;
}

bool UnitPairInGrpPair::findunitpair(GrpInProtein &gip1, GrpInProtein &gip2, std::string &str1, std::string &str2) {
	auto &us1=gip1.grp->units;
	auto &us2=gip2.grp->units;
	bool sat{false};
	double dx=10000.0;
	for(int i1=0;i1<us1.size();i1++) {
		for(int i2=0;i2<us2.size();i2++) {
			int da1=InteractionUnit::getunittype(us1[i1].unitname);
			int da2=InteractionUnit::getunittype(us2[i2].unitname);
			if (da1==1 && da2==1) continue;
			if (da1==-1 && da2==-1) continue;
			int ix1=us1[i1].atoms[0];
			int ix2=us2[i2].atoms[0];
			XYZ &c1=gip1.crd[ix1];
			XYZ &c2=gip2.crd[ix2];
			double dx1=sqrt((c1-c2).squarednorm());
			if(dx1>HBStaDis) continue;
			sat=true;
			if(dx1<dx) {
				dx=dx1;
				str1=us1[i1].unitname;
				str2=us2[i2].unitname;
			}
		}
	}
	return sat;
}

void UnitPairInGrpPair::getgrppaircrd(GrpInProtein &gip1, GrpInProtein &gip2, AtomsPairCrd &apc) {
	int n1=-1,n2=-1;
	if(!findunitpair(gip1,gip2,n1,n2)) return;
	std::string un1=gip1.grp->units[n1].unitname;
	std::string un2=gip2.grp->units[n2].unitname;
	std::pair<std::string,std::string> up=std::make_pair(un1,un2);
	if(un2<un1) up=std::make_pair(un2,un1);
	if(apc.find(up)==apc.end())
		apc.insert({up,std::vector<std::pair<std::vector<XYZ>,std::vector<XYZ>>>()});
	if(un2<un1) apc.at(up).push_back({gip2.crd,gip1.crd});
	else apc.at(up).push_back({gip1.crd,gip2.crd});
}

void UnitPairInGrpPair::getunitpaircrd(GrpInProtein &gip1, GrpInProtein &gip2, AtomsPairCrd &apc) {
	int n1=-1,n2=-1;
	if(!findunitpair(gip1,gip2,n1,n2)) return;
	std::string un1=gip1.grp->units[n1].unitname;
	std::string un2=gip2.grp->units[n2].unitname;
	std::pair<std::string,std::string> up=std::make_pair(un1,un2);
	if(un2<un1) up=std::make_pair(un2,un1);
	if(apc.find(up)==apc.end())
		apc.insert({up,std::vector<std::pair<std::vector<XYZ>,std::vector<XYZ>>>()});
	std::vector<int> ix1=gip1.grp->units[n1].atoms;
	std::vector<int> ix2=gip2.grp->units[n2].atoms;
	std::vector<XYZ> cs1,cs2;
	for(int i=0;i<ix1.size();i++) cs1.push_back(gip1.crd[ix1[i]]);
	for(int i=0;i<ix2.size();i++) cs2.push_back(gip2.crd[ix2[i]]);
	if(un2<un1) apc.at(up).push_back({cs2,cs1});
	else apc.at(up).push_back({cs1,cs2});
}

std::vector<double> UnitPairInGrpPair::upmatrix(std::vector<XYZ> &cs1, std::vector<XYZ> &cs2) {
	std::vector<double> um;
	for(XYZ &c1:cs1) {
		for(XYZ &c2:cs2) {
			um.push_back(sqrt((c2-c1).squarednorm()));
		}
	}
	assert(um.size()==9);
	assert(um[0]<=HBStaDis);
	return um;
}

void UnitPairInGrpPair::getunitpairmatrix(GrpInProtein gip1, GrpInProtein gip2, UnitPairMatrix &upm) {
	int n1=-1,n2=-1;
	if(!findunitpair(gip1,gip2,n1,n2)) return;
	std::string un1=gip1.grp->units[n1].unitname;
	std::string un2=gip2.grp->units[n2].unitname;
	std::pair<std::string,std::string> up=std::make_pair(un1,un2);
	if(un2<un1) up=std::make_pair(un2,un1);
	if(upm.find(up)==upm.end())
		upm.insert({up,std::vector<std::vector<double>>()});
	std::vector<int> ix1=gip1.grp->units[n1].atoms;
	std::vector<int> ix2=gip2.grp->units[n2].atoms;
	std::vector<XYZ> cs1,cs2;
	for(int i=0;i<ix1.size();i++) cs1.push_back(gip1.crd[ix1[i]]);
	for(int i=0;i<ix2.size();i++) cs2.push_back(gip2.crd[ix2[i]]);
	if(un2<un1) upm.at(up).push_back(upmatrix(cs2,cs1));
	else upm.at(up).push_back(upmatrix(cs1,cs2));
}

void UnitPairInGrpPair::randomunitpair(GrpInProtein &gip1, GrpInProtein &gip2, std::vector<ChemGrpCrd> &cgcs, int ntimes) {
	//std::set<std::pair<std::string,std::string>> ups=unitpairingrppair(gip1,gip2);
	int n1,n2;
	if(!findunitpair(gip1,gip2,n1,n2)) return;
	std::pair<std::string,std::string> up={gip1.grp->units[n1].unitname,gip2.grp->units[n2].unitname};
	/*std::vector<XYZ> cs1,cs2;
	std::vector<int> i1=gip1.grp->units[n1].atoms;
	std::vector<int> i2=gip2.grp->units[n2].atoms;
	for(int i=0;i<i1.size();i++) cs1.push_back(gip1.crd[i1[i]]);
	for(int i=0;i<i2.size();i++) cs2.push_back(gip2.crd[i2[i]]);*/
	ChemGrpCrd &cs1=gip1.crd;
	ChemGrpCrd &cs2=gip2.crd;
	std::set<int> s1,s2;
	for(int i=0;i<gip1.grp->units.size();i++) {
		if (gip1.grp->units[i].unitname==gip1.grp->units[n1].unitname) {
			s1.insert(gip1.grp->units[i].atoms[0]);
		}
	}
	for(int i=0;i<gip2.grp->units.size();i++) {
		if (gip2.grp->units[i].unitname==gip2.grp->units[n2].unitname) {
			s2.insert(gip2.grp->units[i].atoms[0]);
		}
	}
	int nx=0;
	auto & rng = NSPdstl::RandomEngine<>::getinstance().realrng(0.0,1.0);
	auto & rngi = NSPdstl::RandomEngine<>::getinstance();
	while(nx<ntimes) {
		QuaternionCrd qc(rng,0);
		Rotation rt(qc,cs2[0]);
		for(auto &c:cs2) rt.apply(&c);
		int j1=rngi.intrng(0,s1.size()-1)();
		int j2=rngi.intrng(0,s2.size()-1)();
		XYZ rand(rng,HBStaDis);
		XYZ dis=cs1[j1]+rand-cs2[j2];
		for(auto &c:cs2) c=c+dis;
		assert((cs2[j2]-cs1[j1]).squarednorm()<HBStaDis*HBStaDis+0.01);
		int m1,m2;
		if(!findunitpair(gip1,gip2,m1,m2)) continue;
		int k1=gip1.grp->units[m1].atoms[0];
		int k2=gip2.grp->units[m2].atoms[0];
		if(s1.find(k1)==s1.end()) continue;
		if(s2.find(k2)==s2.end()) continue;
		cgcs.push_back(cs2);
		nx++;
	}
}

int ChemGrpsInProtein::mcunitinaa(FullSite &site, std::string chemgrpname, std::string molname) {
	GrpInProtein gip(0,0,chemgrpname,molname,true);
	auto it=gip.grp->modulesandatoms.find(molname);
	std::vector<std::string> atoms=it->second;
	bool success=true;
	for(auto &a:atoms) {
		if(site.hasatomcrd(a)) gip.crd.push_back(site.getcrd(a));
		else {
			success=false;
			break;
		}
	}
	if(success) {
		gips.push_back(gip);
		return 1;
	}
	return 0;
}

int ChemGrpsInProtein::mcunitlinkaa(FullSite &site,FullSite &siteb,std::vector<int> &pb,
		std::string chemgrpname,std::string molname) {
	GrpInProtein gip(0,0,chemgrpname,molname,true);
	auto it=gip.grp->modulesandatoms.find(molname);
	std::vector<std::string> atoms=it->second;
	for(int i=0;i<atoms.size();i++) {
		if(pb[i]==0) {
			if(site.hasatomcrd(atoms[i])) gip.crd.push_back(site.getcrd(atoms[i]));
			else return 0;
		} else {
			if(siteb.hasatomcrd(atoms[i])) gip.crd.push_back(siteb.getcrd(atoms[i]));
			else return 0;
		}
	}
	gips.push_back(gip);
	return 1;
}

int ChemGrpsInProtein::grpsmainchain(std::vector<FullSite> &chain,int posi){
	FullSite & site=chain[posi];
	bool ispro=(site.resname()=="PRO");
	int ngrp=0;
	std::string chemgrpname;
	std::string molname;
	if(chain.size()==1) {
		chemgrpname="amidogen";
		molname="terminalN";
		if(!ispro) ngrp+=mcunitinaa(site,chemgrpname,molname);
		chemgrpname="carboxyl";
		molname="terminalC";
		ngrp+=mcunitinaa(site,chemgrpname,molname);
	} else if(posi==0) {
		chemgrpname="amidogen";
		molname="terminalN";
		if(!ispro) ngrp+=mcunitinaa(site,chemgrpname,molname);

		chemgrpname="peptide_bond_O";
		molname="mainchain";
		ngrp+=mcunitinaa(site,chemgrpname,molname);
	} else {
		if(!ispro) {
			FullSite & sitep=chain[posi-1];
			chemgrpname="peptide_bond_N";
			molname="mainchain";
			std::vector<int> pb{0,1,1};
			ngrp+=mcunitlinkaa(sitep,site,pb,chemgrpname,molname);
		}
		if(posi==chain.size()-1) {
			chemgrpname="carboxyl";
			molname="terminalC";
		} else {
			chemgrpname="peptide_bond_O";
			molname="mainchain";
		}
		ngrp+=mcunitinaa(site,chemgrpname,molname);
	}
	return ngrp;
}

int ChemGrpsInProtein::grpssidechain(FullSite &site) {
	std::string resname=site.resname();
	if(SideChainChemGrp.find(resname)==SideChainChemGrp.end()) return 0;
	std::string chemgrpname=SideChainChemGrp.at(resname);
	GrpInProtein gip(0,0,chemgrpname,resname);
	auto it=gip.grp->modulesandatoms.find(resname);
	std::vector<std::string> atoms=it->second;
	bool success=true;
	for(auto &a:atoms) {
		if(site.hasatomcrd(a))gip.crd.push_back(site.getcrd(a));
		else {
			success=false;
			break;
		}
	}
	if(success) {
		gips.push_back(gip);
		return 1;
	}
	return 0;
}

void ChemGrpsInProtein::init(std::vector<std::vector<FullSite>> & chains){
	for(int cid=0;cid<chains.size();++cid){
		for(int resid=0;resid<chains[cid].size();++resid){
			int ngrps=grpsmainchain(chains[cid],resid);
			ngrps+=grpssidechain(chains[cid][resid]);
			if(ngrps==0) continue;
			for(auto it=gips.end()-1;it>=gips.end()-ngrps;--it){
				it->chainid=cid;
				it->resid=resid;
			}
		}
	}
}

void ChemGrpsInProtein::getgrppaircrd(AtomsPairCrd &apc) {
	for(auto &gp1:gips) {
		for(auto &gp2:gips) {
			if(gp1.chainid>gp2.chainid) continue;
			if(gp1.chainid==gp2.chainid && gp1.resid>=gp2.resid-1) continue;
			if(gp1.ismainchain && gp2.ismainchain) continue;
			UnitPairInGrpPair::getgrppaircrd(gp1,gp2,apc);
		}
	}
}

void ChemGrpsInProtein::getunitpairmatrix(UnitPairMatrix &upm) {
	for(auto &gp1:gips) {
		for(auto &gp2:gips) {
			if(gp1.chainid>gp2.chainid) continue;
			if(gp1.chainid==gp2.chainid && gp1.resid>=gp2.resid-1) continue;
			if(gp1.ismainchain && gp2.ismainchain) continue;
			UnitPairInGrpPair::getunitpairmatrix(gp1,gp2,upm);
		}
	}
}

void ChemGrpsInProtein::getunitpaircrd(AtomsPairCrd &apc) {
	for(auto &gp1:gips) {
		for(auto &gp2:gips) {
			if(gp1.chainid>gp2.chainid) continue;
			if(gp1.chainid==gp2.chainid && gp1.resid>=gp2.resid-1) continue;
			if(gp1.ismainchain && gp2.ismainchain) continue;
			UnitPairInGrpPair::getunitpaircrd(gp1,gp2,apc);
		}
	}
}

void UnitPairInGrpPair::printgrppaircrd(AtomsPairCrd &apc, std::string dir) {
	if(dir.back()!='/') dir+='/';
	std::ofstream ofs;
	for(auto &ac:apc) {
		std::string fn=dir+ac.first.first+"__"+ac.first.second;
		ofs.open(fn,std::ofstream::app);
		for(auto &cs:ac.second) {
			for(XYZ &c:cs.first) ofs <<c.x_ <<' ' <<c.y_ <<' ' <<c.z_ <<' ';
			ofs <<std::endl;
			for(XYZ &c:cs.second) ofs <<c.x_ <<' ' <<c.y_ <<' ' <<c.z_ <<' ';
			ofs <<std::endl;
		}
		ofs.close();
	}
}

void UnitPairInGrpPair::printunitpaircrd(AtomsPairCrd &apc, std::string dir) {
	if(dir.back()!='/') dir+='/';
	std::ofstream ofs;
	for(auto &ac:apc) {
		std::string fn=dir+ac.first.first+"__"+ac.first.second;
		ofs.open(fn,std::ofstream::app);
		for(auto &cs:ac.second) {
			for(XYZ &c:cs.first)
				ofs <<std::fixed <<std::setprecision(3) <<c.x_ <<'\t'
					<<std::fixed <<std::setprecision(3) <<c.y_ <<'\t'
					<<std::fixed <<std::setprecision(3) <<c.z_ <<'\t';
			for(XYZ &c:cs.second)
				ofs <<std::fixed <<std::setprecision(3) <<c.x_ <<'\t'
					<<std::fixed <<std::setprecision(3) <<c.y_ <<'\t'
					<<std::fixed <<std::setprecision(3) <<c.z_ <<'\t';
			ofs <<std::endl;
		}
		ofs.close();
	}
}

void UnitPairInGrpPair::printunitpairmatrix(UnitPairMatrix &upm, std::string dir) {
	if(dir.back()!='/') dir+='/';
	std::ofstream ofs;
	for(auto &um:upm) {
		std::string fn=dir+um.first.first+"__"+um.first.second;
		ofs.open(fn,std::ofstream::app);
		for(auto &vd:um.second) {
			for(double &d:vd) ofs <<d <<' ';
			ofs <<std::endl;
		}
		ofs.close();
	}
}




















std::vector<GrpInProtein> NSPproteinrep::grpsinpdb(
	std::vector<std::vector<FullSite>> & chains){
	std::vector<GrpInProtein> result;
	for(int cid=0;cid<chains.size();++cid){
		for(int resid=0;resid<chains[cid].size();++resid){
			int ngrps=grpsmainchain(chains[cid],resid,result);
			ngrps+=grpssidechain(chains[cid][resid],result);
			if(ngrps==0) continue;
			for(auto it=result.end()-1;it>=result.end()-ngrps;--it){
				it->chainid=cid;
				it->resid=resid;
			}
		}
	}
	return result;
}

std::vector<GrpInProtein> NSPproteinrep::grpsinpdb(std::string pdbfile) {
	std::vector<std::vector<FullSite>> chains=readfullsitesfrompdb(pdbfile,true);
	return grpsinpdb(chains);
}

int NSPproteinrep::mcunitinaa(FullSite &site,std::vector<GrpInProtein> &grps,
		std::string chemgrpname,std::string molname) {
	GrpInProtein gip(0,0,chemgrpname,molname,true);
	auto it=gip.grp->modulesandatoms.find(molname);
	std::vector<std::string> atoms=it->second;
	bool success=true;
	for(auto &a:atoms) {
		if(site.hasatomcrd(a)) gip.crd.push_back(site.getcrd(a));
		else {
			success=false;
			break;
		}
	}
	if(success) {
		grps.push_back(gip);
		return 1;
	}
	return 0;
}

int NSPproteinrep::mcunitlinkaa(FullSite &site,FullSite &siteb,std::vector<int> &pb,
		std::vector<GrpInProtein> &grps,std::string chemgrpname,std::string molname) {
	GrpInProtein gip(0,0,chemgrpname,molname,true);
	auto it=gip.grp->modulesandatoms.find(molname);
	std::vector<std::string> atoms=it->second;
	for(int i=0;i<atoms.size();i++) {
		if(pb[i]==0) {
			if(site.hasatomcrd(atoms[i])) gip.crd.push_back(site.getcrd(atoms[i]));
			else return 0;
		} else {
			if(siteb.hasatomcrd(atoms[i])) gip.crd.push_back(siteb.getcrd(atoms[i]));
			else return 0;
		}
	}
	grps.push_back(gip);
	return 1;
}

int NSPproteinrep::grpsmainchain(std::vector<FullSite> &chain,int posi,
		std::vector<GrpInProtein> &grps){
	FullSite & site=chain[posi];
	bool ispro=(site.resname()=="PRO");
	int ngrp=0;
	std::string chemgrpname;
	std::string molname;
	if(chain.size()==1) {
		chemgrpname="amidogen";
		molname="terminalN";
		if(!ispro) ngrp+=mcunitinaa(site,grps,chemgrpname,molname);
		chemgrpname="carboxyl";
		molname="terminalC";
		ngrp+=mcunitinaa(site,grps,chemgrpname,molname);
	} else if(posi==0) {
		chemgrpname="amidogen";
		molname="terminalN";
		if(!ispro) ngrp+=mcunitinaa(site,grps,chemgrpname,molname);

		chemgrpname="peptide_bond_O";
		molname="mainchain";
		ngrp+=mcunitinaa(site,grps,chemgrpname,molname);
	} else {
		if(!ispro) {
			FullSite & sitep=chain[posi-1];
			chemgrpname="peptide_bond_N";
			molname="mainchain";
			std::vector<int> pb{1,0,1};
			ngrp+=mcunitlinkaa(sitep,site,pb,grps,chemgrpname,molname);
		}
		if(posi==chain.size()-1) {
			chemgrpname="carboxyl";
			molname="terminalC";
		} else {
			chemgrpname="peptide_bond_O";
			molname="mainchain";
		}
		ngrp+=mcunitinaa(site,grps,chemgrpname,molname);
	}
	return ngrp;
}

int NSPproteinrep::grpssidechain(FullSite &site,std::vector<GrpInProtein> &grps) {
	std::string resname=site.resname();
	if(SideChainChemGrp.find(resname)==SideChainChemGrp.end()) return 0;
	std::string chemgrpname=SideChainChemGrp.at(resname);
	GrpInProtein gip(0,0,chemgrpname,resname);
	auto it=gip.grp->modulesandatoms.find(resname);
	std::vector<std::string> atoms=it->second;
	bool success=true;
	for(auto &a:atoms) {
		if(site.hasatomcrd(a))gip.crd.push_back(site.getcrd(a));
		else {
			success=false;
			break;
		}
	}
	if(success) {
		grps.push_back(gip);
		return 1;
	}
	return 0;
}

OneGrpInfo GrpInProtein::unitsingroup() {
	OneGrpInfo ogi;
	auto &uits=grp->units;
	std::vector<std::string> atoms=grp->modulesandatoms.at(molname);
	for(auto &ut:uits) {
		UnitInfo ui;
		ui.first=ut.unitname;
		std::vector<int> ix=ut.atoms;
		for(int i=0;i<ix.size();i++) {
			ui.second.push_back({atoms[ix[i]],crd[ix[i]]});
		}
		ogi.push_back(ui);
	}
	return ogi;
}

void NSPproteinrep::atom3inpdb(std::vector<GrpInProtein> &gips, GrpsInfo &gis) {
	for(auto &gp:gips) {
		if(gis.find(gp.grpname)==gis.end())
			gis.insert({gp.grpname,std::vector<OneGrpInfo>()});
		gis.at(gp.grpname).push_back(gp.unitsingroup());
	}
}

void NSPproteinrep::atom3inpdb(std::vector<std::vector<FullSite>> & chains, GrpsInfo &gis) {
	std::vector<GrpInProtein> gips=grpsinpdb(chains);
	atom3inpdb(gips,gis);
}
void NSPproteinrep::atom3inpdb(std::string pdbfile, GrpsInfo &gis) {
	std::vector<std::vector<FullSite>> chains=readfullsitesfrompdb(pdbfile,true);
	atom3inpdb(chains,gis);
}

void NSPproteinrep::getunitpaircrd(std::vector<GrpInProtein> &gips, UnitPairCrd &upcs) {
	for(auto &gp1:gips) {
		OneGrpInfo og1=gp1.unitsingroup();
		for(auto &gp2:gips) {
			if(gp1.chainid>gp2.chainid) continue;
			if(gp1.chainid==gp2.chainid && gp1.resid>=gp2.resid-1) continue;
			if(gp1.ismainchain && gp2.ismainchain) continue;
			OneGrpInfo og2=gp2.unitsingroup();
			bool sat{false};
			int n1=-1,n2=-1;
			double dx=10000.0;
			for(int i1=0;i1<og1.size();i1++) {
				for(int i2=0;i2<og2.size();i2++) {
					int da1=InteractionUnit::getunittype(og1[i1].first);
					int da2=InteractionUnit::getunittype(og2[i2].first);
					if (da1==1 && da2==1) continue;
					if (da1==-1 && da2==-1) continue;
					XYZ &c1=og1[i1].second[0].second;
					XYZ &c2=og2[i2].second[0].second;
					double dx1=sqrt((c1-c2).squarednorm());
					if(dx1>HBStaDis) continue;
					sat=true;
					if(dx1<dx) {
						dx=dx1;
						n1=i1;
						n2=i2;
					}
				}
			}
			if(!sat) continue;
			std::string &un1=og1[n1].first;
			std::string &un2=og2[n2].first;
			std::pair<std::string,std::string> up=std::make_pair(un1,un2);
			if(un2<un1) up=std::make_pair(un2,un1);
			if(upcs.find(up)==upcs.end())
				upcs.insert({up,std::vector<std::vector<XYZ>>()});
			std::vector<XYZ> cs;
			if(un2<un1) {
				for(int i=0;i<gp2.crd.size();i++) cs.push_back(gp2.crd[i]);
				for(int i=0;i<gp1.crd.size();i++) cs.push_back(gp1.crd[i]);
			} else {
				for(int i=0;i<gp1.crd.size();i++) cs.push_back(gp1.crd[i]);
				for(int i=0;i<gp2.crd.size();i++) cs.push_back(gp2.crd[i]);
			}
			upcs.at(up).push_back(cs);
		}
	}
}

bool NSPproteinrep::getupmx(UnitInfo &u1, UnitInfo &u2, std::vector<double> &ds) {
	int i1=InteractionUnit::getunittype(u1.first);
	int i2=InteractionUnit::getunittype(u2.first);
	assert(i1==-1||i1==0||i1==1);
	assert(i2==-1||i2==0||i2==1);
	if(i1==1 && i2==1) return false;
	if(i1==-1 && i2==-1) return false;
	for(int i=0;i<u1.second.size();i++) {
		for(int j=0;j<u2.second.size();j++) {
			XYZ &c1=u1.second[i].second;
			XYZ &c2=u2.second[j].second;
			ds.push_back(sqrt((c1-c2).squarednorm()));
			if(ds[0]>HBStaDis) return false;
		}
	}
	assert(ds.size()==9);
	return true;
}

void NSPproteinrep::getunitpairmatrix(std::vector<GrpInProtein> &gips, UnitPairMatrix &upms) {
	std::vector<OneGrpInfo> ogis;
	for(auto &gip:gips) ogis.push_back(gip.unitsingroup());
	for(int i=0;i<ogis.size();i++) {
		for(int j=i+1;j<ogis.size();j++) {
			std::string u1,u2;
			std::vector<double> ms;
			bool sat{false};
			for(int g=0;g<ogis[i].size();g++) {
				for(int h=0;h<ogis[j].size();h++) {
					std::string &ig=ogis[i][g].first;
					std::string &jh=ogis[j][h].first;
					std::vector<double> ds;
					bool getted;
					if(ig<jh) getted=getupmx(ogis[i][g],ogis[j][h],ds);
					else getted=getupmx(ogis[j][h],ogis[i][g],ds);
					if(!getted) continue;
					sat=true;
					if(ms.empty() || ds[0]<ms[0]) {
						ms=ds;
						if(ig<jh) {
							u1=ig;
							u2=jh;
						} else {
							u1=jh;
							u2=ig;
						}
					}
				}
			}
			if(!sat) continue;
			std::pair<std::string,std::string> up=std::make_pair(u1,u2);
			if(upms.find(up)==upms.end())
				upms.insert({up,std::vector<std::vector<double>>()});
			upms.at(up).push_back(ms);
		}
	}
}

void NSPproteinrep::printupmatrix(UnitPairMatrix &upms, std::ostream &os) {
	for(auto &up:upms) {
		os <<up.first.first <<' ' <<up.first.second <<' ' <<up.second.size() <<std::endl;
		for(auto &vd:up.second) {
			for(double &d:vd) os <<d <<' ';
			os <<std::endl;
		}
	}
}

/*
std::map<std::string,std::vector<std::vector<UnitCrd>>> NSPproteinrep::atom3inpdb(
		std::vector<GrpInProtein> &gips) {
	std::map<std::string,std::vector<std::vector<UnitCrd>>> mods;
	for(auto &gp:gips) {
		if(mods.find(gp.grpname)==mods.end())
			mods.insert({gp.grpname,std::vector<std::vector<UnitCrd>>()});
		mods.at(gp.grpname).push_back(gp.unitsingroup());
	}
	return mods;
}

std::map<std::string,std::vector<std::vector<UnitCrd>>> NSPproteinrep::atom3inpdb(
		std::vector<std::vector<FullSite>> & chains) {
	std::vector<GrpInProtein> gips=grpsinpdb(chains);
	return atom3inpdb(gips);
}
std::map<std::string,std::vector<std::vector<UnitCrd>>> NSPproteinrep::atom3inpdb(
		std::string pdbfile) {
	std::vector<std::vector<FullSite>> chains=readfullsitesfrompdb(pdbfile,true);
	return atom3inpdb(chains);
}*/
/*
UnitPair NSPproteinrep::makeunitpair(const GrpInProtein &gip1, const GrpInProtein &gip2){
	if(gip1.chainid==gip2.chainid){
		if(gip1.resid==gip2.resid) return UnitPair();
		if(gip1.ismainchain ){
			if(gip2.ismainchain || gip1.resid==gip2.resid+1) return UnitPair();
		} else {
			if(gip2.ismainchain && gip2.resid==gip1.resid+1) return UnitPair();
		}
	}
	return UnitPair(*(gip1.grp),gip1.crd,*(gip2.grp),gip2.crd);
}*/
