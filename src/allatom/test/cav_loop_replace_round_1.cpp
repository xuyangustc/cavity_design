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

std::vector<XYZ> CaCrd(std::vector<Residue>&ch1, int st, int en) {
	std::vector<XYZ> cas;
	for(int i=0;i<6;i++) {
		cas.push_back(ch1[st+i].rds().at("CA").crd);
		cas.push_back(ch1[en-i].rds().at("CA").crd);
	}
	return cas;
}

std::vector<std::vector<Residue>> ReadHTHLib(std::string listfile, std::string dir) {
	//std::cout <<"HTH reading..." <<std::endl;
	if(dir.back()!='/') dir += '/';
	std::vector<std::string> fs =ReadLine(listfile);
	//std::cout <<fs.size() <<std::endl;
	std::vector<std::vector<Residue>> fgs(fs.size());
	for(int i=0;i<fs.size();i++) {
		//std::cout <<i <<std::endl;
		std::vector<std::vector<Residue>> chs=ReadPDB(dir+fs[i],true);
		int sz = chs[0].size()-1;
		for(int j=0;j<6;j++) {
			chs[0][j].changeaa("LEU");
			chs[0][sz-j].changeaa("LEU");
		}
		for(int j=6;j<=sz-6;j++) {
			std::string rn=chs[0][j].resname();
			if(rn=="PRO") continue;
			else if(rn=="GLY") chs[0][j].changeaa("ALA");
			else chs[0][j].changeaa("GLY");
		}
		fgs[i]=chs[0];
		//std::cout <<fgs[i].size() <<std::endl;
	}
	return fgs;
}

std::vector<std::vector<int>> Findsten(std::string ssfile) {
	std::vector<std::string> ls = ReadLine(ssfile);
	std::vector<std::vector<int>> vv;
	ls[1].push_back('C');
	int st, en;
	for(int j=0;j<ls[1].size();j++) {
		if(ls[1][j]!='H') continue;
		st=j;
		for(j++;j<ls[1].size();j++) {
			if(ls[1][j]=='H') continue;
			en=j-1;
			break;
		}
		if(en-st<7) continue;
		vv.push_back({st,en});
		//std::cout <<st <<' ' <<en <<std::endl;
	}
	//exit(1);
	std::vector<std::vector<int>> hth;
	for(int i=1;i<vv.size();i++) {
		hth.push_back({vv[i-1][1]-5, vv[i][0]+5});
	}
	return hth;
}

std::vector<std::vector<XYZ>> ExtractCrd(std::vector<std::vector<Residue>>&lps) {
	std::vector<std::vector<XYZ>> crd;
	for(int j=0;j<lps.size();j++) {
		std::vector<XYZ> cas;
		int sz=lps[j].size()-1;
		for(int i=0;i<6;i++) {
			cas.push_back(lps[j][i].rds().at("CA").crd);
			cas.push_back(lps[j][sz-i].rds().at("CA").crd);
		}
		crd.push_back(cas);
	}
	return crd;
}

std::vector<int> Halign(std::vector<std::vector<XYZ>>&lps,
		std::vector<Residue>&ch1, int st, int en) {
	std::vector<XYZ> cas;
	for(int i=0;i<6;i++) {
		cas.push_back(ch1[st+i].rds().at("CA").crd);
		cas.push_back(ch1[en-i].rds().at("CA").crd);
	}
	std::vector<int> rt;
	for(int i=0;i<lps.size();i++) {
		if(RMSD(cas,lps[i])<0.8) rt.push_back(i);
	}
	return rt;
}

void test(std::vector<std::vector<XYZ>>&lps, std::vector<std::vector<Residue>>&lpsr,
		std::vector<Residue>&ch1, int st, int en) {
	std::vector<XYZ> cas;
	for(int i=0;i<6;i++) {
		cas.push_back(ch1[st+i].rds().at("CA").crd);
		cas.push_back(ch1[en-i].rds().at("CA").crd);
	}
	for(int i=0;i<lps.size();i++) {
		NSPgeometry::QuatFit qf;
		double r2=qf.setup(cas,lps[i]);
		if(sqrt(r2)>0.8) continue;
		RigidTransform rt =qf.getRigidTransform();
		for(Residue &r:lpsr[i]) {
			for(auto &p:r.rds()) rt.apply(&(p.second.crd));
		}
		std::vector<std::vector<Residue>> ch1{lpsr[i]};
		residue2pdb("out",ch1);
		exit(1);
	}
}

int main(int argc, char ** argv) {
	NSPdstl::RandomEngine<>::getinstance().reseed(5918);
	std::vector<std::vector<Residue>> chs=ReadPDB(std::string(argv[1]), true);
	std::vector<std::vector<int>> sten=Findsten(std::string(argv[2]));
	std::vector<std::vector<Residue>> lps=ReadHTHLib(std::string(argv[3]),std::string(argv[4]));
	std::vector<std::vector<XYZ>> lpcrd=ExtractCrd(lps);
	std::vector<std::vector<int>> rts;
	for(int i=0;i<sten.size();i++) {
		std::vector<int> rt=Halign(lpcrd, chs[0], sten[i][0], sten[i][1]);
		std::cout <<rt.size() <<' ' <<sten.size() <<std::endl;
		if(rt.empty()) return 1;
		rts.push_back(rt);
		//test(lpcrd, lps, chs[0], sten[i][0], sten[i][1]);
	}
	std::set<std::vector<int>> combs;
	for(int i=0;i<10;i++) {
		std::vector<int> v;
		for(int j=0;j<rts.size();j++) {
			int c=NSPdstl::RandomEngine<>::getinstance().intrng(0,rts[j].size()-1)();
			v.push_back(c);
		}
		combs.insert(v);
	}
	//std::cout <<combs.size() <<std::endl;
	//for(auto &v:combs) {
	//	for(int i:v) std::cout <<i <<' ';
	//	std::cout <<std::endl;
	//}
	std::vector<std::vector<int>> cbs;
	for(auto &v:combs) cbs.push_back(v);
	std::string outpdb=std::string(argv[5]);
	std::string outres=std::string(argv[6]);
	std::string outlist=std::string(argv[7]);
	for(Residue &r:chs[0]) r.changeaa("LEU");
	std::ofstream ofs(outlist);
	for(int i=0;i<cbs.size();i++) {
		std::vector<std::vector<Residue>> vlp;
		for(int j=0;j<sten.size();j++) {
			std::vector<XYZ> c0 = CaCrd(chs[0],sten[j][0],sten[j][1]);
			std::vector<Residue> lp=lps[cbs[i][j]];
			std::vector<XYZ> c1 = CaCrd(lp,0,lp.size()-1);
			NSPgeometry::QuatFit qf;
			double r2=qf.setup(c0,c1);
			RigidTransform rt =qf.getRigidTransform();
			for(int e=0;e<lp.size();e++)
				for(auto &p:lp[e].rds()) rt.apply(&(p.second.crd));
			std::vector<Residue> lp1;
			for(int k=3;k<lp.size()-3;k++) lp1.push_back(lp[k]);
			vlp.push_back(lp1);
		}
		//std::cout <<2 <<std::endl;
		std::vector<std::vector<int>> loop_added;
		//std::cout <<chs[0].size() <<std::endl;
		std::vector<std::vector<int>> sten1;
		for(auto v:sten) sten1.push_back({0,v[0]+3,0,v[1]-3});
		std::vector<std::vector<Residue>> chs1=LOOP::NewChains(chs,vlp,sten1,loop_added);
		//std::cout <<1 <<std::endl;
		std::vector<int> gly;
		std::vector<int> pro;
		for(int j=0;j<chs1[0].size();j++) {
			if(chs1[0][j].resname()=="ALA") {
				chs1[0][j].changeaa("GLY");
				gly.push_back(j);
			}
			if(chs1[0][j].resname()=="ALA") pro.push_back(j);
		}
		residue2pdb(outpdb+"_"+std::to_string(i),chs1);
		std::ofstream ofs1(outres+"_"+std::to_string(i));
		ofs1 <<"default allButCys" <<std::endl;
		for(int k:gly) ofs1 <<"A " <<k <<" G" <<std::endl;
		for(int k:pro) ofs1 <<"A " <<k <<" P" <<std::endl;
		ofs1.close();
		ofs <<i <<std::endl;
	}
	ofs.close();
}

