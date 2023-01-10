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
	for(int i=st;i<=en;i++) {
		cas.push_back(ch1[i].rds().at("CA").crd);
	}
	return cas;
}

std::vector<std::vector<Residue>> ReadHTHLib(std::string listfile, std::string dir) {
	if(dir.back()!='/') dir += '/';
	std::vector<std::string> fs =ReadLine(listfile);
	std::vector<std::vector<Residue>> fgs(fs.size());
	for(int i=0;i<fs.size();i++) {
		std::vector<std::vector<Residue>> chs=ReadPDB(dir+fs[i],true);
		int sz = chs[0].size()-1;
		fgs[i]=chs[0];
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
	}
	std::vector<std::vector<int>> hth;
	for(int i=1;i<vv.size();i++) {
		hth.push_back({vv[i-1][1]-5, vv[i][0]+5});
	}
	return hth;
}

std::vector<std::vector<XYZ>> ExtractCrd(std::vector<std::vector<Residue>>&lps) {
	std::vector<std::vector<XYZ>> crd;
	for(int j=0;j<lps.size();j++) {
		crd.push_back(CaCrd(lps[j],0,lps[j].size()-1));
	}
	return crd;
}

bool Halign(std::vector<std::vector<XYZ>>&lps,
		std::vector<Residue>&ch1, int st, int en) {
	//std::cout <<st <<' ' <<en <<std::endl;
	std::vector<XYZ> cas=CaCrd(ch1,st,en);
	//std::cout <<cas.size() <<std::endl;
	std::vector<XYZ> cas3;
	for(int i=3;i<cas.size()-3;i++) cas3.push_back(cas[i]);
	std::vector<int> rt;
	//std::cout <<lps.size() <<std::endl;
	for(int i=0;i<lps.size();i++) {
		if(lps[i].size()!=cas.size()) continue;
		//std::cout <<lps[i].size() <<std::endl;
		//std::cout <<RMSD(cas,lps[i]) <<std::endl;
		if(RMSD(cas,lps[i])>0.8) continue;
		std::vector<XYZ> c3;
		for(int j=3;j<lps[i].size()-3;j++) c3.push_back(lps[i][j]);
		if(RMSD(c3,cas3)<0.6) return true;
	}
	return false;
}

int main(int argc, char ** argv) {
	NSPdstl::RandomEngine<>::getinstance().reseed(5918);
	std::vector<std::vector<Residue>> chs=ReadPDB(std::string(argv[1]), true);
	std::vector<std::vector<int>> sten=Findsten(std::string(argv[2]));
	std::vector<std::vector<Residue>> lps=ReadHTHLib(std::string(argv[3]),std::string(argv[4]));
	std::vector<std::vector<XYZ>> lpcrd=ExtractCrd(lps);
	std::vector<std::vector<int>> rts;
	for(int i=0;i<sten.size();i++) {
		if(!Halign(lpcrd, chs[0], sten[i][0], sten[i][1])) {
			//std::cout <<i <<" none" <<std::endl;
			return 1;
		}
		//std::cout <<i <<" true" <<std::endl;
	}
	std::string outpdb=std::string(argv[5]);
	residue2pdb(outpdb,chs);
}

