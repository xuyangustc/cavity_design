/*
 * basic_info.cpp
 *
 *  Created on: Nov 18, 2018
 *      Author: xuyang
 */
#include "allatom/basic_info.h"
#include "dataio/datapaths.h"
#include "geometry/xyz.h"
#include "dstl/randomengine.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>
using namespace NSPallatom;

std::map<std::string,char> AA_Atom::aa31() {
	static std::map<std::string,char> aa;
#pragma omp critical
	{
		if(aa.empty()) {
			std::map<std::string,char> aa1 {
				{"GLY",'G'},
				{"ALA",'A'},
				{"VAL",'V'},
				{"LEU",'L'},
				{"ILE",'I'},
				{"SER",'S'},
				{"THR",'T'},
				{"CYS",'C'},
				{"MET",'M'},
				{"ASP",'D'},
				{"GLU",'E'},
				{"ASN",'N'},
				{"GLN",'Q'},
				{"LYS",'K'},
				{"ARG",'R'},
				{"HIS",'H'},
				{"PRO",'P'},
				{"PHE",'F'},
				{"TYR",'Y'},
				{"TRP",'W'}
			};
			aa=aa1;
		}
	}
	return aa;
}

std::map<char,std::string> AA_Atom::aa13() {
	static std::map<char,std::string> aa;
#pragma omp critical
	{
		if(aa.empty()) {
			std::map<char,std::string> aa1 {
				{'G',"GLY"},
				{'A',"ALA"},
				{'V',"VAL"},
				{'L',"LEU"},
				{'I',"ILE"},
				{'S',"SER"},
				{'T',"THR"},
				{'C',"CYS"},
				{'M',"MET"},
				{'D',"ASP"},
				{'E',"GLU"},
				{'N',"ASN"},
				{'Q',"GLN"},
				{'K',"LYS"},
				{'R',"ARG"},
				{'H',"HIS"},
				{'P',"PRO"},
				{'F',"PHE"},
				{'Y',"TYR"},
				{'W',"TRP"}
			};
			aa=aa1;
		}
	}
	return aa;
}

std::set<std::string> AA_Atom::mainchain() {
	static std::set<std::string> a{"N","CA","C","O"};
	return a;
}

std::map<char,std::set<std::string>> AA_Atom::sidechain() {
	static std::map<char,std::set<std::string>> sc;
	std::map<std::string,char> a31;
	if(sc.empty()) {
		a31=aa31();
	}
#pragma omp critical
	{
		if(sc.empty()) {
			std::map<std::string,std::vector<std::vector<std::string>>> sc1{
				{"GLY",{  }},
				{"ALA",{ {"C","CA","CB"} }},
				{"SER",{ {"C","CA","CB","OG"} }},
				{"CYS",{ {"C","CA","CB","SG"} }},
				{"MET",{ {"C","CA","CB","CG","SD","CE"} }},
				{"LYS",{ {"C","CA","CB","CG","CD","CE","NZ"} }},

				{"VAL",{ {"C","CA","CB","CG1"},{"CA","CB","CG2"} }},
				{"THR",{ {"C","CA","CB","OG1"},{"CA","CB","CG2"} }},
				{"LEU",{ {"C","CA","CB","CG","CD1"},{"CB","CG","CD2"} }},
				{"ASP",{ {"C","CA","CB","CG","OD1"},{"CB","CG","OD2"} }},
				{"ASN",{ {"C","CA","CB","CG","OD1"},{"CB","CG","ND2"} }},
				{"GLU",{ {"C","CA","CB","CG","CD","OE1"},{"CG","CD","OE2"} }},
				{"GLN",{ {"C","CA","CB","CG","CD","OE1"},{"CG","CD","NE2"} }},
				{"ARG",{ {"C","CA","CB","CG","CD","NE","CZ","NH1"},{"NE","CZ","NH2"} }},
				{"ILE",{ {"C","CA","CB","CG1","CD1"},{"CA","CB","CG2"} }},

				{"PRO",{ {"C","CA","CB","CG","CD"} }},
				{"HIS",{ {"C","CA","CB","CG","ND1","CE1","NE2","CD2"} }},
				{"PHE",{ {"C","CA","CB","CG","CD1","CE1","CZ"},{"CB","CG","CD2","CE2"} }},
				{"TYR",{ {"C","CA","CB","CG","CD1","CE1","CZ","OH"},{"CB","CG","CD2","CE2"} }},
				{"TRP",{ {"C","CA","CB","CG","CD1","NE1","CE2","CZ2","CH2","CZ3","CE3","CD2"} }},
			};
			//std::map<std::string,char> a31=aa31();
			for(auto &s:sc1) {
				char c=a31.at(s.first);
				std::set<std::string> s1;
				for(auto &s2:s.second) {
					for(auto &s3:s2) {
						if(s3=="C") continue;
						if(s3=="CA") continue;
						s1.insert(s3);
					}
				}
				sc.insert({c,s1});
			}
		}
	}
	return sc;
}

std::map<std::string,std::set<std::string>> AA_Atom::sidechain3() {
	static std::map<std::string,std::set<std::string>> scs;
	std::map<char,std::set<std::string>> sc1;
	std::map<char,std::string> a1 = aa13();
	if(scs.empty()) {
		sc1 = sidechain();
		a1 = aa13();
	}
#pragma omp critical
	{
		if(scs.empty()) {
			//std::map<char,std::set<std::string>> sc1 = sidechain();
			//std::map<char,std::string> a1 = aa13();
			for(auto &s:sc1) {
				scs.insert({a1.at(s.first),s.second});
			}
		}
	}
	return scs;
}

std::map<char,std::map<std::string,double>> AA_Atom::radius() {
	static std::map<char,std::map<std::string,double>> rs;
	std::map<std::string,char> a31;
	if(rs.empty()) {
		a31=aa31();
	}
#pragma omp critical
	{
		if(rs.empty()) {
			//std::cout <<11 <<std::endl;
			//std::map<std::string,char> a31=aa31();
			//std::cout <<2 <<std::endl;
			//std::set<std::string> mc=mainchain();
			//std::map<char,std::set<std::string>> as=sidechain();
			//for(auto &a:as) {
			//	for(auto &s:mc) a.second.insert(s);
			//}
			std::string fn=NSPdataio::datafilename("resLib.dat");
			std::ifstream ifs(fn);
			std::string line;
			std::stringstream ss;
			while(std::getline(ifs,line)) {
				std::string s1,s2;
				int nline;
				ss << line;
				ss >>s1 >>s2 >>nline;
				ss.clear();
				char aaname=a31.at(s2);
				std::map<std::string,double> scrs;
				for(int i=0;i<nline;i++) {
					std::getline(ifs,line);
					std::string s3, s4;
					double d1;
					char c1;
					ss <<line;
					ss >>s3 >>s4 >>d1 >>c1 >>c1 >>c1;
					if(c1=='a') {
						ss >>s3 >>s3;
					}
					ss.clear();
					scrs.insert({s4,d1});
				}
				rs.insert({aaname,scrs});
			}
			ifs.close();
		}
	}
	return rs;
}

std::map<std::string,std::map<std::string,double>> AA_Atom::radius3() {
	static std::map<std::string,std::map<std::string,double>> rs;
	std::map<char,std::map<std::string,double>> rs1;
	std::map<char,std::string> a1;
	if(rs.empty()) {
		rs1 = radius();
		a1 = aa13();
	}
#pragma omp critical
	{
		if(rs.empty()) {
			//std::map<char,std::map<std::string,double>> rs1 = radius();
//#pragma omp critical
	//{
			//std::map<char,std::string> a1 = aa13();
	//}
			for(auto &r:rs1) {
				rs.insert({a1.at(r.first),r.second});
			}
		}
	}
	return rs;
}

std::set<std::pair<char,std::string>> AA_Atom::donors() {
	static std::set<std::pair<char,std::string>> as;
	std::map<char,std::set<std::string>> sc;
	if(as.empty()) {
		std::map<char,std::set<std::string>> sc=sidechain();
	}
#pragma omp critical
	{
		if(as.empty()) {
			//std::map<char,std::set<std::string>> sc=sidechain();
			for(auto &sc1:sc) {
				if(sc1.first!='P') as.insert({sc1.first,"N"});
				for(auto &s:sc1.second) if(s[0]=='N'||s[0]=='O') as.insert({sc1.first,s});
			}
		}
	}
	return as;
}

std::set<std::pair<char,std::string>> AA_Atom::accepters() {
	static std::set<std::pair<char,std::string>> as;
	std::map<char,std::set<std::string>> sc;
	if(as.empty()) {
		std::map<char,std::set<std::string>> sc=sidechain();
	}
#pragma omp critical
	{
		if(as.empty()) {
			std::map<char,std::set<std::string>> sc=sidechain();
			for(auto &sc1:sc) {
				as.insert({sc1.first,"O"});
				for(auto &s:sc1.second) if(s[0]=='N'||s[0]=='O') as.insert({sc1.first,s});
			}
		}
	}
	return as;
}

std::map<std::pair<char,char>,int> AA_Atom::BLOSUM62() {
	static std::map<std::pair<char,char>,int> mt;
#pragma omp critical
	{
		if(mt.empty()) {
			std::map<char,std::vector<int>> dt;
			std::string fn=NSPdataio::datafilename("blosum.mat");
			std::ifstream ifs(fn);
			std::string line;
			std::stringstream ss;
			for(int i=0;i<20;i++) {
				std::getline(ifs,line);
				ss << line;
				char c;
				ss >> c;
				std::vector<int> vi(i+1);
				for(int j=0;j<vi.size();j++) ss >> vi[j];
				ss.clear();
				dt.insert({c,vi});
			}
			std::getline(ifs,line);
			ss << line;
			std::vector<char> ress(20);
			for(int i=0;i<20;i++) ss >> ress[i];
			for(auto &d:dt) {
				std::vector<int> & vi=d.second;
				for(int i=0;i<vi.size();i++) {
					mt[{d.first,ress[i]}] = vi[i];
				}
			}
			ifs.close();
		}
	}
	return mt;
}

std::vector<NSPgeometry::XYZ> AA_Atom::sphere256() {
	static std::vector<NSPgeometry::XYZ> crds;
#pragma omp critical
	{
		if(crds.empty()) {
			std::string fn=NSPdataio::datafilename("atomSasaPoints");
			std::stringstream ss;
			std::string readline;
			std::ifstream ifs(fn);
			while(std::getline(ifs,readline)) {
				double d0, d1, d2;
				ss << readline;
				ss >>d0 >>d1 >>d2;
				ss.clear();
				crds.push_back(NSPgeometry::XYZ(d0,d1,d2));
			}
			ifs.close();
		}
	}
	return crds;
}

std::set<std::string> AA_Atom::polarresidues() {
	static std::set<std::string> a{"SER","THR","ASP","ASN","GLU","GLN","LYS","ARG","HIS","TYR","TRP"};
	return a;
}

std::map<std::string,std::map<std::string,std::set<std::vector<std::string>>>> AA_Atom::sidechainHBunits() {
	static std::map<std::string,std::map<std::string,std::set<std::vector<std::string>>>> scus;
#pragma omp critical
	{
		if(!scus.empty()) {
			std::string fn = NSPdataio::datafilename("SideChainHBUnit");
			std::string readline;
			std::stringstream ss;
			std::ifstream ifs(fn);
			while(std::getline(ifs,readline)) {
				if(readline.empty()) continue;
				if(readline[0]=='#') continue;
				std::vector<std::string> vs(4);
				ss << readline;
				ss >> vs[0] >> vs[1] >> vs[2] >> vs[3];
				ss.clear();
				if(scus.find(vs[0])==scus.end()) scus.insert({vs[0],std::map<std::string,std::set<std::vector<std::string>>>()});
				if(scus.at(vs[0]).find(vs[1])==scus.at(vs[0]).end()) scus.at(vs[0]).insert({vs[1],std::set<std::vector<std::string>>()});
				scus.at(vs[0]).at(vs[1]).insert({vs[2],vs[3]});
			}
			ifs.close();
		}
	}
	return scus;
}

std::map<std::string,std::vector<std::vector<std::string>>> AA_Atom::groups(int index) {
	int ninfo=3;
	static std::vector<std::map<std::string,std::vector<std::vector<std::string>>>> gps;
	//std::cout <<gps.size() <<std::endl;
#pragma omp critical
	{
		if(gps.empty()) {
			gps.resize(ninfo);
			std::string fn = NSPdataio::datafilename("aa2group");
			std::string readline;
			std::stringstream ss;
			std::ifstream ifs(fn);
			while(std::getline(ifs,readline)) {
				if(readline.empty()) continue;
				if(readline[0]=='#') continue;
				std::string rn;
				int ng;
				//std::cout <<readline <<std::endl;
				ss <<readline;
				ss >> rn >> ng;
				ss.clear();
				for(int i=0;i<ninfo;i++) {
					gps[i].insert({rn,std::vector<std::vector<std::string>>()});
				}
				for(int i=0;i<ng;i++) {
					int nl;
					for(int j=0;j<ninfo;j++) {
						std::vector<std::string> ls;
						std::getline(ifs,readline);
						ss <<readline;
						ss >> nl;
						for(int k=0;k<nl;k++) {
							std::string s;
							ss >>s;
							ls.push_back(s);
						}
						ss.clear();
						gps[j].at(rn).push_back(ls);
					}
				}
			}
			ifs.close();
		}
	}
	//std::cout <<gps.size() <<' ' <<index <<std::endl;
	return gps[index];
}

std::map<std::string,std::vector<std::set<std::string>>> AA_Atom::groups() {
	static std::map<std::string,std::vector<std::set<std::string>>> gps;
#pragma omp critical
	{
		if(!gps.empty()) {
			std::map<std::string,std::vector<std::vector<std::string>>> dts=groups(0);
			for(auto &dt:dts) {
				gps.insert({dt.first,std::vector<std::set<std::string>>()});
				for(auto &v:dt.second) {
					gps.at(dt.first).push_back(std::set<std::string>());
					for(auto &s:v) gps.at(dt.first).back().insert(s);
				}
			}
		}
	}
	return gps;
}

std::map<std::string,std::vector<std::vector<std::string>>> AA_Atom::groups_v() {
	return groups(0);
}

std::map<std::string,std::vector<std::set<std::string>>> AA_Atom::groups_center() {
	static std::map<std::string,std::vector<std::set<std::string>>> gps;
#pragma omp critical
	{
		if(!gps.empty()) {
			std::map<std::string,std::vector<std::vector<std::string>>> dts=groups(1);
			for(auto &dt:dts) {
				gps.insert({dt.first,std::vector<std::set<std::string>>()});
				for(auto &v:dt.second) {
					gps.at(dt.first).push_back(std::set<std::string>());
					for(auto &s:v) gps.at(dt.first).back().insert(s);
				}
			}
		}
	}
	return gps;
}

std::map<std::string,std::vector<std::vector<std::string>>> AA_Atom::groups_crd() {
	//return groups(2);
	return groups(1);
}

std::map<std::string,std::vector<std::vector<std::string>>> AA_Atom::groups_force() {
	return groups(3);
}

/*
std::map<std::string,std::vector<std::set<std::string>>> AA_Atom::groups() {
	static std::map<std::string,std::vector<std::set<std::string>>> gps;
	if(!gps.empty()) return gps;
	std::string fn = NSPdataio::datafilename("aa2group");
	std::string readline;
	std::stringstream ss;
	std::ifstream ifs(fn);
	while(std::getline(ifs,readline)) {
		if(readline.empty()) continue;
		if(readline[0]=='#') continue;
		std::string rn;
		int ng;
		//std::cout <<readline <<std::endl;
		ss <<readline;
		ss >> rn >> ng;
		ss.clear();
		gps.insert({rn,std::vector<std::set<std::string>>()});
		for(int i=0;i<ng;i++) {
			std::set<std::string> ls1, ls2;
			std::vector<std::string> ls3;
			int nl;
			std::getline(ifs,readline);
			//std::cout <<'\t' <<readline <<std::endl;
			ss <<readline;
			ss >>nl;
			for(int j=0;j<nl;j++) {
				std::string s;
				ss >> s;
				ls1.insert(s);
			}
			ss.clear();
			std::getline(ifs,readline);
			//std::cout <<'\t' <<readline <<std::endl;
			ss <<readline;
			ss >>nl;
			for(int j=0;j<nl;j++) {
				std::string s;
				ss >> s;
				ls2.insert(s);
			}
			ss.clear();
			std::getline(ifs,readline);
			//std::cout <<'\t' <<readline <<std::endl;
			ss <<readline;
			ss >>nl;
			for(int j=0;j<nl;j++) {
				std::string s;
				ss >> s;
				ls3.push_back(s);
			}
			ss.clear();
			gps.at(rn).push_back(ls1);
		}
	}
	ifs.close();
	return gps;
}

std::map<std::string,std::vector<std::set<std::string>>> AA_Atom::groups_center() {
	static std::map<std::string,std::vector<std::set<std::string>>> gps;
	if(!gps.empty()) return gps;
	std::string fn = NSPdataio::datafilename("aa2group");
	std::string readline;
	std::stringstream ss;
	std::ifstream ifs(fn);
	while(std::getline(ifs,readline)) {
		if(readline.empty()) continue;
		if(readline[0]=='#') continue;
		std::string rn;
		int ng;
		//std::cout <<readline <<std::endl;
		ss <<readline;
		ss >> rn >> ng;
		ss.clear();
		gps.insert({rn,std::vector<std::set<std::string>>()});
		for(int i=0;i<ng;i++) {
			std::set<std::string> ls1, ls2;
			std::vector<std::string> ls3;
			int nl;
			std::getline(ifs,readline);
			//std::cout <<'\t' <<readline <<std::endl;
			ss <<readline;
			ss >>nl;
			for(int j=0;j<nl;j++) {
				std::string s;
				ss >> s;
				ls1.insert(s);
			}
			ss.clear();
			std::getline(ifs,readline);
			//std::cout <<'\t' <<readline <<std::endl;
			ss <<readline;
			ss >>nl;
			for(int j=0;j<nl;j++) {
				std::string s;
				ss >> s;
				ls2.insert(s);
			}
			ss.clear();
			std::getline(ifs,readline);
			//std::cout <<'\t' <<readline <<std::endl;
			ss <<readline;
			ss >>nl;
			for(int j=0;j<nl;j++) {
				std::string s;
				ss >> s;
				ls3.push_back(s);
			}
			ss.clear();
			gps.at(rn).push_back(ls2);
		}
	}
	ifs.close();
	return gps;
}

std::map<std::string,std::vector<std::vector<std::string>>> AA_Atom::groups_crd() {
	static std::map<std::string,std::vector<std::vector<std::string>>> gps;
	if(!gps.empty()) return gps;
	std::string fn = NSPdataio::datafilename("aa2group");
	std::string readline;
	std::stringstream ss;
	std::ifstream ifs(fn);
	while(std::getline(ifs,readline)) {
		if(readline.empty()) continue;
		if(readline[0]=='#') continue;
		std::string rn;
		int ng;
		//std::cout <<readline <<std::endl;
		ss <<readline;
		ss >> rn >> ng;
		ss.clear();
		gps.insert({rn,std::vector<std::vector<std::string>>()});
		for(int i=0;i<ng;i++) {
			std::set<std::string> ls1, ls2;
			std::vector<std::string> ls3;
			int nl;
			std::getline(ifs,readline);
			//std::cout <<'\t' <<readline <<std::endl;
			ss <<readline;
			ss >>nl;
			for(int j=0;j<nl;j++) {
				std::string s;
				ss >> s;
				ls1.insert(s);
			}
			ss.clear();
			std::getline(ifs,readline);
			//std::cout <<'\t' <<readline <<std::endl;
			ss <<readline;
			ss >>nl;
			for(int j=0;j<nl;j++) {
				std::string s;
				ss >> s;
				ls2.insert(s);
			}
			ss.clear();
			std::getline(ifs,readline);
			//std::cout <<'\t' <<readline <<std::endl;
			ss <<readline;
			ss >>nl;
			for(int j=0;j<nl;j++) {
				std::string s;
				ss >> s;
				ls3.push_back(s);
			}
			ss.clear();
			gps.at(rn).push_back(ls3);
		}
	}
	ifs.close();
	return gps;
}
*/


std::map<std::vector<std::vector<std::string>>,double> AA_Atom::atompairclashdis() {
	static std::map<std::vector<std::vector<std::string>>,double> apc;
#pragma omp critical
	{
		if(apc.empty()) {
			std::string fn = NSPdataio::datafilename("AtomPairClash");
			std::ifstream ifs(fn);
			std::string readline;
			std::stringstream ss;
			while(std::getline(ifs,readline)) {
				if(readline.empty()) continue;
				if(readline[0]=='#') continue;
				std::string s0, s1, s2, s3;
				double dis;
				ss << readline;
				ss >> s0 >> s1 >> s2 >> s3 >> dis;
				ss.clear();
				std::vector<std::string> s1s{s1}, s3s{s3};
				if(s0=="VAL" && s1=="CG") {
					s1s.clear();
					s1s.push_back("CG1");
					s1s.push_back("CG2");
				} else if(s0=="LEU" && s1=="CD") {
					s1s.clear();
					s1s.push_back("CD1");
					s1s.push_back("CD2");
				} else if(s0=="ASP" && s1=="OD") {
					s1s.clear();
					s1s.push_back("OD1");
					s1s.push_back("OD2");
				} else if(s0=="GLU" && s1=="OE") {
					s1s.clear();
					s1s.push_back("OE1");
					s1s.push_back("OE2");
				} else if(s0=="ARG" && s1=="NH") {
					s1s.clear();
					s1s.push_back("NH1");
					s1s.push_back("NH2");
				} else if(s0=="PHE") {
					if(s1=="CD") {
						s1s.clear();
						s1s.push_back("CD1");
						s1s.push_back("CD2");
					}
					if(s1=="CE") {
						s1s.clear();
						s1s.push_back("CE1");
						s1s.push_back("CE2");
					}
				} else if(s0=="TYR") {
					if(s1=="CD") {
						s1s.clear();
						s1s.push_back("CD1");
						s1s.push_back("CD2");
					}
					if(s1=="CE") {
						s1s.clear();
						s1s.push_back("CE1");
						s1s.push_back("CE2");
					}
				}
				if(s2=="VAL" && s3=="CG") {
					s3s.clear();
					s3s.push_back("CG1");
					s3s.push_back("CG2");
				} else if(s2=="LEU" && s3=="CD") {
					s3s.clear();
					s3s.push_back("CD1");
					s3s.push_back("CD2");
				} else if(s2=="ASP" && s3=="OD") {
					s3s.clear();
					s3s.push_back("OD1");
					s3s.push_back("OD2");
				} else if(s2=="GLU" && s3=="OE") {
					s3s.clear();
					s3s.push_back("OE1");
					s3s.push_back("OE2");
				} else if(s2=="ARG" && s3=="NH") {
					s3s.clear();
					s3s.push_back("NH1");
					s3s.push_back("NH2");
				} else if(s2=="PHE") {
					if(s3=="CD") {
						s3s.clear();
						s3s.push_back("CD1");
						s3s.push_back("CD2");
					}
					if(s3=="CE") {
						s3s.clear();
						s3s.push_back("CE1");
						s3s.push_back("CE2");
					}
				} else if(s2=="TYR") {
					if(s3=="CD") {
						s3s.clear();
						s3s.push_back("CD1");
						s3s.push_back("CD2");
					}
					if(s3=="CE") {
						s3s.clear();
						s3s.push_back("CE1");
						s3s.push_back("CE2");
					}
				}
				for(int i=0;i<s1s.size();i++) {
					for(int j=0;j<s3s.size();j++) {
						std::vector<std::string> v1{s0,s1s[i]}, v2{s2,s3s[j]};
						apc.insert({{v1,v2},dis});
						apc.insert({{v2,v1},dis});
					}
				}
			}
			ifs.close();
		}
	}
	return apc;
}

std::map<std::vector<std::string>,double> AA_Atom::Weight_SAI() {
	static std::map<std::vector<std::string>,double> ws;
#pragma omp critical
	{
		if(ws.empty()) {
			std::string fn = NSPdataio::datafilename("filter/AtomAverageSAI");
			std::ifstream ifs(fn);
			std::string readline;
			std::stringstream ss;
			double low, high;
			std::getline(ifs,readline);
			low=std::stod(readline);
			std::getline(ifs,readline);
			high=std::stod(readline);
			while(std::getline(ifs,readline)) {
				//if(readline.empty()) continue;
				//if(readline[0]=='#') continue;
				std::string s;
				double d;
				ss<<readline;
				ss >> s >> d;
				ss.clear();
				std::vector<std::string> vs;
				vs.push_back(s.substr(0,3));
				vs.push_back(s.substr(4));
				if(d<low) d=1.0;
				else if(d>high) d=0.0;
				else if(d<(low+high)/2.0) {
					double a = 2.0 / (high-low) /(high-low);
					double delta = d - low;
					d = 1.0 - a * delta * delta;
				} else {
					double a = 2.0 / (high-low) /(high-low);
					double delta = high - d;
					d = a * delta * delta;
				}
				ws.insert({vs,d});
			}
			ifs.close();
		}
	}
	return ws;
}

std::map<std::string,std::vector<std::vector<std::string>>> AA_Atom::Atom3() {
	static std::map<std::string,std::vector<std::vector<std::string>>> a3s;
#pragma omp critical
	{
		if(a3s.empty()) {
			std::string fn = NSPdataio::datafilename("atom3/Atom3Combination");
			std::string readline;
			std::stringstream ss;
			std::ifstream ifs(fn);
			while(std::getline(ifs,readline)) {
				if(readline.empty()) continue;
				if(readline[0]=='#') continue;
				std::string aa;
				std::vector<std::string> a3(3);
				ss << readline;
				ss >>aa;
				for(int i=0;i<a3.size();i++) {
					ss >>a3[i];
				}
				ss.clear();
				if(a3s.find(aa)==a3s.end()) a3s.insert({aa,std::vector<std::vector<std::string>>()});
				a3s.at(aa).push_back(a3);
			}
			ifs.close();
		}
	}
	return a3s;
}

std::map<std::string,std::map<std::string,std::vector<int>>> AA_Atom::CADis() {
	static std::map<std::string,std::map<std::string,std::vector<int>>> cads;
#pragma omp critical
	{
		if(cads.empty()) {
			std::string fn = NSPdataio::datafilename("atom3/CADistance");
			std::string readline;
			std::stringstream ss;
			std::ifstream ifs(fn);
			std::getline(ifs,readline);
			std::getline(ifs,readline);
			std::getline(ifs,readline);
			int nl=std::stoi(readline);
			for(int i=0;i<nl;i++) {
				std::getline(ifs,readline);
				ss <<readline;
				int na;
				std::string aa;
				ss >>aa >>na;
				ss.clear();
				std::map<std::string,std::vector<int>> cad;
				for(int j=0;j<na;j++) {
					std::getline(ifs,readline);
					ss <<readline;
					int n;
					std::string a;
					ss >>a >>n;
					ss.clear();
					cad.insert({a,{n,n}});
				}
				cads.insert({aa,cad});
			}
			std::getline(ifs,readline);
			for(int i=0;i<3;i++) {
				std::getline(ifs,readline);
				std::string a;
				int n;
				ss <<readline;
				ss >>a >>n;
				ss.clear();
				for(auto &cad:cads) {
					cad.second.insert({a,{n}});
				}
			}
			std::getline(ifs,readline);
			for(int i=0;i<3;i++) {
				std::getline(ifs,readline);
				std::string a;
				int n;
				ss <<readline;
				ss >>a >>n;
				ss.clear();
				for(auto &cad:cads) {
					cad.second.at(a).push_back(n);
				}
			}
			std::getline(ifs,readline);
			for(int i=0;i<2;i++) {
				std::getline(ifs,readline);
				std::string a;
				int n,c;
				ss <<readline;
				ss >>a >>n >>c;
				ss.clear();
				cads.at("PRO").insert({a,{n,c}});
			}
			assert(!std::getline(ifs,readline));
			ifs.close();
		}
	}
	return cads;
}

std::vector<NSPgeometry::XYZ> AA_Atom::PolyHedron20() {
	static std::vector<NSPgeometry::XYZ> cs;
#pragma omp critical
	{
		if(cs.empty()) {
			std::string fn=NSPdataio::datafilename("PolyHedron20");
			std::ifstream ifs(fn);
			std::string line;
			std::stringstream ss;
			while(std::getline(ifs,line)) {
				NSPgeometry::XYZ c;
				ss << line;
				ss >>c.x_ >>c.y_ >>c.z_;
				ss.clear();
				cs.push_back(c);
			}
			ifs.close();
		}
	}
	return cs;
}

std::vector<NSPgeometry::XYZ> AA_Atom::PolyHedron12() {
	static std::vector<NSPgeometry::XYZ> cs;
#pragma omp critical
	{
		if(cs.empty()) {
			std::string fn=NSPdataio::datafilename("PolyHedron12");
			std::ifstream ifs(fn);
			std::string line;
			std::stringstream ss;
			while(std::getline(ifs,line)) {
				NSPgeometry::XYZ c;
				ss << line;
				ss >>c.x_ >>c.y_ >>c.z_;
				ss.clear();
				cs.push_back(c);
			}
			ifs.close();
		}
	}
	return cs;
}


std::vector<std::vector<int>> NSPallatom::GetPermutation_repeat(int n) {
	std::vector<std::vector<int>> vv;
	for(int i=0;i<n;i++) vv.push_back({i});
	for(int i=1;i<n;i++) {
		std::vector<std::vector<int>> vv1=vv;
		vv.clear();
		for(int j=0;j<vv1.size();j++) {
			for(int k=0;k<n;k++) {
				std::vector<int> v=vv1[j];
				v.push_back(k);
				vv.push_back(v);
			}
		}
	}
	return vv;
}

std::vector<std::vector<int>> NSPallatom::GetPermutation_norepeat(int n) {
	std::vector<std::vector<int>> vv = GetPermutation_repeat(n);
	std::vector<std::vector<int>> vs;
	for(int i=0;i<vv.size();i++) {
		std::set<int> si;
		for(int j:vv[i]) si.insert(j);
		if(si.size()==n) vs.push_back(vv[i]);
	}
	return vs;
}

std::vector<std::vector<int>> NSPallatom::ElementCombination(const std::vector<int> &v) {
	for(int i:v) {
		if(i<=0) {
			std::cout <<"Element Combination Wrong!" <<std::endl;
			exit(1);
		}
	}
	std::vector<std::vector<int>> vv;
	for(int i=0;i<v[0];i++) vv.push_back({i});
	for(int i=1;i<v.size();i++) {
		auto vv1 = vv;
		vv.clear();
		for(int j=0;j<vv1.size();j++) {
			for(int k=0;k<v[i];k++) {
				std::vector<int> v1=vv1[j];
				v1.push_back(k);
				vv.push_back(v1);
			}
		}
		//if(vv.empty()) vv = vv1;
	}
	return vv;
}

std::vector<int> NSPallatom::RandomOrder(int n) {
	std::set<int> rcd;
	for(int i=0;i<n;i++) rcd.insert(i);
	std::vector<int> vi;
	auto &rng = NSPdstl::RandomEngine<>::getinstance();
	for(int i=0;i<n;i++) {
		std::vector<int> v;
		for(int j:rcd) v.push_back(j);
		int k;
#pragma omp critical
		{
			k=rng.intrng(0,v.size()-1)();
		}
		vi.push_back(v[k]);
		rcd.erase(v[k]);
	}
	assert(vi.size()==n);
	assert(rcd.size()==0);
	return vi;
}






