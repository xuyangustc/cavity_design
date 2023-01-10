/*
 * filter.cpp
 *
 *  Created on: Jan 14, 2019
 *      Author: xuyang
 */
#include "allatom/filter.h"
#include "dataio/controlfile.h"
#include "sd/genchain.h"
#include <cstdio>
using namespace NSPallatom;
using namespace NSPgeometry;

/*std::vector<std::string> getenergynames() {
	return {"total","bond","ang","impdih","steric","scconf","phipsi","local","s1","s2","scpacking","localbbhb","sitepair"};
}*/
std::vector<std::string> EnergyFilter::eneterm = {"total","bond","ang","impdih",
		"steric","scconf","phipsi","local","s1","s2","scpacking","localbbhb","sitepair"};
std::set<std::string> EnergyFilter::ipdint = {"total","scconf","phipsi","local","s1","localbbhb"};
std::map<std::string,int> EnergyFilter::etseq = {{"total",0},{"bond",1},{"ang",2},{"impdih",3},{"steric",4},
		{"scconf",5},{"phipsi",6},{"local",7},{"s1",8},{"s2",9},{"scpacking",10},{"localbbhb",11},{"sitepair",12}};

StrideRecord::StrideRecord(std::string line) {
	resname = line.substr(5,3);
	chainid = line[9];
	resnumber = std::stoi(line.substr(11,4));
	resorder = std::stoi(line.substr(16,4));
	sscode = line[24];
	for(int i=26;i<=38;i++)
		if(line[i]!=' ')
			ssname.push_back(line[i]);
	phi = std::stod(line.substr(42,7));
	psi = std::stod(line.substr(52,7));
	sasa = std::stod(line.substr(64,5));
}

std::vector<StrideRecord> NSPallatom::PDB2Stride(std::string pdbfile) {
	std::vector<StrideRecord> rds;
	std::string cmd = "stride "+pdbfile;
	FILE *pp = popen(cmd.c_str(),"r");
	if(!pp) {
		std::cout <<"making strfile wrong!" <<std::endl;
		return rds;
	}
	std::vector<std::string> lines;
	char tmp[100];
	while(!feof(pp)) {
		fgets(tmp,100,pp);
		lines.push_back(std::string(tmp));
	}
	pclose(pp);
	for(std::string & line:lines) {
		if(line.substr(0,3)=="ASG")
			rds.push_back(StrideRecord(line));
	}
	return rds;
}

std::map<std::string,std::vector<int>> NSPallatom::getssfromstride(const std::vector<StrideRecord> &srs) {
	std::map<std::string,std::vector<int>> sss{{"loop",std::vector<int>()},{"helix",std::vector<int>()},{"strand",std::vector<int>()}};
	if(srs.empty()) return sss;
	std::vector<std::string> sslb(1);
	char cid=srs[0].chainid;
	for(int i=0;i<srs.size();i++) {
		if(cid!=srs[i].chainid) {
			sslb.push_back(std::string());
			cid=srs[i].chainid;
		}
		char sid=srs[i].sscode;
		sid=(sid=='H'?'H':(sid=='E'?'E':'C'));
		sslb.back().push_back(sid);
	}
	int st,en;
	for(int i=0;i<sslb.size();i++) {
		st=0;
		for(int j=0;j<sslb[i].size();j++) {
			if(sslb[i][j]!=sslb[i][st]) {
				std::string lb=(sslb[i][st]=='H'?"helix":(sslb[i][st]=='E'?"strand":"loop"));
				sss.at(lb).push_back(i);
				sss.at(lb).push_back(st);
				sss.at(lb).push_back(j-1);
				st=j;
			}
			if(j==sslb[i].size()-1) {
				std::string lb=(sslb[i][st]=='H'?"helix":(sslb[i][st]=='E'?"strand":"loop"));
				sss.at(lb).push_back(i);
				sss.at(lb).push_back(st);
				sss.at(lb).push_back(j-1);
				st=j;
			}
		}
	}
	return sss;
}







void EnergyFilter::assign_ene(std::vector<std::vector<std::vector<EAtom>>> &et,
		const std::vector<std::pair<MultiAtom,double>> &es) {
	if(es.size()==0) return;
	et.resize(chains_.size());
	for(int i=0;i<chains_.size();i++) et[i].resize(chains_[i].size());
	int size=es[0].first.size();
	for(const auto &e:es) {
		for(int a=0;a<size;a++) {
			int i=e.first[a].first[0];
			int j=e.first[a].first[1];
			et[i][j].push_back(e);
		}
	}
}

void EnergyFilter::assign_ene(std::vector<std::vector<std::vector<EMultiResidue>>> &et,
		const std::vector<std::pair<MultiResidue,double>> &es) {
	if(es.size()==0) return;
	et.resize(chains_.size());
	for(int i=0;i<chains_.size();i++) et[i].resize(chains_[i].size());
	int pairsize=2;
	assert(pairsize==es[0].first.size());
	for(const auto &e:es) {
		const std::vector<std::vector<int>> & lb = e.first;
		int i0=lb[0][0];
		int j0=lb[0][1];
		int i1=lb[1][0];
		int j1=lb[1][1];
		et[i0][j0].push_back({e});
		et[i1][j1].push_back({e});
	}
}

void EnergyFilter::assign_ene(std::vector<std::vector<double>> &et,
		const std::vector<std::pair<SingleResidue,double>> &es) {
	et.resize(chains_.size());
	for(int i=0;i<chains_.size();i++) et[i].assign(chains_[i].size(),BigValue);
	for(const auto &e:es) {
		int i0=e.first[0];
		int j0=e.first[1];
		et[i0][j0] = e.second;
	}
}

void EnergyFilter::sitetotalene() {
	etot_.resize(chains_.size());
	for(int i=0;i<chains_.size();i++) {
		etot_[i].assign(chains_[i].size(),0.0);
		for(int j=0;j<chains_[i].size();j++) {
			if(!ebond_.empty())
			for(int k=0;k<ebond_[i][j].size();k++) {
				etot_[i][j] += ebond_[i][j][k].second / 2.0;
			}
			if(!eang_.empty())
			for(int k=0;k<eang_[i][j].size();k++) {
				etot_[i][j] += eang_[i][j][k].second / 3.0;
			}
			if(!eimpdih_.empty())
			for(int k=0;k<eimpdih_[i][j].size();k++) {
				etot_[i][j] += eimpdih_[i][j][k].second / 4.0;
			}
			if(!esteric_.empty())
			for(int k=0;k<esteric_[i][j].size();k++) {
				etot_[i][j] += esteric_[i][j][k].second / 2.0;
			}
			if(!elocalbbhb_.empty())
			for(int k=0;k<elocalbbhb_[i][j].size();k++) {
				etot_[i][j] += elocalbbhb_[i][j][k].second / 2.0;
			}
			if(!esitepair_.empty())
			for(int k=0;k<esitepair_[i][j].size();k++) {
				etot_[i][j] += esitepair_[i][j][k].second / 2.0;
			}
			if(!escpacking_.empty())
			for(int k=0;k<escpacking_[i][j].size();k++) {
				etot_[i][j] += escpacking_[i][j][k].second / 2.0;
			}
			if(!es2_.empty())
			for(int k=0;k<es2_[i][j].size();k++) {
				etot_[i][j] += es2_[i][j][k].second / 2.0;
			}
			if(!ephipsi_.empty()) etot_[i][j] += ephipsi_[i][j];
			if(!elocal_.empty()) etot_[i][j] += elocal_[i][j];
			if(!escconf_.empty()) etot_[i][j] += escconf_[i][j];
			if(!es1_.empty()) etot_[i][j] += es1_[i][j];
		}
	}
}

void EnergyFilter::mcene() {
	std::string defaultffparfile=NSPdataio::datafilename("DefaultFF.par");
	NSPdataio::ControlFile cf;
	cf.readfile(defaultffparfile);
	std::vector<std::string> ffcontrolines=cf.getcontrolines("ForceField");
	std::string ffcontrolname{"control_ff"};
	NSPsd::defineforcefieldcontrol(ffcontrolname,ffcontrolines);
	std::vector<std::vector<NSPproteinrep::BackBoneSite>> bsss(chains_.size());
	for(int i=0;i<chains_.size();i++) {
		for(int j=0;j<chains_[i].size();j++) {
			bsss[i].push_back(chains_[i][j].getbackbonesite());
		}
	}
	NSPsd::ForceField ff=NSPsd::make_forcefield_allatom(bsss, ffcontrolname);
	std::string allfixed, mcfixed;
	NSPsd::ActiveSelections acts(&ff,allfixed,mcfixed);
	std::vector<std::vector<NSPproteinrep::FullSite>> fsss(chains_.size());
	for(int i=0;i<chains_.size();i++) {
		for(int j=0;j<chains_[i].size();j++) {
			fsss[i].push_back(chains_[i][j].getfullsite());
		}
	}
	std::vector<std::pair<std::vector<int>,std::string>> atomseq;
	std::vector<double> cs=NSPsd::GenChain::getcrd(fsss,atomseq);
	for(double & c:cs) c=c*A2NM;
	NSPsd::NeighborList nbl(cs,ff);
	/*std::vector<std::pair<std::vector<std::pair<std::vector<int>,std::string>>,double>> esubbond;
	std::vector<std::pair<std::vector<std::pair<std::vector<int>,std::string>>,double>> esubang;
	std::vector<std::pair<std::vector<std::pair<std::vector<int>,std::string>>,double>> esubimpdih;
	std::vector<std::pair<std::vector<std::pair<std::vector<int>,std::string>>,double>> esubsteric;
	std::vector<std::pair<std::vector<int>,double>> esubphipsi;
	std::vector<std::pair<std::vector<int>,double>> esublocal;
	std::vector<std::pair<std::vector<int>,double>> esubscconf;
	std::vector<std::pair<std::vector<std::vector<int>>,double>> esubsitepair;
	std::vector<std::pair<std::vector<std::vector<int>>,double>> esubscpacking;
	std::vector<std::pair<std::vector<std::vector<int>>,double>> esublocalbbhb;*/
	std::vector<double> pots;
	std::vector<int> sizes;
	for(int i=0;i<chains_.size();i++) sizes.push_back(chains_[i].size());
	NSPsd::EnergyTerm etm;
	ff.forces_subentry(cs, nbl, &pots, acts, atomseq, etm, sizes);
	std::vector<std::pair<std::vector<std::pair<std::vector<int>,std::string>>,double>> &esubbond = etm.ebond;
	std::vector<std::pair<std::vector<std::pair<std::vector<int>,std::string>>,double>> &esubang = etm.eang;
	std::vector<std::pair<std::vector<std::pair<std::vector<int>,std::string>>,double>> &esubimpdih = etm.eimpdih;
	std::vector<std::pair<std::vector<std::pair<std::vector<int>,std::string>>,double>> &esubsteric = etm.esteric;
	std::vector<std::pair<std::vector<int>,double>> &esubphipsi = etm.ephipsi;
	std::vector<std::pair<std::vector<int>,double>> &esublocal = etm.elocal;
	std::vector<std::pair<std::vector<int>,double>> &esubscconf = etm.escconf;
	std::vector<std::pair<std::vector<std::vector<int>>,double>> &esubsitepair = etm.esitepair;
	std::vector<std::pair<std::vector<std::vector<int>>,double>> &esubscpacking = etm.escpacking;
	std::vector<std::pair<std::vector<std::vector<int>>,double>> &esublocalbbhb = etm.elocalbbhb;
	assign_ene(ebond_,esubbond);
	assign_ene(eang_,esubang);
	assign_ene(eimpdih_,esubimpdih);
	assign_ene(esteric_,esubsteric);
	assign_ene(ephipsi_,esubphipsi);
	assign_ene(elocal_,esublocal);
	assign_ene(escconf_,esubscconf);
	assign_ene(esitepair_,esubsitepair);
	assign_ene(escpacking_,esubscpacking);
	assign_ene(elocalbbhb_,esublocalbbhb);





	//for check
	/*std::vector<double> ps;
	ff.forces(cs,nbl,&ps,acts);
	std::cout <<'\t' <<pots.size() <<"   " <<ps.size() <<std::endl;
	std::cout <<"Bond:\t" <<pots[NSPsd::ForceField::EBOND] <<"   " <<ps[NSPsd::ForceField::EBOND] <<std::endl;
	std::cout <<"Ang:\t" <<pots[NSPsd::ForceField::EANG] <<"   " <<ps[NSPsd::ForceField::EANG] <<std::endl;
	std::cout <<"ImpDih:\t" <<pots[NSPsd::ForceField::EIMPDIH] <<"   " <<ps[NSPsd::ForceField::EIMPDIH] <<std::endl;
	std::cout <<"Steric:\t" <<pots[NSPsd::ForceField::ESTERIC] <<"   " <<ps[NSPsd::ForceField::ESTERIC] <<std::endl;
	std::cout <<"PhiPsi:\t" <<pots[NSPsd::ForceField::EPHIPSI] <<"   " <<ps[NSPsd::ForceField::EPHIPSI] <<std::endl;
	std::cout <<"Local:\t" <<pots[NSPsd::ForceField::ELOCALSTRUCTURE] <<"   " <<ps[NSPsd::ForceField::ELOCALSTRUCTURE] <<std::endl;
	std::cout <<"SitePair:\t" <<pots[NSPsd::ForceField::ESITEPAIRS] <<"   " <<ps[NSPsd::ForceField::ESITEPAIRS] <<std::endl;
	//std::cout <<"SCConf:\t" <<pots[NSPsd::ForceField::ESCCONF] <<"   " <<ps[NSPsd::ForceField::ESCCONF] <<std::endl;
	//std::cout <<"SCPacking:\t" <<pots[NSPsd::ForceField::ESCPACKING] <<"   " <<ps[NSPsd::ForceField::ESCPACKING] <<std::endl;
	std::cout <<"LocalBBHB:\t" <<pots[NSPsd::ForceField::ELOCALHB] <<"   " <<ps[NSPsd::ForceField::ELOCALHB] <<std::endl;
	std::cout <<std::endl;*/
}

void EnergyFilter::scene(std::string abacusfile, bool readss) {
	std::vector<std::string> sslb;
	std::vector<std::vector<double>> s1s;
	std::vector<std::vector<std::vector<EMultiResidue>>> s2s;
	std::ifstream ifs(abacusfile);
	std::string readline;
	std::stringstream ss;
	std::getline(ifs,readline);
	int nch = std::stoi(readline);
	sslb.resize(nch);
	s1s.resize(nch);
	s2s.resize(nch);
	for(int i=0;i<nch;i++) {
		std::getline(ifs,readline);
		int nres = std::stoi(readline);
		sslb[i].resize(nres);
		s1s[i].resize(nres);
		s2s[i].resize(nres);
		for(int j=0;j<nres;j++) {
			std::getline(ifs,readline);
			std::string aaname;
			char sl;
			double s11;
			int s2size;
			std::vector<EMultiResidue> s21;
			ss << readline;
			ss >> aaname >>sl >>s11 >>s2size;
			for(int k=0;k<s2size;k++) {
				int i0, j0;
				double d0;
				ss >> i0 >> j0 >> d0;
				s21.push_back({{{i,j},{i0,j0}},d0});
			}
			ss.clear();
			sslb[i][j] = sl;
			s1s[i][j] = s11;
			s2s[i][j] = s21;
		}
	}
	ifs.close();
	assert(sslb.size()==chains_.size());
	assert(s1s.size()==chains_.size());
	assert(s2s.size()==chains_.size());
	for(int i=0;i<chains_.size();i++) {
		assert(sslb[i].size()==chains_[i].size());
		assert(s1s[i].size()==chains_[i].size());
		assert(s2s[i].size()==chains_[i].size());
	}
	es2_ = s2s;
	es1_ = s1s;
	if(readss) ss_ = sslb;
}

void EnergyFilter::enelist(std::map<SingleResidue,double>&el, const std::vector<std::vector<double>>&es) {
	el.clear();
	for(int i=0;i<es.size();i++) {
		for(int j=0;j<es[i].size();j++) {
			if(fabs(es[i][j]-BigValue)<0.001) continue;
			el.insert({{i,j},es[i][j]});
		}
	}
}

void EnergyFilter::enelist(std::map<MultiResidue,double>&el,
		const std::vector<std::vector<std::vector<EMultiResidue>>>&es) {
	el.clear();
	for(const auto &vve:es) {
		for(const auto &ve:vve) {
			for(const auto &emr:ve) {
				MultiResidue mr=emr.first;
				std::sort(mr.begin(),mr.end(),[](std::vector<int>v1,
						std::vector<int>v2)->bool{return (v1[0]<v2[0]?true:
								(v1[0]>v2[0]?false:(v1[1]<v2[1])));});
				el.insert({mr,emr.second});
			}
		}
	}
}

void EnergyFilter::enelist(std::map<MultiAtom,double>&el,const std::vector<std::vector<std::vector<EAtom>>>&es) {
	el.clear();
	for(const auto &vve:es) {
		for(const auto &ve:vve) {
			for(const auto &emr:ve) {
				MultiAtom ma=emr.first;
				std::sort(ma.begin(),ma.end(),[](SingleAtom sa1, SingleAtom sa2)->bool{
					return (sa1.first[0]<sa2.first[0]?true:(sa1.first[0]>sa2.first[0]?false:(
							sa1.first[1]<sa2.first[1]?true:(sa1.first[1]>sa2.first[1]?false:(
									sa1.second<sa2.second)))));});
				el.insert({ma,emr.second});
			}
		}
	}
}

void EnergyFilter::enelist() {
	enelist(lbond_,ebond_);
	enelist(lang_,eang_);
	enelist(limpdih_,eimpdih_);
	enelist(lsteric_,esteric_);
	enelist(llocal_,elocal_);
	enelist(lphipsi_,ephipsi_);
	enelist(ls1_,es1_);
	enelist(lsitepair_,esitepair_);
	enelist(llocalbbhb_,elocalbbhb_);
	enelist(ls2_,es2_);
	enelist(lscconf_,escconf_);
	enelist(lscpacking_,escpacking_);
	enelist(ltot_,etot_);
}

void EnergyFilter::island_atom(std::vector<IAtom>&il,const std::map<MultiAtom,double>&el,double ec,double ed) {
	il.clear();
	//std::set<IAtom> il1;
	ed *= ed;
	for(const auto &mp:el) {
		if(mp.second<ec) continue;
		const MultiAtom &ma=mp.first;
		//find assemble
		std::vector<NSPgeometry::XYZ> cs;
		for(const SingleAtom &sa:ma) {
			cs.push_back(getcrd(sa));
		}
		std::vector<int> vi;
		for(int i=0;i<il.size();i++) {
			bool in{false};
			for(MultiAtom &ma1:il[i]) {
				for(SingleAtom &sa:ma1) {
					NSPgeometry::XYZ c=getcrd(sa);
					for(NSPgeometry::XYZ &c1:cs) {
						if((c-c1).squarednorm()<ed) {
							in=true;
							break;
						}
					}
					if(in) break;
				}
				if(in) break;
			}
			if(in) vi.push_back(i);
		}
		//merge the common
		if(vi.empty()) {
			il.push_back(IAtom());
			il.back().push_back(ma);
		} else {
			for(int i=1;i<vi.size();i++) {
				for(auto &ma1:il[vi[i]]) {
					il[vi[0]].push_back(ma1);
				}
				il[vi[i]].clear();
			}
			il[vi[0]].push_back(ma);
		}
	}
	std::vector<IAtom> via;
	for(IAtom &ia:il) {
		if(ia.empty()) continue;
		via.push_back(ia);
	}
	il=via;
}

void EnergyFilter::island_singleresidue(std::vector<IResidue>&il,const std::map<SingleResidue,double>&el,double ec,int ed) {
	il.clear();
	for(const auto &mp:el) {
		if(mp.second<ec) continue;
		const SingleResidue &sr = mp.first;
		//find islands
		std::vector<int> vi;
		for(int i=0;i<il.size();i++) {
			bool in{false};
			for(SingleResidue &sr1:il[i]) {
				if(sr[0]!=sr1[0]) continue;
				if(fabs(sr[1]-sr1[1])>ed) continue;
				in=true;
				break;
			}
			if(in) vi.push_back(i);
		}
		//merge islands
		if(vi.empty()) {
			il.push_back(IResidue());
			il.back().push_back(sr);
		} else {
			for(int i=1;i<vi.size();i++) {
				for(SingleResidue &sr1:il[vi[i]]) {
					il[vi[0]].push_back(sr1);
				}
				il[vi[i]].clear();
			}
			il[vi[0]].push_back(sr);
		}
	}
	std::vector<IResidue> il1;
	for(IResidue &ir:il) {
		if(ir.empty()) continue;
		il1.push_back(ir);
	}
	il=il1;
}

void EnergyFilter::island_multiresidue_int(std::vector<IMultiResidue>&il,const std::map<MultiResidue,double>&el,double ec,int ed) {
	il.clear();
	for(const auto &mp:el) {
		if(mp.second<ec) continue;
		const MultiResidue &mr = mp.first;
		//find islands
		std::vector<int> vi;
		for(int i=0;i<il.size();i++) {
			bool in{false};
			for(MultiResidue &mr1:il[i]) {
				for(SingleResidue &sr1:mr1) {
					for(const SingleResidue &sr:mr) {
						if(sr[0]==sr1[0] && fabs(sr[1]-sr1[1])<=ed) {
							in=true;
							break;
						}
						if(in) break;
					}
					if(in) break;
				}
				if(in) break;
			}
			if(in) vi.push_back(i);
		}
		//merge islands
		if(vi.empty()) {
			il.push_back(IMultiResidue());
			il.back().push_back(mr);
		} else {
			for(int i=1;i<vi.size();i++) {
				for(MultiResidue &mr1:il[vi[i]]) {
					il[vi[0]].push_back(mr1);
				}
				il[vi[i]].clear();
			}
			il[vi[0]].push_back(mr);
		}
	}
	std::vector<IMultiResidue> il1;
	for(IMultiResidue &ir:il) {
		if(ir.empty()) continue;
		il1.push_back(ir);
	}
	il=il1;
}

void EnergyFilter::island_multiresidue_double(
		std::vector<IMultiResidue>&il,const std::map<MultiResidue,double>&el,double ec,double ed,bool onlymc) {
	il.clear();
	ed *= ed;
	for(const auto &mp:el) {
		if(mp.second<ec) continue;
		const MultiResidue &mr = mp.first;
		//find islands
		std::vector<NSPgeometry::XYZ> cs;
		for(const SingleResidue &sr:mr) {
			std::vector<NSPgeometry::XYZ> cs1 = getcrds(sr,onlymc);
			for(auto &c:cs1) cs.push_back(c);
		}
		std::vector<int> vi;
		for(int i=0;i<il.size();i++) {
			bool in{false};
			for(MultiResidue &mr1:il[i]) {
				for(SingleResidue &sr1:mr1) {
					std::vector<NSPgeometry::XYZ> cs1=getcrds(sr1,onlymc);
					for(auto &c1:cs1) {
						for(auto &c:cs) {
							if((c-c1).squarednorm()<ed) {
								in = true;
								break;
							}
						}
						if(in) break;
					}
					if(in) break;
				}
				if(in) break;
			}
			if(in) vi.push_back(i);
		}
		//merge islands
		if(vi.empty()) {
			il.push_back(IMultiResidue());
			il.back().push_back(mr);
		} else {
			for(int i=1;i<vi.size();i++) {
				for(MultiResidue &mr1:il[vi[i]]) {
					il[vi[0]].push_back(mr1);
				}
				il[vi[i]].clear();
			}
			il[vi[0]].push_back(mr);
		}
	}
	std::vector<IMultiResidue> il1;
	for(IMultiResidue &ir:il) {
		if(ir.empty()) continue;
		il1.push_back(ir);
	}
	il=il1;
}

void EnergyFilter::island() {
	island_atom(ibond_,lbond_,parene_[ENE::BOND],pardisdou_[ENE::BOND]);
	island_atom(iang_,lang_,parene_[ENE::ANG],pardisdou_[ENE::ANG]);
	island_atom(iimpdih_,limpdih_,parene_[ENE::IMPDIH],pardisdou_[ENE::IMPDIH]);
	island_atom(isteric_,lsteric_,parene_[ENE::STERIC],pardisdou_[ENE::STERIC]);
	island_multiresidue_int(ilocalbbhb_,llocalbbhb_,parene_[ENE::LOCALBBHB],pardisint_[ENE::LOCALBBHB]);
	island_multiresidue_double(isitepair_,lsitepair_,parene_[ENE::SITEPAIR],pardisdou_[ENE::SITEPAIR],true);
	island_multiresidue_double(iscpacking_,lscpacking_,parene_[ENE::SCPACKING],pardisdou_[ENE::SCPACKING],false);
	island_multiresidue_double(is2_,ls2_,parene_[ENE::S2],pardisdou_[ENE::S2],false);
	island_singleresidue(iphipsi_,lphipsi_,parene_[ENE::PHIPSI],pardisint_[ENE::PHIPSI]);
	island_singleresidue(ilocal_,llocal_,parene_[ENE::LOCAL],pardisint_[ENE::LOCAL]);
	island_singleresidue(iscconf_,lscconf_,parene_[ENE::SCCONF],pardisint_[ENE::SCCONF]);
	island_singleresidue(is1_,ls1_,parene_[ENE::S1],pardisint_[ENE::S1]);
	island_singleresidue(itot_,ltot_,parene_[ENE::TOTAL],pardisint_[ENE::TOTAL]);
}

void EnergyFilter::setislandpar() {
	parene_.clear();
	pardisint_.clear();
	pardisdou_.clear();
	maxisland_.clear();
	parene_.assign(ENE::ETERM,10000.0);
	pardisint_.assign(ENE::ETERM,0);
	pardisdou_.assign(ENE::ETERM,0.0);
	maxisland_.assign(ENE::ETERM,100000);
	std::string parfile=NSPdataio::datafilename("IslandCutOffMC95");
	std::ifstream ifs(parfile);
	std::string readline;
	std::stringstream ss;
	// distance, maxnum, energy
	while(std::getline(ifs,readline)) {
		if(readline.empty()) continue;
		std::string lb;
		int di,mi;
		double dd, ed;
		ss << readline;
		ss >> lb;
		if(etseq.find(lb)==etseq.end()) {
			std::cout <<"can not recognize label :  " <<lb <<std::endl;
			exit(1);
		}
		if(ipdint.find(lb)==ipdint.end()) ss >> dd;
		else ss >> di;
		ss >> mi >> ed;
		ss.clear();
		int ix=etseq.at(lb);
		parene_[ix] = ed;
		maxisland_[ix] = mi;
		if(ipdint.find(lb)==ipdint.end()) {
			pardisdou_[ix] = dd;
		} else {
			pardisint_[ix] = di;
		}
	}
	ifs.close();
}

std::vector<std::string> EnergyFilter::enetermfilter() {
	std::vector<std::string> lbs;
	for(const auto &ie:ibond_) {
		if(ie.size()<maxisland_[ENE::BOND]) continue;
		lbs.push_back(eneterm[ENE::BOND]);
		break;
	}
	for(const auto &ie:iang_) {
		if(ie.size()<maxisland_[ENE::ANG]) continue;
		lbs.push_back(eneterm[ENE::ANG]);
		break;
	}
	for(const auto &ie:iimpdih_) {
		if(ie.size()<maxisland_[ENE::IMPDIH]) continue;
		lbs.push_back(eneterm[ENE::IMPDIH]);
		break;
	}
	for(const auto &ie:isteric_) {
		if(ie.size()<maxisland_[ENE::STERIC]) continue;
		lbs.push_back(eneterm[ENE::STERIC]);// std::cout <<ie.size() <<std::endl;
		break;
	}
	for(const auto &ie:iscconf_) {
		if(ie.size()<maxisland_[ENE::SCCONF]) continue;
		lbs.push_back(eneterm[ENE::SCCONF]);
		break;
	}
	for(const auto &ie:iphipsi_) {
		if(ie.size()<maxisland_[ENE::PHIPSI]) continue;
		lbs.push_back(eneterm[ENE::PHIPSI]);
		break;
	}
	for(const auto &ie:ilocal_) {
		if(ie.size()<maxisland_[ENE::LOCAL]) continue;
		lbs.push_back(eneterm[ENE::LOCAL]);
		break;
	}
	for(const auto &ie:is1_) {
		if(ie.size()<maxisland_[ENE::S1]) continue;
		lbs.push_back(eneterm[ENE::S1]);
		break;
	}
	for(const auto &ie:ilocalbbhb_) {
		if(ie.size()<maxisland_[ENE::LOCALBBHB]) continue;
		lbs.push_back(eneterm[ENE::LOCALBBHB]);
		break;
	}
	for(const auto &ie:is2_) {
		if(ie.size()<maxisland_[ENE::S2]) continue;
		lbs.push_back(eneterm[ENE::S2]);
		break;
	}
	for(const auto &ie:isitepair_) {
		if(ie.size()<maxisland_[ENE::SITEPAIR]) continue;
		lbs.push_back(eneterm[ENE::SITEPAIR]);// std::cout <<ie.size() <<std::endl;
		break;
	}
	for(const auto &ie:iscpacking_) {
		if(ie.size()<maxisland_[ENE::SCPACKING]) continue;
		lbs.push_back(eneterm[ENE::SCPACKING]);
		break;
	}
	for(const auto &ie:itot_) {
		if(ie.size()<maxisland_[ENE::TOTAL]) continue;
		lbs.push_back(eneterm[ENE::TOTAL]);
		break;
	}
	return lbs;
}

std::map<std::string,std::pair<int,double>> EnergyFilter::highene() {
	std::map<std::string,std::pair<int,double>> hes;
	hes.insert({eneterm[ENE::BOND],{0,0.0}});
	for(auto &e:lbond_) {
		if(e.second < parene_[ENE::BOND]) continue;
		hes.at(eneterm[ENE::BOND]).first++;
		hes.at(eneterm[ENE::BOND]).second += e.second;
	}
	hes.insert({eneterm[ENE::ANG],{0,0.0}});
	for(auto &e:lang_) {
		if(e.second < parene_[ENE::ANG]) continue;
		hes.at(eneterm[ENE::ANG]).first++;
		hes.at(eneterm[ENE::ANG]).second += e.second;
	}
	hes.insert({eneterm[ENE::IMPDIH],{0,0.0}});
	for(auto &e:limpdih_) {
		if(e.second < parene_[ENE::IMPDIH]) continue;
		hes.at(eneterm[ENE::IMPDIH]).first++;
		hes.at(eneterm[ENE::IMPDIH]).second += e.second;
	}
	hes.insert({eneterm[ENE::STERIC],{0,0.0}});
	for(auto &e:lsteric_) {
		if(e.second < parene_[ENE::STERIC]) continue;
		hes.at(eneterm[ENE::STERIC]).first++;
		hes.at(eneterm[ENE::STERIC]).second += e.second;
		//std::cout <<e.first[0].first[0] <<'\t' <<e.first[0].first[1] <<'\t' <<e.first[0].second <<'\t'
		//		<<e.first[1].first[0] <<'\t' <<e.first[1].first[1] <<'\t' <<e.first[1].second <<'\t' <<e.second <<std::endl;
	}
	hes.insert({eneterm[ENE::PHIPSI],{0,0.0}});
	for(auto &e:lphipsi_) {
		if(e.second < parene_[ENE::PHIPSI]) continue;
		hes.at(eneterm[ENE::PHIPSI]).first++;
		hes.at(eneterm[ENE::PHIPSI]).second += e.second;
	}
	hes.insert({eneterm[ENE::LOCAL],{0,0.0}});
	for(auto &e:llocal_) {
		if(e.second < parene_[ENE::LOCAL]) continue;
		hes.at(eneterm[ENE::LOCAL]).first++;
		hes.at(eneterm[ENE::LOCAL]).second += e.second;
	}
	hes.insert({eneterm[ENE::SCCONF],{0,0.0}});
	for(auto &e:lscconf_) {
		if(e.second < parene_[ENE::SCCONF]) continue;
		hes.at(eneterm[ENE::SCCONF]).first++;
		hes.at(eneterm[ENE::SCCONF]).second += e.second;
	}
	hes.insert({eneterm[ENE::S1],{0,0.0}});
	for(auto &e:ls1_) {
		if(e.second < parene_[ENE::S1]) continue;
		hes.at(eneterm[ENE::S1]).first++;
		hes.at(eneterm[ENE::S1]).second += e.second;
	}
	hes.insert({eneterm[ENE::SITEPAIR],{0,0.0}});
	for(auto &e:lsitepair_) {
		if(e.second < parene_[ENE::SITEPAIR]) continue;
		hes.at(eneterm[ENE::SITEPAIR]).first++;
		hes.at(eneterm[ENE::SITEPAIR]).second += e.second;
	}
	hes.insert({eneterm[ENE::SCPACKING],{0,0.0}});
	for(auto &e:lscpacking_) {
		if(e.second < parene_[ENE::SCPACKING]) continue;
		hes.at(eneterm[ENE::SCPACKING]).first++;
		hes.at(eneterm[ENE::SCPACKING]).second += e.second;
	}
	hes.insert({eneterm[ENE::LOCALBBHB],{0,0.0}});
	for(auto &e:llocalbbhb_) {
		if(e.second < parene_[ENE::LOCALBBHB]) continue;
		hes.at(eneterm[ENE::LOCALBBHB]).first++;
		hes.at(eneterm[ENE::LOCALBBHB]).second += e.second;
	}
	hes.insert({eneterm[ENE::S2],{0,0.0}});
	for(auto &e:ls2_) {
		if(e.second < parene_[ENE::S2]) continue;
		hes.at(eneterm[ENE::S2]).first++;
		hes.at(eneterm[ENE::S2]).second += e.second;
	}
	return hes;
}

std::map<std::string,double> EnergyFilter::totalene() {
	std::map<std::string,std::pair<int,double>> hes;
	hes.insert({eneterm[ENE::BOND],{0,0.0}});
	for(auto &e:lbond_) {
		//if(e.second < parene_[ENE::BOND]) continue;
		hes.at(eneterm[ENE::BOND]).first++;
		hes.at(eneterm[ENE::BOND]).second += e.second;
	}
	hes.insert({eneterm[ENE::ANG],{0,0.0}});
	for(auto &e:lang_) {
		//if(e.second < parene_[ENE::ANG]) continue;
		hes.at(eneterm[ENE::ANG]).first++;
		hes.at(eneterm[ENE::ANG]).second += e.second;
	}
	hes.insert({eneterm[ENE::IMPDIH],{0,0.0}});
	for(auto &e:limpdih_) {
		//if(e.second < parene_[ENE::IMPDIH]) continue;
		hes.at(eneterm[ENE::IMPDIH]).first++;
		hes.at(eneterm[ENE::IMPDIH]).second += e.second;
	}
	hes.insert({eneterm[ENE::STERIC],{0,0.0}});
	for(auto &e:lsteric_) {
		//if(e.second < parene_[ENE::STERIC]) continue;
		hes.at(eneterm[ENE::STERIC]).first++;
		hes.at(eneterm[ENE::STERIC]).second += e.second;
	}
	hes.insert({eneterm[ENE::PHIPSI],{0,0.0}});
	for(auto &e:lphipsi_) {
		//if(e.second < parene_[ENE::PHIPSI]) continue;
		hes.at(eneterm[ENE::PHIPSI]).first++;
		hes.at(eneterm[ENE::PHIPSI]).second += e.second;
	}
	hes.insert({eneterm[ENE::LOCAL],{0,0.0}});
	for(auto &e:llocal_) {
		//if(e.second < parene_[ENE::LOCAL]) continue;
		hes.at(eneterm[ENE::LOCAL]).first++;
		hes.at(eneterm[ENE::LOCAL]).second += e.second;
	}
	hes.insert({eneterm[ENE::SCCONF],{0,0.0}});
	for(auto &e:lscconf_) {
		//if(e.second < parene_[ENE::SCCONF]) continue;
		hes.at(eneterm[ENE::SCCONF]).first++;
		hes.at(eneterm[ENE::SCCONF]).second += e.second;
	}
	hes.insert({eneterm[ENE::S1],{0,0.0}});
	for(auto &e:ls1_) {
		//if(e.second < parene_[ENE::S1]) continue;
		hes.at(eneterm[ENE::S1]).first++;
		hes.at(eneterm[ENE::S1]).second += e.second;
	}
	hes.insert({eneterm[ENE::SITEPAIR],{0,0.0}});
	for(auto &e:lsitepair_) {
		//if(e.second < parene_[ENE::SITEPAIR]) continue;
		hes.at(eneterm[ENE::SITEPAIR]).first++;
		hes.at(eneterm[ENE::SITEPAIR]).second += e.second;
	}
	hes.insert({eneterm[ENE::SCPACKING],{0,0.0}});
	for(auto &e:lscpacking_) {
		//if(e.second < parene_[ENE::SCPACKING]) continue;
		hes.at(eneterm[ENE::SCPACKING]).first++;
		hes.at(eneterm[ENE::SCPACKING]).second += e.second;
	}
	hes.insert({eneterm[ENE::LOCALBBHB],{0,0.0}});
	for(auto &e:llocalbbhb_) {
		//if(e.second < parene_[ENE::LOCALBBHB]) continue;
		hes.at(eneterm[ENE::LOCALBBHB]).first++;
		hes.at(eneterm[ENE::LOCALBBHB]).second += e.second;
	}
	hes.insert({eneterm[ENE::S2],{0,0.0}});
	for(auto &e:ls2_) {
		//if(e.second < parene_[ENE::S2]) continue;
		hes.at(eneterm[ENE::S2]).first++;
		hes.at(eneterm[ENE::S2]).second += e.second;
	}
	std::map<std::string,double> tes;
	for(auto &mp:hes) {
		if(mp.second.first==0) continue;
		tes.insert({mp.first,mp.second.second});
	}
	return tes;
}

std::vector<std::vector<int>> Filter::looplthfilter(std::vector<std::vector<int>>&acts, std::string pdbfile) {
	std::set<std::vector<int>> ats;
	for(int i=0;i<acts.size();i++) {
		for(int j=acts[i][1];j<=acts[i][2];j++) ats.insert({acts[i][0],j});
	}
	std::vector<StrideRecord> p2s=PDB2Stride(pdbfile);
	std::map<std::string,std::vector<int>> sss=getssfromstride(p2s);
	std::vector<std::vector<int>> lps;
	std::vector<int> &vi=sss.at("loop");
	for(int i=0;i<vi.size();i+=3) {
		bool isact{false};
		for(int j=vi[i+1];j<=vi[i+2];j++) {
			if(ats.find({vi[i],j})==ats.end()) continue;
			isact=true;
			break;
		}
		if(!isact) continue;
		lps.push_back({vi[i],vi[i+1],vi[i+2]});
	}
	return lps;
}

void Filter::gethbps(bool diffpro) {
	const auto chs_ = chains_;
	std::set<std::pair<std::string,std::string>> drs, acs;
	std::set<std::pair<char,std::string>> drs1=AA_Atom::donors();
	if(!diffpro) drs1.insert({'P',"N"});
	std::set<std::pair<char,std::string>> acs1=AA_Atom::accepters();
	std::map<char,std::string> aa13=AA_Atom::aa13();
	for(auto &p:drs1) drs.insert({aa13.at(p.first),p.second});
	for(auto &p:acs1) acs.insert({aa13.at(p.first),p.second});
	std::map<std::pair<std::string,std::string>,double> rads;
	std::map<char,std::map<std::string,double>> rads1=AA_Atom::radius();
	for(auto &p1:rads1) {
		for(auto &p2:p1.second) {
			rads.insert({{aa13.at(p1.first),p2.first},p2.second});
		}
	}
	hbps_.resize(chs_.size());
	for(int i=0;i<chs_.size();i++) hbps_[i].resize(chs_[i].size());
	for(int i=0;i<chs_.size();i++) {
		for(int j=0;j<chs_[i].size();j++) {
			hbps_[i][j].chainid=i;
			hbps_[i][j].resid=j;
			hbps_[i][j].resname=chs_[i][j].resname();
			std::map<std::string,AtomDataInPDB> rds=chs_[i][j].rds();
			for(auto &a:rds) {
				if(drs.find({hbps_[i][j].resname,a.first})!=drs.end()) {
					hbps_[i][j].hbs.insert({a.first,std::set<std::pair<std::vector<int>,std::string>>()});
					continue;
				}
				if(acs.find({hbps_[i][j].resname,a.first})!=acs.end()) {
					hbps_[i][j].hbs.insert({a.first,std::set<std::pair<std::vector<int>,std::string>>()});
					continue;
				}
			}
		}
	}
	//check hb in the same residue, one from sidechain, another from mainchain
	std::set<std::string> mcp{"N","O"};
	for(int i0=0;i0<hbps_.size();i0++) {
		for(int j0=0;j0<hbps_[i0].size();j0++) {
			for(auto &a:hbps_[i0][j0].hbs) {
				if(mcp.find(a.first)!=mcp.end()) continue;
				NSPgeometry::XYZ c1=chs_[i0][j0].rds().at(a.first).crd;
				double d1=rads.at({hbps_[i0][j0].resname,a.first});
				for(auto s:mcp) {
					if(hbps_[i0][j0].hbs.find(s)==hbps_[i0][j0].hbs.end()) continue;
					NSPgeometry::XYZ c2=chs_[i0][j0].rds().at(s).crd;
					double d2=rads.at({hbps_[i0][j0].resname,s});
					double d3=d1+d2+rad_add_;
					if((c1-c2).squarednorm()>d3*d3) continue;
					hbps_[i0][j0].hbs.at(a.first).insert({{i0,j0},s});
					hbps_[i0][j0].hbs.at(s).insert({{i0,j0},a.first});
				}
			}
		}
	}
	//check hb in different residue
	for(int i0=0;i0<hbps_.size();i0++) {
		for(int j0=0;j0<hbps_[i0].size();j0++) {
			for(auto &h0:hbps_[i0][j0].hbs) {
				for(int i1=i0;i1<hbps_.size();i1++) {
					for(int j1=0;j1<hbps_[i1].size();j1++) {
						if(i0==i1 && j0>=j1) continue;
						for(auto &h1:hbps_[i1][j1].hbs) {
							if(i0==i1) {
								if(j0+1==j1) {
									if(h0.first=="N" && h1.first=="O") continue;
									if(h0.first=="O" && h1.first=="N") continue;
								}
							}
							if( (drs.find({hbps_[i0][j0].resname,h0.first})!=drs.end() &&
									acs.find({hbps_[i1][j1].resname,h1.first})!=acs.end()) || (
											acs.find({hbps_[i0][j0].resname,h0.first})!=acs.end() &&
													drs.find({hbps_[i1][j1].resname,h1.first})!=drs.end()) ) {
								double d2=(chs_[i0][j0].rds().at(h0.first).crd-chs_[i1][j1].rds().at(h1.first).crd).squarednorm();
								double d1=rads.at({hbps_[i0][j0].resname,h0.first})+rads.at({hbps_[i1][j1].resname,h1.first})+rad_add_;
								if(d2>d1*d1) continue;
								h0.second.insert({{i1,j1},h1.first});
								h1.second.insert({{i0,j0},h0.first});
							}
						}
					}
				}
			}
		}
	}
}

std::vector<int> Filter::hbinhelix(int cid, int st, int en) {
	int nt=0, nhb=0;
	for(int i=st;i<en-4;i++) {
		std::map<std::string,std::set<std::pair<std::vector<int>,std::string>>> &hbs1=hbps_[cid][i].hbs;
		std::map<std::string,std::set<std::pair<std::vector<int>,std::string>>> &hbs2=hbps_[cid][i+4].hbs;
		if(hbs1.find("O")==hbs1.end()) continue;
		if(hbs2.find("N")==hbs2.end()) continue;
		nt++;
		if(hbs1.at("O").find({{cid,i+4},"N"})==hbs1.at("O").end()) continue;
		nhb++;
	}
	return {nt,nhb};
}

std::vector<int> Filter::hbinstrand(int c0, int s0, int e0, int c1, int s1, int e1) {
	bool previssml=(e0-s0<e1-s1?true:false);
	int nt=(previssml?e0-s0:e1-s1);
	int nhb=0;
	int i0,j0,k0,i1,j1,k1;
	if(previssml) {
		i0=c0;
		j0=s0;
		k0=e0;
		i1=c1;
		j1=s1;
		k1=e1;
	} else {
		i0=c1;
		j0=s1;
		k0=e1;
		i1=c0;
		j1=s0;
		k1=e0;
	}
	for(int i=j0;i<k0;i++) {
		std::set<std::pair<std::vector<int>,std::string>> &hbs1=hbps_[i0][i].hbs.at("N");
		bool finded{false};
		for(auto &p:hbs1) {
			if(p.first[0]==i1) {
				if(p.first[1]>=j1 && p.first[1]<k1) finded=true;
			}
		}
		if(finded) nhb++;
		std::set<std::pair<std::vector<int>,std::string>> &hbs2=hbps_[i0][i].hbs.at("O");
		finded=false;
		for(auto &p:hbs2) {
			if(p.first[0]==i1) {
				if(p.first[1]>=j1 && p.first[1]<k1) finded=true;
			}
		}
		if(finded) nhb++;
	}
	return {nt,nhb};
}

void Filter::mainchainHBfilter(const std::vector<std::vector<int>>&helixs, std::vector<std::vector<int>> &hbh,
		const std::vector<std::vector<int>>&strands, std::vector<std::vector<int>> &hbs) {
	for(auto &p:helixs) {
		hbh.push_back(hbinhelix(p[0],p[1],p[2]));
	}
	for(auto &p:strands) {
		hbs.push_back(hbinstrand(p[0],p[1],p[2],p[3],p[4],p[5]));
	}
}

std::vector<std::vector<double>> Filter::mainchainHBfilter(const std::vector<std::vector<int>>&helixs,
		const std::vector<std::vector<int>>&strands) {
	if(hbps_.empty()) gethbps();
	std::vector<std::vector<int>> hbh,hbs;
	mainchainHBfilter(helixs,hbh,strands,hbs);
	std::vector<std::vector<double>> prop(2);
	//std::cout <<"Helix  " <<hbh.size() <<std::endl;
	for(auto &vi:hbh) {
		double d=(double)vi[1]/(double)vi[0];
		prop[0].push_back(d);
		//if(d<0.7) return false;
		//std::cout <<vi[1] <<'\t' <<vi[0] <<'\t' <<d <<std::endl;
	}
	//std::cout <<"Strand " <<hbs.size() <<std::endl;
	for(auto &vi:hbs) {
		double d=(double)vi[1]/(double)vi[0];
		prop[1].push_back(d);
		//if(d<0.2) return false;
		//std::cout <<vi[1] <<'\t' <<vi[0] <<'\t' <<d <<std::endl;
	}
	return prop;
}

DECOY Filter::psudosidechain(const DECOY &chs1) {
	DECOY chs=chs1;
	std::set<std::string> mcas{"C","CA","N","O"};
	for(int i=0;i<chs.size();i++) {
		for(int j=0;j<chs[i].size();j++) {
			XYZ n=chs[i][j].rds().at("N").crd;
			XYZ ca=chs[i][j].rds().at("CA").crd;
			XYZ c=chs[i][j].rds().at("C").crd;
			XYZ can = ~(n-ca);
			XYZ cac = ~(c-ca);
			XYZ z = ~(can^cac);
			XYZ x = ~(can+cac);
			XYZ y = ~(z^x);
			LocalFrame lf;
			lf.origin_ = ca;
			lf.axis_.push_back(x);
			lf.axis_.push_back(y);
			lf.axis_.push_back(z);
			std::vector<XYZ> psdList(3);
			psdList[0] = XYZ(-1.040, 0.006, 1.345); /* CB */
			psdList[1] = XYZ(-2.332, -0.663, 0.981); /* CG */
			psdList[2] = XYZ(-3.372, -0.657, 2.325); /* CD */
			XYZ CB = lf.local2globalcrd(psdList[0]);
			XYZ CG = lf.local2globalcrd(psdList[1]);
			XYZ CD = lf.local2globalcrd(psdList[2]);
			AtomDataInPDB acb(CB,'C');
			AtomDataInPDB acg(CG,'C');
			AtomDataInPDB acd(CD,'C');
			std::map<std::string,AtomDataInPDB> rds;
			auto &rds1 = chs[i][j].rds();
			for(auto &mc:mcas) rds.insert({mc,rds1.at(mc)});
			rds.insert({"CB",acb});
			rds.insert({"CG",acg});
			rds.insert({"CD",acd});
			chs[i][j].rds() = rds;
		}
	}
	return chs;
}

void Filter::mcpolarexpose() {
	std::vector<std::vector<Residue>> chs=chains_;
	std::set<std::string> mcas{"C","CA","N","O"};
	std::vector<std::pair<XYZ,SingleAtom>> allas;
	std::vector<std::pair<XYZ,SingleAtom>> polas;
	for(int i=0;i<chs.size();i++) {
		for(int j=0;j<chs[i].size();j++) {
			XYZ n=chs[i][j].rds().at("N").crd;
			XYZ ca=chs[i][j].rds().at("CA").crd;
			XYZ c=chs[i][j].rds().at("C").crd;
			XYZ can = ~(n-ca);
			XYZ cac = ~(c-ca);
			XYZ z = ~(can^cac);
			XYZ x = ~(can+cac);
			XYZ y = ~(z^x);
			LocalFrame lf;
			lf.origin_ = ca;
			lf.axis_.push_back(x);
			lf.axis_.push_back(y);
			lf.axis_.push_back(z);
			std::vector<XYZ> psdList(3);
			psdList[0] = XYZ(-1.040, 0.006, 1.345); /* CB */
			psdList[1] = XYZ(-2.332, -0.663, 0.981); /* CG */
			psdList[2] = XYZ(-3.372, -0.657, 2.325); /* CD */
			XYZ CB = lf.local2globalcrd(psdList[0]);
			XYZ CG = lf.local2globalcrd(psdList[1]);
			XYZ CD = lf.local2globalcrd(psdList[2]);

			allas.push_back({chs[i][j].rds().at("N").crd,{{i,j},"N"}});
			allas.push_back({chs[i][j].rds().at("C").crd,{{i,j},"C"}});
			allas.push_back({chs[i][j].rds().at("CA").crd,{{i,j},"CA"}});
			allas.push_back({chs[i][j].rds().at("O").crd,{{i,j},"O"}});
			allas.push_back({CB,{{i,j},"CB"}});
			allas.push_back({CG,{{i,j},"CG"}});
			allas.push_back({CD,{{i,j},"CD"}});
			polas.push_back({chs[i][j].rds().at("N").crd,{{i,j},"N"}});
			polas.push_back({chs[i][j].rds().at("O").crd,{{i,j},"O"}});
		}
	}
	double rad = 3.5;
	std::vector<XYZ> ps = AA_Atom::sphere256();
	for(XYZ &p:ps) p = p * rad;
	double r2 = rad * rad;
	for(int i=0;i<polas.size();i++) {
		SingleAtom sa = polas[i].second;
		expose_.insert({sa,0});
		std::vector<XYZ> ps1 = ps;
		for(auto &p:ps1) p = p + polas[i].first;
		for(XYZ &c:ps1) {
			bool out{true};
			for(auto &p:allas) {
				if(p.second == sa) continue;
				if((p.first-c).squarednorm()>r2) continue;
				out = false;
				break;
			}
			if(out) expose_.at(sa)++;
		}
	}
}

void Filter::allatomexpose() {
	expose_.clear();
	double wrad = 1.6;
	std::map<std::string,std::map<std::string,double>> rs = AA_Atom::radius3();
	for(auto &r:rs) {
		for(auto &p:r.second) p.second = 1.9;
	}
	std::vector<XYZ> crds;
	std::vector<double> rads;
	std::vector<std::string> rns;
	std::vector<std::string> ans;
	std::vector<int> chs;
	std::vector<int> rids;
	for(int i=0;i<chains_.size();i++) {
		for(int j=0;j<chains_[i].size();j++) {
			for(auto &p:chains_[i][j].rds()) {
				std::string rn = chains_[i][j].resname();
				std::string an = p.first;
				if(rs.find(rn)==rs.end()) continue;
				if(rs.at(rn).find(an)==rs.at(rn).end()) continue;
				crds.push_back(p.second.crd);
				//std::cout <<rn <<'\t' <<an <<std::endl;
				rads.push_back(rs.at(rn).at(an)+wrad);
				rns.push_back(rn);
				ans.push_back(an);
				chs.push_back(i);
				rids.push_back(j);
				expose_.insert({{{i,j},an},0});
			}
		}
	}
	std::vector<XYZ> ps = AA_Atom::sphere256();
	for(int i=0;i<crds.size();i++) {
		std::vector<XYZ> cs = ps;
		for(XYZ &c:cs) c = c*rads[i] + crds[i];
		SingleAtom sa = {{chs[i],rids[i]},ans[i]};
		for(XYZ &c:cs) {
			bool out{true};
			for(int j=0;j<crds.size();j++) {
				if(i==j) continue;
				double dis2 = rads[j] * rads[j];
				if((c-crds[j]).squarednorm()>dis2) continue;
				out = false;
				break;
			}
			if(out) expose_.at(sa) ++;
		}
	}
}

void Filter::allatomexpose_psudosc() {
	expose_.clear();
	std::vector<XYZ> crds;
	std::vector<std::string> rns;
	std::vector<std::string> ans;
	std::vector<int> chs;
	std::vector<int> rids;
	for(int i=0;i<chains_.size();i++) {
		for(int j=0;j<chains_[i].size();j++) {
			for(auto &p:chains_[i][j].rds()) {
				std::string rn = chains_[i][j].resname();
				std::string an = p.first;
				crds.push_back(p.second.crd);
				rns.push_back(rn);
				ans.push_back(an);
				chs.push_back(i);
				rids.push_back(j);
				expose_.insert({{{i,j},an},0});
			}
		}
	}

	std::vector<double> rads(crds.size(),3.5);
	std::vector<XYZ> ps = AA_Atom::sphere256();
	for(int i=0;i<crds.size();i++) {
		std::vector<XYZ> cs = ps;
		for(XYZ &c:cs) c = c*rads[i] + crds[i];
		SingleAtom sa = {{chs[i],rids[i]},ans[i]};
		for(XYZ &c:cs) {
			bool out{true};
			for(int j=0;j<crds.size();j++) {
				if(i==j) continue;
				double dis2 = rads[j] * rads[j];
				if((c-crds[j]).squarednorm()>dis2) continue;
				out = false;
				break;
			}
			if(out) expose_.at(sa) ++;
		}
	}
}

void Filter::expose_remove_self() {
	expose_.clear();
	double wrad = 1.6;
	std::map<std::string,std::map<std::string,double>> rs = AA_Atom::radius3();
	for(auto &r:rs) {
		for(auto &p:r.second) p.second = 1.9;
	}
	std::map<int,std::map<int,std::map<std::pair<std::string,double>,XYZ>>> dts;
	for(int i=0;i<chains_.size();i++) {
		dts.insert({i,std::map<int,std::map<std::pair<std::string,double>,XYZ>>()});
		for(int j=0;j<chains_[i].size();j++) {
			dts.at(i).insert({j,std::map<std::pair<std::string,double>,XYZ>()});
			for(auto &p:chains_[i][j].rds()) {
				std::string rn = chains_[i][j].resname();
				std::string an = p.first;
				if(rs.find(rn)==rs.end()) continue;
				if(rs.at(rn).find(an)==rs.at(rn).end()) continue;
				double rad=rs.at(rn).at(an)+wrad;
				dts.at(i).at(j).insert({{an,rad},p.second.crd});
			}
		}
	}

	std::vector<XYZ> ps = AA_Atom::sphere256();

	for(auto &p1:dts) {
		for(auto &p2:p1.second) {
			for(auto &p3:p2.second) {
				std::vector<std::pair<XYZ,double>> vp;
				for(auto &p4:p2.second) {
					if(p4.first.first==p3.first.first) continue;
					vp.push_back({p4.second,p4.first.second});
				}
				std::set<std::string> mcs{"C","CA","N","O"};
				if(p1.second.find(p2.first-1)!=p1.second.end()) {
					auto &mp=p1.second.at(p2.first-1);
					for(auto &p4:mp) {
						if(mcs.find(p4.first.first)==mcs.end()) continue;
						vp.push_back({p4.second,p4.first.second});
					}
				}
				if(p1.second.find(p2.first+1)!=p1.second.end()) {
					auto &mp=p1.second.at(p2.first+1);
					for(auto &p4:mp) {
						if(mcs.find(p4.first.first)==mcs.end()) continue;
						vp.push_back({p4.second,p4.first.second});
					}
				}
				std::vector<XYZ> cs = ps;
				for(XYZ &c:cs) c = c*p3.first.second + p3.second;
				SingleAtom sa = {{p1.first,p2.first},p3.first.first};
				expose_.insert({sa,0});
				//std::cout <<vp.size() <<std::endl;
				for(XYZ &c:cs) {
					bool out{true};
					for(auto &p4:vp) {
						double dis2=(c-p4.first).squarednorm();
						double dis22=p4.second*p4.second;
						if(dis2>dis22) continue;
						out=false;
						break;
					}
					if(!out) continue;
					expose_.at(sa)++;
				}
				//std::cout <<expose_.at(sa) <<std::endl;
			}
		}
	}
}

//std::vector<std::map<std::string,std::vector<int>>> Filter::mcpolar() {
std::vector<std::map<std::string,std::vector<std::pair<SingleAtom,int>>>> Filter::mcpolar() {
	hbps_.clear();
	gethbps(false);
	expose_.clear();
	mcpolarexpose();
	//if(hbps_.empty()) gethbps(false);
	//if(expose_.empty()) mcpolarexpose();
	std::set<SingleAtom> sas;
	for(auto &vh:hbps_) {
		for(HB &h:vh) {
			for(auto &p:h.hbs) {
				if(p.second.empty()) continue;
				sas.insert({{h.chainid,h.resid},p.first});
				for(auto &sa1:p.second) sas.insert(sa1);
			}
		}
	}
	std::map<std::string,std::vector<int>> mp;
	mp.insert({"N",std::vector<int>()});
	mp.insert({"O",std::vector<int>()});
	std::vector<std::map<std::string,std::vector<int>>> mcp(2,mp);

	std::map<std::string,std::vector<std::pair<SingleAtom,int>>> mpa;
	mpa.insert({"N",std::vector<std::pair<SingleAtom,int>>()});
	mpa.insert({"O",std::vector<std::pair<SingleAtom,int>>()});
	std::vector<std::map<std::string,std::vector<std::pair<SingleAtom,int>>>> mcpa(2,mpa);

	for(auto &p:expose_) {
		if(sas.find(p.first)==sas.end()) {
			mcp[0].at(p.first.second).push_back(p.second);
			mcpa[0].at(p.first.second).push_back(p);
		} else {
			mcp[1].at(p.first.second).push_back(p.second);
			mcpa[1].at(p.first.second).push_back(p);
		}
	}
	return mcpa;
}

std::set<std::vector<SingleAtom>> Filter::gethbpair(bool diffpro) {
	if(hbps_.empty()) gethbps(diffpro);
	std::set<std::vector<SingleAtom>> sas;
	for(auto &vh:hbps_) {
		for(HB &h:vh) {
			for(auto &p:h.hbs) {
				SingleAtom sa1{{h.chainid,h.resid},p.first};
				for(auto &sa2:p.second) sas.insert({sa1,sa2});
			}
		}
	}
	return sas;
}

/*std::vector<std::map<std::vector<std::string>,std::vector<int>>> Filter::expose(bool considerHB) {
	if(considerHB) {
		hbps_.clear();
		gethbps();
	}
	expose_.clear();
	allatomexpose();
	std::map<SingleResidue,std::string> rns;
	for(int i=0;i<chains_.size();i++) {
		for(int j=0;j<chains_[i].size();j++) {
			rns.insert({{i,j},chains_[i][j].resname()});
		}
	}
	std::set<SingleAtom> sas;
	for(auto &vh:hbps_) {
		for(HB &h:vh) {
			for(auto &p:h.hbs) {
				if(p.second.empty()) continue;
				sas.insert({{h.chainid,h.resid},p.first});
			}
		}
	}

	std::map<std::vector<std::string>,std::vector<int>> mp;
	std::map<std::string,std::set<std::string>> scs = AA_Atom::sidechain3();
	std::set<std::string> mcs{"C","CA","N","O"};
	for(auto &ps:scs) {
		for(auto &pm:mcs) mp.insert({{ps.first,pm},std::vector<int>()});
		for(auto &s:ps.second) mp.insert({{ps.first,s},std::vector<int>()});
	}
	std::vector<std::map<std::vector<std::string>,std::vector<int>>> mcp(1,mp);
	if(considerHB) mcp.push_back(mp);

	//std::cout <<mcp[0].size() <<'\t' <<mcp[1].size() <<std::endl;

	for(auto &p:expose_) {
		std::vector<std::string> vs{rns.at(p.first.first),p.first.second};
		//std::cout <<vs[0] <<'\t' <<vs[1] <<std::endl;
		if(considerHB) {
			if(sas.find(p.first)==sas.end()) { // no HB
				mcp[0].at(vs).push_back(p.second);
				//mcpa[0].at(p.first.second).push_back(p);
			} else {
				mcp[1].at(vs).push_back(p.second);
				//mcpa[1].at(p.first.second).push_back(p);
			}
		} else {
			mcp[0].at(vs).push_back(p.second);
		}
	}
	return mcp;
}*/

bool NSPallatom::mainchainHBfilter(const std::vector<std::vector<Residue>>&chs,
		const std::vector<std::vector<int>>&helixs, const std::vector<std::vector<int>>&strands) {
	Filter fl(chs);
	std::vector<std::vector<double>> mchb=fl.mainchainHBfilter(helixs, strands);
	for(int j=0;j<mchb[0].size();j++) {
		if(mchb[0][j]>0.7) continue;
		return false;
		break;
	}
	for(int j=0;j<mchb[1].size();j++) {
		if(mchb[1][j]>0.2) continue;
		return false;
		break;
	}
	return true;
}
















bool EnergyTermRecord::allatomismc(const MultiAtom &ma) {
	static std::set<std::string> mc{"C", "CA", "N", "O"};
	for(const SingleAtom &sa:ma) {
		if(mc.find(sa.second)==mc.end()) return false;
	}
	return true;
}

std::pair<MultiAtom,double> EnergyTermRecord::readmultiatom(std::stringstream &ss, int na) {
	std::pair<MultiAtom,double> p;
	for(int i=0;i<na;i++) {
		int c, r;
		std::string a;
		ss >> c >> r >> a;
		p.first.push_back({{c,r},a});
	}
	ss >> p.second;
	return p;
}

std::pair<MultiResidue,double> EnergyTermRecord::readmultiresidue(std::stringstream &ss) {
	std::pair<MultiResidue,double> p;
	SingleResidue s1(2), s2(2);
	ss >> s1[0] >> s1[1] >> s2[0] >> s2[1];
	p.first.push_back(s1);
	p.first.push_back(s2);
	ss >> p.second;
	return p;
}

std::pair<SingleResidue,double> EnergyTermRecord::readsingleresidue(std::stringstream &ss) {
	std::pair<SingleResidue,double> p;
	p.first.resize(2);
	ss >> p.first[0] >> p.first[1] >> p.second;
	return p;
}

void EnergyTermRecord::readenergy(std::string enefile, bool mconly) {
	std::ifstream ifs(enefile);
	std::string readline;
	std::stringstream ss;
	while(std::getline(ifs,readline)) {
		std::string label;
		int nline;
		ss << readline;
		ss >> label >> nline;
		ss.clear();
		if(label=="bond") {
			for(int i=0;i<nline;i++) {
				std::getline(ifs,readline);
				ss << readline;
				std::pair<MultiAtom,double> md = readmultiatom(ss,2);
				ss.clear();
				if(mconly) {
					if(!allatomismc(md.first)) continue;
				}
				ebond.insert(md);
			}
		} else if(label=="ang") {
			for(int i=0;i<nline;i++) {
				std::getline(ifs,readline);
				ss << readline;
				std::pair<MultiAtom,double> md = readmultiatom(ss,3);
				//eang.insert(readmultiatom(ss,3));
				ss.clear();
				if(mconly) {
					if(!allatomismc(md.first)) continue;
				}
				eang.insert(md);
			}
		} else if(label=="impdih") {
			for(int i=0;i<nline;i++) {
				std::getline(ifs,readline);
				ss << readline;
				std::pair<MultiAtom,double> md = readmultiatom(ss,4);
				//eimpdih.insert(readmultiatom(ss,4));
				ss.clear();
				if(mconly) {
					if(!allatomismc(md.first)) continue;
				}
				eimpdih.insert(md);
			}
		} else if(label=="steric") {
			for(int i=0;i<nline;i++) {
				std::getline(ifs,readline);
				ss << readline;
				std::pair<MultiAtom,double> md = readmultiatom(ss,2);
				//esteric.insert(readmultiatom(ss,2));
				ss.clear();
				if(mconly) {
					if(!allatomismc(md.first)) continue;
				}
				esteric.insert(md);
			}
		} else if(label=="sitepair") {
			for(int i=0;i<nline;i++) {
				std::getline(ifs,readline);
				ss << readline;
				esitepair.insert(readmultiresidue(ss));
				ss.clear();
			}
		} else if(label=="localbbhb") {
			for(int i=0;i<nline;i++) {
				std::getline(ifs,readline);
				ss << readline;
				elocalbbhb.insert(readmultiresidue(ss));
				ss.clear();
			}
		} else if(label=="scpacking") {
			for(int i=0;i<nline;i++) {
				std::getline(ifs,readline);
				ss << readline;
				escpacking.insert(readmultiresidue(ss));
				ss.clear();
			}
		} else if(label=="phipsi") {
			for(int i=0;i<nline;i++) {
				std::getline(ifs,readline);
				ss << readline;
				ephipsi.insert(readsingleresidue(ss));
				ss.clear();
			}
		} else if(label=="local") {
			for(int i=0;i<nline;i++) {
				std::getline(ifs,readline);
				ss << readline;
				elocal.insert(readsingleresidue(ss));
				ss.clear();
			}
		} else if(label=="scconf") {
			for(int i=0;i<nline;i++) {
				std::getline(ifs,readline);
				ss << readline;
				escconf.insert(readsingleresidue(ss));
				ss.clear();
			}
		} else if(label=="total") {
			for(int i=0;i<nline;i++) {
				std::getline(ifs,readline);
				ss << readline;
				etot.insert(readsingleresidue(ss));
				ss.clear();
			}
		} else {
			std::cout <<"Can Not Recognize Energy :   " <<label <<std::endl;
		}
	}
	ifs.close();
}

std::map<SingleResidue,double> EnergyTermRecord::singleresidueenergy(const std::map<MultiResidue,double>&mes) {
	std::map<SingleResidue,double> ses;
	for(const auto &p:mes) {
		double e = p.second / 2.0;
		const SingleResidue& sr0 = p.first[0];
		const SingleResidue& sr1 = p.first[1];
		if(ses.find(sr0)==ses.end()) ses.insert({sr0,0.0});
		if(ses.find(sr1)==ses.end()) ses.insert({sr1,0.0});
		ses.at(sr0) += e;
		ses.at(sr1) += e;
	}
	return ses;
}

void EnergyTermRecord::printene(std::string outene) {
	std::ofstream ofs(outene);
	ofs <<"bond\t" <<ebond.size() <<std::endl;
	for(auto & p : ebond) {
		auto & ma = p.first;
		for(auto & sa : ma) {
			ofs <<sa.first[0] <<'\t' <<sa.first[1] <<'\t' <<sa.second <<'\t';
		}
		ofs <<p.second <<std::endl;
	}
	ofs <<"ang\t" <<eang.size() <<std::endl;
	for(auto & p : eang) {
		auto & ma = p.first;
		for(auto & sa : ma) {
			ofs <<sa.first[0] <<'\t' <<sa.first[1] <<'\t' <<sa.second <<'\t';
		}
		ofs <<p.second <<std::endl;
	}
	ofs <<"impdih\t" <<eimpdih.size() <<std::endl;
	for(auto & p : eimpdih) {
		auto & ma = p.first;
		for(auto & sa : ma) {
			ofs <<sa.first[0] <<'\t' <<sa.first[1] <<'\t' <<sa.second <<'\t';
		}
		ofs <<p.second <<std::endl;
	}
	ofs <<"steric\t" <<esteric.size() <<std::endl;
	for(auto & p : esteric) {
		auto & ma = p.first;
		for(auto & sa : ma) {
			ofs <<sa.first[0] <<'\t' <<sa.first[1] <<'\t' <<sa.second <<'\t';
		}
		ofs <<p.second <<std::endl;
	}
	ofs <<"sitepair\t" <<esitepair.size() <<std::endl;
	for(auto & p : esitepair) {
		auto & mr = p.first;
		for(auto & sr : mr) {
			ofs <<sr[0] <<'\t' <<sr[1] <<'\t';
		}
		ofs <<p.second <<std::endl;
	}
	ofs <<"localbbhb\t" <<elocalbbhb.size() <<std::endl;
	for(auto & p : elocalbbhb) {
		auto & mr = p.first;
		for(auto & sr : mr) {
			ofs <<sr[0] <<'\t' <<sr[1] <<'\t';
		}
		ofs <<p.second <<std::endl;
	}
	ofs <<"scpacking\t" <<escpacking.size() <<std::endl;
	for(auto & p : escpacking) {
		auto & mr = p.first;
		for(auto & sr : mr) {
			ofs <<sr[0] <<'\t' <<sr[1] <<'\t';
		}
		ofs <<p.second <<std::endl;
	}
	ofs <<"phipsi\t" <<ephipsi.size() <<std::endl;
	for(auto & p : ephipsi) {
		ofs <<p.first[0] <<'\t' <<p.first[1] <<'\t' <<p.second <<std::endl;
	}
	ofs <<"local\t" <<elocal.size() <<std::endl;
	for(auto & p : elocal) {
		ofs <<p.first[0] <<'\t' <<p.first[1] <<'\t' <<p.second <<std::endl;
	}
	ofs <<"scconf\t" <<escconf.size() <<std::endl;
	for(auto & p : escconf) {
		ofs <<p.first[0] <<'\t' <<p.first[1] <<'\t' <<p.second <<std::endl;
	}
	ofs <<"total\t" <<etot.size() <<std::endl;
	for(auto & p : etot) {
		ofs <<p.first[0] <<'\t' <<p.first[1] <<'\t' <<p.second <<std::endl;
	}
	ofs.close();
}

void EnergyTermRecord::eneinit(NSPsd::EnergyTerm &et) {
	for(int i=0;i<decoy.size();i++) {
		for(int j=0;j<decoy[i].size();j++) {
			etot.insert({{i,j},0.0});
		}
	}
	for(auto &e:et.ebond) {
		rank(e.first);
		ebond.insert(e);
	}
	for(auto &e:et.eang) {
		rank(e.first);
		eang.insert(e);
	}
	for(auto &e:et.eimpdih) {
		rank(e.first);
		eimpdih.insert(e);
	}
	for(auto &e:et.esteric) {
		rank(e.first);
		esteric.insert(e);
	}
	for(auto &e:et.esitepair) {
		rank(e.first);
		esitepair.insert(e);
	}
	for(auto &e:et.elocalbbhb) {
		rank(e.first);
		elocalbbhb.insert(e);
	}
	for(auto &e:et.escpacking) {
		rank(e.first);
		escpacking.insert(e);
	}
	for(auto &e:et.ephipsi) ephipsi.insert(e);
	for(auto &e:et.elocal) elocal.insert(e);
	for(auto &e:et.escconf) escconf.insert(e);

	for(auto &p:ebond) {
		const MultiAtom &ma=p.first;
		for(const SingleAtom &sa:ma) {
			if(etot.find(sa.first)==etot.end()) etot.insert({sa.first,0.0});
			etot.at(sa.first) += p.second / 2.0;
		}
	}
	for(auto &p:eang) {
		const MultiAtom &ma=p.first;
		for(const SingleAtom &sa:ma) {
			if(etot.find(sa.first)==etot.end()) etot.insert({sa.first,0.0});
			etot.at(sa.first) += p.second / 3.0;
		}
	}
	for(auto &p:eimpdih) {
		const MultiAtom &ma=p.first;
		for(const SingleAtom &sa:ma) {
			if(etot.find(sa.first)==etot.end()) etot.insert({sa.first,0.0});
			etot.at(sa.first) += p.second / 4.0;
		}
	}
	for(auto &p:esteric) {
		const MultiAtom &ma=p.first;
		for(const SingleAtom &sa:ma) {
			if(etot.find(sa.first)==etot.end()) etot.insert({sa.first,0.0});
			etot.at(sa.first) += p.second / 2.0;
		}
	}
	for(auto &p:esitepair) {
		for(const SingleResidue &sr:p.first) {
			if(etot.find(sr)==etot.end()) etot.insert({sr,0.0});
			etot.at(sr) += p.second / 2.0;
		}
	}
	for(auto &p:elocalbbhb) {
		for(const SingleResidue &sr:p.first) {
			if(etot.find(sr)==etot.end()) etot.insert({sr,0.0});
			etot.at(sr) += p.second / 2.0;
		}
	}
	for(auto &p:escpacking) {
		for(const SingleResidue &sr:p.first) {
			if(etot.find(sr)==etot.end()) etot.insert({sr,0.0});
			etot.at(sr) += p.second / 2.0;
		}
	}
	for(auto &p:ephipsi) {
		if(etot.find(p.first)==etot.end()) etot.insert({p.first,0.0});
		etot.at(p.first) += p.second;
	}
	for(auto &p:elocal) {
		if(etot.find(p.first)==etot.end()) etot.insert({p.first,0.0});
		etot.at(p.first) += p.second;
	}
	for(auto &p:escconf) {
		if(etot.find(p.first)==etot.end()) etot.insert({p.first,0.0});
		etot.at(p.first) += p.second;
	}
}

double EnergyTermRecord::totalene() {
	double te = 0.0;
	for(auto &p:ebond) te += p.second;
	for(auto &p:eang) te += p.second;
	for(auto &p:eimpdih) te += p.second;
	for(auto &p:esteric) te += p.second;
	for(auto &p:ephipsi) te += p.second;
	for(auto &p:elocal) te += p.second;
	for(auto &p:escconf) te += p.second;
	for(auto &p:esitepair) te += p.second;
	for(auto &p:elocalbbhb) te += p.second;
	for(auto &p:escpacking) te += p.second;
	return te;
}

EnergyTermRecord::EnergyTermRecord(const DECOY &dy):decoy(dy) {
	std::string defaultffparfile=NSPdataio::datafilename("DefaultFF.par");
	NSPdataio::ControlFile cf;
	cf.readfile(defaultffparfile);
	std::vector<std::string> ffcontrolines=cf.getcontrolines("ForceField");
	std::string ffcontrolname{"control_ff"};
	NSPsd::defineforcefieldcontrol(ffcontrolname,ffcontrolines);
	std::vector<std::vector<NSPproteinrep::FullSite>> fsss(decoy.size());
	for(int i=0;i<decoy.size();i++) {
		for(int j=0;j<decoy[i].size();j++) {
			fsss[i].push_back(decoy[i][j].getfullsite());
		}
	}
	NSPsd::ForceField ff=NSPsd::make_forcefield_allatom(fsss, ffcontrolname);

	std::string allfixed, mcfixed;
	NSPsd::ActiveSelections acts(&ff,allfixed,mcfixed);

	std::vector<std::pair<std::vector<int>,std::string>> atomseq;
	std::vector<double> cs=NSPsd::GenChain::getcrd(fsss,atomseq);
	for(double & c:cs) c=c*A2NM;
	NSPsd::NeighborList nbl(cs,ff);
	std::vector<double> pots;
	std::vector<int> sizes;
	for(int i=0;i<decoy.size();i++) sizes.push_back(decoy[i].size());
	NSPsd::EnergyTerm etm;
	ff.forces_subentry(cs, nbl, &pots, acts, atomseq, etm, sizes);
	eneinit(etm);
}
/*
void EnergyTermRecord::reads1s2(std::string efile) {
	std::string readline;
	std::stringstream ss;
	std::ifstream ifs(efile);
	std::getline(ifs,readline);
	while(std::getline(ifs,readline)) {
		if(readline[0]=='T') break;
		std::string aa;
		char sec;
		double s1,s2,d;
		ss <<readline;
		ss >>aa >>aa >>sec >>s1 >>s1 >>s2 >>d >>d;
		ss.clear();
		es1.push_back({{sec,aa},s1});
		es2.push_back({{sec,aa},s2});*/
		/*if(es1.find(sec)==es1.end()) es1.at(sec).insert({sec,std::map<std::string,std::vector<double>>()});
		if(es1.at(sec).find(aa)==es1.at(sec).end()) es1.at(sec).insert({aa,std::vector<double>()});
		es1.at(sec).at(aa).push_back(s1);
		if(es2.find(sec)==es2.end()) es2.at(sec).insert({sec,std::map<std::string,std::vector<double>>()});
		if(es2.at(sec).find(aa)==es2.at(sec).end()) es2.at(sec).insert({aa,std::vector<double>()});
		es2.at(sec).at(aa).push_back(s2);*/
//	}
//	ifs.close();
//}








void TrajectoryEnergy::changeweight(std::vector<std::string>&ctrl, std::string label, double wgt) {
	int sz = label.size();
	for(auto &l:ctrl) {
		if(l.substr(0,sz)==label) {
			l = label + "\t=\t" + std::to_string(wgt);
			return;
		}
	}
	std::cout <<"Can Not Find Energy Term in Control Lines!" <<std::endl;
	std::cout <<label <<std::endl;
}

void TrajectoryEnergy::dosd(int nstep, int seed,
		std::map<std::string,double> &sdctrl, std::map<std::string,double> &ffctrl) {
	NSPdstl::RandomEngine<>::getinstance().reseed(seed);
	std::vector<std::vector<NSPproteinrep::FullSite>> fsss;
	std::vector<int> sizes;
	for(auto &ch:decoy_) {
		fsss.push_back(std::vector<NSPproteinrep::FullSite>());
		sizes.push_back(ch.size());
		for(auto &r:ch) {
			fsss.back().push_back(r.getfullsite());
		}
	}
	std::string controlname{"sdffcontrol"};
	std::string controlfile=NSPdataio::datafilename("FilterSDControl.par");
	NSPdataio::ControlFile cf;
	cf.readfile(controlfile);
	std::vector<std::string> sdcontrolines=cf.getcontrolines("SD");
	std::vector<std::string> ffcontrolines=cf.getcontrolines("ForceField");
	for(auto &p:sdctrl) {
		if(p.first == "temperature") changeweight(sdcontrolines,"Temperatures",p.second);
		else if(p.first == "friction") changeweight(sdcontrolines,"FrictionCoeff",p.second);
		else {
			std::cout <<"Can Not Read Energy Term in Input Line!" <<std::endl;
			std::cout <<p.first <<std::endl;
		}
	}
	for(auto &p:ffctrl) {
		if(p.first == "steric") changeweight(ffcontrolines,"StericWeight",p.second);
		else if(p.first == "phipsi") changeweight(ffcontrolines,"PhiPsiWeight",p.second);
		else if(p.first == "local") changeweight(ffcontrolines,"LSWeight",p.second);
		else if(p.first == "scconf") changeweight(ffcontrolines,"SideConfWeight",p.second);
		else if(p.first == "sitepair") changeweight(ffcontrolines,"PairPackingWeight",p.second);
		else if(p.first == "localbbhb") changeweight(ffcontrolines,"LocalHBWeight",p.second);
		else if(p.first == "scpacking") changeweight(ffcontrolines,"SCPackingWeight",p.second);
		else {
			std::cout <<"Can Not Read Energy Term in Input Line!" <<std::endl;
			std::cout <<p.first <<std::endl;
		}
	}
	//for(auto &s:sdcontrolines) std::cout <<s <<std::endl;
	//for(auto &s:ffcontrolines) std::cout <<s <<std::endl;
	NSPsd::definesdcontrol(controlname+"_sd",sdcontrolines);
	NSPsd::defineforcefieldcontrol(controlname+"_ff",ffcontrolines);
	NSPsd::GenChain genchain(controlname+"_genchain",fsss);
	std::vector<std::pair<std::vector<int>,std::string>> atomseq;
	std::vector<double> initcrd = genchain.getcrd(fsss,atomseq);
	for (auto &c : initcrd) c *= A2NM;
	NSPsd::SDRun::SDRunIn sdrunin(initcrd,controlname);
	std::vector<std::set<int>> notcis(fsss.size());
	std::string allfixed,mcfixed;
	for(auto &f:fixed_) {
		std::string s = "chain" + std::to_string(f[0])+" "+std::to_string(f[1])+"-"+std::to_string(f[2])+" ";
		allfixed = allfixed + s;
	}
	//std::cout <<allfixed <<std::endl;
	NSPsd::SDRun sdrun = genchain.make_sdrun(sdrunin,seed,allfixed,mcfixed,notcis);
	int iterval = 20;
	for(int i=0;i<nstep/iterval;i++) {
		std::cout <<"SD :   " <<iterval*i <<"  /  " <<nstep <<std::endl;
		NSPsd::EnergyTerm et;
		if(!(sdrun.runsteps_subentry(iterval,atomseq,et,sizes))){
			std::cout<<"Shake failure occurred."<<std::endl;
			exit(1);
		}
		std::vector<std::vector<NSPproteinrep::FullSite>> fs1 = genchain.crd2fullsite(sdrun.state().crd, 1.0/A2NM);
		DECOY dy;
		for(auto &fss:fs1) {
			dy.push_back(std::vector<Residue>());
			for(auto &fs:fss) {
				dy.back().push_back(Residue(fs));
			}
		}
		eneterm_.push_back(EnergyTermRecord(dy,et));
	}
}

void TrajectoryEnergy::energyanalysis(int mintime) {
	HighEnergyTime & het = het_;
	for(int i=0;i<eneterm_.size();i++) {
		for(auto &p:eneterm_[i].ebond) {
			auto & e = em_.ebond;
			if(e.find(p.first)==e.end()) e.insert({p.first,std::vector<double>()});
			e.at(p.first).push_back(p.second);
			if(p.second<cut0[ENE::BOND]) continue;
			auto & t = het.tbond;
			if(t.find(p.first)==t.end()) t.insert({p.first,std::vector<int>()});
			t.at(p.first).push_back(i);
		}
		for(auto &p:eneterm_[i].eang) {
			auto & e = em_.eang;
			if(e.find(p.first)==e.end()) e.insert({p.first,std::vector<double>()});
			e.at(p.first).push_back(p.second);
			if(p.second<cut0[ENE::ANG]) continue;
			auto & t = het.tang;
			if(t.find(p.first)==t.end()) t.insert({p.first,std::vector<int>()});
			t.at(p.first).push_back(i);
		}
		for(auto &p:eneterm_[i].eimpdih) {
			auto & e = em_.eimpdih;
			if(e.find(p.first)==e.end()) e.insert({p.first,std::vector<double>()});
			e.at(p.first).push_back(p.second);
			if(p.second<cut0[ENE::IMPDIH]) continue;
			auto & t = het.timpdih;
			if(t.find(p.first)==t.end()) t.insert({p.first,std::vector<int>()});
			t.at(p.first).push_back(i);
		}
		for(auto &p:eneterm_[i].esteric) {
			auto & e = em_.esteric;
			if(e.find(p.first)==e.end()) e.insert({p.first,std::vector<double>()});
			e.at(p.first).push_back(p.second);
			if(p.second<cut0[ENE::STERIC]) continue;
			auto & t = het.tsteric;
			if(t.find(p.first)==t.end()) t.insert({p.first,std::vector<int>()});
			t.at(p.first).push_back(i);
		}
		for(auto &p:eneterm_[i].ephipsi) {
			auto & e = em_.ephipsi;
			if(e.find(p.first)==e.end()) e.insert({p.first,std::vector<double>()});
			e.at(p.first).push_back(p.second);
			if(p.second<cut0[ENE::PHIPSI]) continue;
			auto & t = het.tphipsi;
			if(t.find(p.first)==t.end()) t.insert({p.first,std::vector<int>()});
			t.at(p.first).push_back(i);
		}
		for(auto &p:eneterm_[i].elocal) {
			auto & e = em_.elocal;
			if(e.find(p.first)==e.end()) e.insert({p.first,std::vector<double>()});
			e.at(p.first).push_back(p.second);
			if(p.second<cut0[ENE::LOCAL]) continue;
			auto & t = het.tlocal;
			if(t.find(p.first)==t.end()) t.insert({p.first,std::vector<int>()});
			t.at(p.first).push_back(i);
		}
		for(auto &p:eneterm_[i].escconf) {
			auto & e = em_.escconf;
			if(e.find(p.first)==e.end()) e.insert({p.first,std::vector<double>()});
			e.at(p.first).push_back(p.second);
			if(p.second<cut0[ENE::SCCONF]) continue;
			auto & t = het.tscconf;
			if(t.find(p.first)==t.end()) t.insert({p.first,std::vector<int>()});
			t.at(p.first).push_back(i);
		}
		for(auto &p:eneterm_[i].esitepair) {
			auto & e = em_.esitepair;
			if(e.find(p.first)==e.end()) e.insert({p.first,std::vector<double>()});
			e.at(p.first).push_back(p.second);
			if(p.second<cut0[ENE::SITEPAIR]) continue;
			auto & t = het.tsitepair;
			if(t.find(p.first)==t.end()) t.insert({p.first,std::vector<int>()});
			t.at(p.first).push_back(i);
		}
		for(auto &p:eneterm_[i].elocalbbhb) {
			auto & e = em_.elocalbbhb;
			if(e.find(p.first)==e.end()) e.insert({p.first,std::vector<double>()});
			e.at(p.first).push_back(p.second);
			if(p.second<cut0[ENE::LOCALBBHB]) continue;
			auto & t = het.tlocalbbhb;
			if(t.find(p.first)==t.end()) t.insert({p.first,std::vector<int>()});
			t.at(p.first).push_back(i);
		}
		for(auto &p:eneterm_[i].escpacking) {
			auto & e = em_.escpacking;
			if(e.find(p.first)==e.end()) e.insert({p.first,std::vector<double>()});
			e.at(p.first).push_back(p.second);
			if(p.second<cut0[ENE::SCPACKING]) continue;
			auto & t = het.tscpacking;
			if(t.find(p.first)==t.end()) t.insert({p.first,std::vector<int>()});
			t.at(p.first).push_back(i);
		}
		for(auto &p:eneterm_[i].etot) {
			auto & e = em_.etot;
			if(e.find(p.first)==e.end()) e.insert({p.first,std::vector<double>()});
			e.at(p.first).push_back(p.second);
			if(p.second<cut0[ENE::TOTAL]) continue;
			auto & t = het.ttot;
			if(t.find(p.first)==t.end()) t.insert({p.first,std::vector<int>()});
			t.at(p.first).push_back(i);
		}
	}

	for(auto &p:em_.ebond) {
		if(p.second.size()<mintime) continue;
		double sum = 0.0;
		for(double d:p.second) sum += d;
		sum /= (double)(p.second.size());
		ea_.ebond.insert({p.first,sum});
		if(sum>cut0[ENE::BOND]) he_.ebond.insert({p.first,sum});
	}
	for(auto &p:em_.eang) {
		if(p.second.size()<mintime) continue;
		double sum = 0.0;
		for(double d:p.second) sum += d;
		sum /= (double)(p.second.size());
		ea_.eang.insert({p.first,sum});
		if(sum>cut0[ENE::ANG]) he_.eang.insert({p.first,sum});
	}
	for(auto &p:em_.eimpdih) {
		if(p.second.size()<mintime) continue;
		double sum = 0.0;
		for(double d:p.second) sum += d;
		sum /= (double)(p.second.size());
		ea_.eimpdih.insert({p.first,sum});
		if(sum>cut0[ENE::IMPDIH]) he_.eimpdih.insert({p.first,sum});
	}
	for(auto &p:em_.esteric) {
		if(p.second.size()<mintime) continue;
		double sum = 0.0;
		for(double d:p.second) sum += d;
		sum /= (double)(p.second.size());
		ea_.esteric.insert({p.first,sum});
		if(sum>cut0[ENE::STERIC]) he_.esteric.insert({p.first,sum});
	}
	for(auto &p:em_.ephipsi) {
		if(p.second.size()<mintime) continue;
		double sum = 0.0;
		for(double d:p.second) sum += d;
		sum /= (double)(p.second.size());
		ea_.ephipsi.insert({p.first,sum});
		if(sum>cut0[ENE::PHIPSI]) he_.ephipsi.insert({p.first,sum});
	}
	for(auto &p:em_.elocal) {
		if(p.second.size()<mintime) continue;
		double sum = 0.0;
		for(double d:p.second) sum += d;
		sum /= (double)(p.second.size());
		ea_.elocal.insert({p.first,sum});
		if(sum>cut0[ENE::LOCAL]) he_.elocal.insert({p.first,sum});
	}
	for(auto &p:em_.escconf) {
		if(p.second.size()<mintime) continue;
		double sum = 0.0;
		for(double d:p.second) sum += d;
		sum /= (double)(p.second.size());
		ea_.escconf.insert({p.first,sum});
		if(sum>cut0[ENE::SCCONF]) he_.escconf.insert({p.first,sum});
	}
	for(auto &p:em_.etot) {
		if(p.second.size()<mintime) continue;
		double sum = 0.0;
		for(double d:p.second) sum += d;
		sum /= (double)(p.second.size());
		ea_.etot.insert({p.first,sum});
		if(sum>cut0[ENE::TOTAL]) he_.etot.insert({p.first,sum});
	}
	for(auto &p:em_.esitepair) {
		if(p.second.size()<mintime) continue;
		double sum = 0.0;
		for(double d:p.second) sum += d;
		sum /= (double)(p.second.size());
		ea_.esitepair.insert({p.first,sum});
		if(sum>cut0[ENE::SITEPAIR]) he_.esitepair.insert({p.first,sum});
	}
	for(auto &p:em_.escpacking) {
		if(p.second.size()<mintime) continue;
		double sum = 0.0;
		for(double d:p.second) sum += d;
		sum /= (double)(p.second.size());
		ea_.escpacking.insert({p.first,sum});
		if(sum>cut0[ENE::SCPACKING]) he_.escpacking.insert({p.first,sum});
	}
	for(auto &p:em_.elocalbbhb) {
		if(p.second.size()<mintime) continue;
		double sum = 0.0;
		for(double d:p.second) sum += d;
		sum /= (double)(p.second.size());
		ea_.elocalbbhb.insert({p.first,sum});
		if(sum>cut0[ENE::LOCALBBHB]) he_.elocalbbhb.insert({p.first,sum});
	}
}

void TrajectoryEnergy::readcut(std::string parfile) {
	cut0.assign(ENE::ETERM,BigValue);
	cut1.assign(ENE::ETERM,BigValue);
	cut2.assign(ENE::ETERM,BigValue);
	std::ifstream ifs(parfile);
	std::string readline;
	std::stringstream ss;
	while(std::getline(ifs,readline)) {
		if(readline.empty()) continue;
		if(readline[0]=='#') continue;
		std::string lb;
		double d0, d1, d2;
		ss << readline;
		ss >> lb >> d0 >> d1 >> d2;
		ss.clear();
		int ix;
		if(lb=="bond") ix = ENE::BOND;
		else if(lb=="ang") ix = ENE::ANG;
		else if(lb=="impdih") ix = ENE::IMPDIH;
		else if(lb=="steric") ix = ENE::STERIC;
		else if(lb=="phipsi") ix = ENE::PHIPSI;
		else if(lb=="local") ix = ENE::LOCAL;
		else if(lb=="scconf") ix = ENE::SCCONF;
		else if(lb=="scpacking") ix = ENE::SCPACKING;
		else if(lb=="sitepair") ix = ENE::SITEPAIR;
		else if(lb=="localbbhb") ix = ENE::LOCALBBHB;
		else if(lb=="total") ix = ENE::TOTAL;
		else {
			std::cout <<"Can Not Recognize Energy Term :  " <<lb <<std::endl;
		}
		cut0[ix] = d0;
		cut1[ix] = d1;
		cut2[ix] = d2;
	}
	ifs.close();
}

std::map<std::string,double> TrajectoryEnergy::ene2rank() {
	std::map<std::string,double> ers;
	std::vector<std::string> ets = EnergyFilter::eneterm;
	for(std::string &s:ets) if(s!="s1" && s!="s2") ers.insert({s,0.0});

	for(auto &e:ea_.ebond) {
		ers.at("bond") += e.second;
	}
	for(auto &e:ea_.eang) {
		ers.at("ang") += e.second;
	}
	for(auto &e:ea_.eimpdih) {
		ers.at("impdih") += e.second;
	}
	for(auto &e:ea_.esteric) {
		ers.at("steric") += e.second;
	}
	for(auto &e:ea_.ephipsi) {
		ers.at("phipsi") += e.second;
	}
	for(auto &e:ea_.elocal) {
		ers.at("local") += e.second;
	}
	for(auto &e:ea_.etot) {
		ers.at("total") += e.second;
	}
	for(auto &e:ea_.escconf) {
		ers.at("scconf") += e.second;
	}
	for(auto &e:ea_.esitepair) {
		ers.at("sitepair") += e.second;
	}
	for(auto &e:ea_.escpacking) {
		ers.at("scpacking") += e.second;
	}
	for(auto &e:ea_.elocalbbhb) {
		ers.at("localbbhb") += e.second;
	}
	return ers;
}

/*
std::map<std::string,double> TrajectoryEnergy::ene2rank() {
	std::map<std::string,double> ers;
	std::vector<std::string> ets = EnergyFilter::eneterm;
	for(std::string &s:ets) if(s!="s1" && s!="s2") ers.insert({s,0.0});
	int nx = ENE::BOND;
	int ix = 0;
	double dx = 0.0;
	for(auto &e:em_.ebond) {
		for(auto &d:e.second) {
			dx += index(cut1[nx],cut2[nx],d);
		}
		ix += e.second.size();
	}
	if(ix!=0) ers.at("bond") = dx / (double)ix;
	nx = ENE::ANG;
	ix = 0;
	dx = 0.0;
	for(auto &e:em_.eang) {
		for(auto &d:e.second) {
			dx += index(cut1[nx],cut2[nx],d);
		}
		ix += e.second.size();
	}
	if(ix!=0) ers.at("ang") = dx / (double)ix;
	nx = ENE::IMPDIH;
	ix = 0;
	dx = 0.0;
	for(auto &e:em_.eimpdih) {
		for(auto &d:e.second) {
			dx += index(cut1[nx],cut2[nx],d);
		}
		ix += e.second.size();
	}
	if(ix!=0) ers.at("impdih") = dx / (double)ix;
	nx = ENE::STERIC;
	ix = 0;
	dx = 0.0;
	for(auto &e:em_.esteric) {
		for(auto &d:e.second) {
			dx += index(cut1[nx],cut2[nx],d);
		}
		ix += e.second.size();
	}
	if(ix!=0) ers.at("steric") = dx / (double)ix;
	nx = ENE::PHIPSI;
	ix = 0;
	dx = 0.0;
	for(auto &e:em_.ephipsi) {
		for(auto &d:e.second) {
			dx += index(cut1[nx],cut2[nx],d);
		}
		ix += e.second.size();
	}
	if(ix!=0) ers.at("phipsi") = dx / (double)ix;
	nx = ENE::LOCAL;
	ix = 0;
	dx = 0.0;
	for(auto &e:em_.elocal) {
		for(auto &d:e.second) {
			dx += index(cut1[nx],cut2[nx],d);
		}
		ix += e.second.size();
	}
	if(ix!=0) ers.at("local") = dx / (double)ix;
	nx = ENE::SCCONF;
	ix = 0;
	dx = 0.0;
	for(auto &e:em_.escconf) {
		for(auto &d:e.second) {
			dx += index(cut1[nx],cut2[nx],d);
		}
		ix += e.second.size();
	}
	if(ix!=0) ers.at("scconf") = dx / (double)ix;
	nx = ENE::SITEPAIR;
	ix = 0;
	dx = 0.0;
	for(auto &e:em_.esitepair) {
		for(auto &d:e.second) {
			dx += index(cut1[nx],cut2[nx],d);
		}
		ix += e.second.size();
	}
	if(ix!=0) ers.at("sitepair") = dx / (double)ix;
	nx = ENE::LOCALBBHB;
	ix = 0;
	dx = 0.0;
	for(auto &e:em_.elocalbbhb) {
		for(auto &d:e.second) {
			dx += index(cut1[nx],cut2[nx],d);
		}
		ix += e.second.size();
	}
	if(ix!=0) ers.at("localbbhb") = dx / (double)ix;
	nx = ENE::SCPACKING;
	ix = 0;
	dx = 0.0;
	for(auto &e:em_.escpacking) {
		for(auto &d:e.second) {
			dx += index(cut1[nx],cut2[nx],d);
		}
		ix += e.second.size();
	}
	if(ix!=0) ers.at("scpacking") = dx / (double)ix;
	nx = ENE::TOTAL;
	ix = 0;
	dx = 0.0;
	for(auto &e:em_.etot) {
		for(auto &d:e.second) {
			dx += index(cut1[nx],cut2[nx],d);
		}
		ix += e.second.size();
	}
	if(ix!=0) ers.at("total") = dx / (double)ix;
	return ers;
}*/

double TrajectoryEnergy::index(double c1, double c2, double v) {
	// default :   c1 < c2
	double delta = c2 -c1;
	if(fabs(delta)<0.00001) return 0.0;
	if(v<c1) return 0.0;
	if(v>c2) return 1.0;
	double mid = (c1+c2) /2.0;
	double a = 2.0 / delta / delta;
	if(v<mid) {
		double dis = v - c1;
		return a * dis * dis;
	} else {
		double dis = c2 - v;
		return 1.0 - a * dis * dis;
	}
}

double TrajectoryEnergy::index01(double c1, double c2, double v) {
	// default :   c1 < c2
	double delta = c2 -c1;
	if(fabs(delta)<0.00001) return 0.0;
	if(v<c1) return 0.0;
	if(v>c2) return 1.0;
	double mid = (c1+c2) /2.0;
	double a = 2.0 / delta / delta;
	if(v<mid) {
		double dis = v - c1;
		return a * dis * dis;
	} else {
		double dis = c2 - v;
		return 1.0 - a * dis * dis;
	}
}

void TrajectoryEnergy::EnergyMerge::standard(const std::vector<double> &cut11, const std::vector<double> &cut22) {
	for(auto &p:ephipsi) {
		int ix = ENE::PHIPSI;
		for(auto &d:p.second) d = TrajectoryEnergy::index01(cut11[ix],cut22[ix],d);
	}
	for(auto &p:elocal) {
		int ix = ENE::LOCAL;
		for(auto &d:p.second) d = TrajectoryEnergy::index01(cut11[ix],cut22[ix],d);
	}
	for(auto &p:esitepair) {
		int ix = ENE::SITEPAIR;
		for(auto &d:p.second) d = TrajectoryEnergy::index01(cut11[ix],cut22[ix],d);
	}
	for(auto &p:elocalbbhb) {
		int ix = ENE::LOCALBBHB;
		for(auto &d:p.second) d = TrajectoryEnergy::index01(cut11[ix],cut22[ix],d);
	}
}


std::vector<std::vector<int>> EnergyTermRecord::v1_v2() {
	std::vector<std::vector<int>> v12;
	for(int i=0;i<decoy.size();i++) {
		for(int j=0;j<decoy[i].size();j++) {
			v12.push_back({i,j});
		}
	}
	return v12;
}

void EnergyTermRecord::read_s1_rot(std::string enefile, std::map<SingleResidue,double> &es) {
	std::vector<std::vector<int>> v12 = v1_v2();
	std::ifstream ifs(enefile);
	std::string readline;
	int nc=0;
	while(std::getline(ifs,readline)) {
		double e = std::stod(readline);
		es.insert({v12[nc++],e});
	}
	ifs.close();
}

void EnergyTermRecord::read_s2_pck(std::string enefile, std::map<MultiResidue,double> &es) {
	std::vector<std::vector<int>> v12 = v1_v2();
	std::ifstream ifs(enefile);
	std::string readline;
	int nc=0;
	while(std::getline(ifs,readline)) {
		double i = std::stoi(readline);
		std::getline(ifs,readline);
		double j = std::stoi(readline);
		std::getline(ifs,readline);
		double e = std::stod(readline);
		es.insert({{v12[i],v12[j]},e});
	}
	ifs.close();
}

void EnergyTermRecord::readss(std::string ssfile) {
	std::ifstream ifs(ssfile);
	std::string readline;
	int nc=-1;
	while(std::getline(ifs,readline)) {
		nc++;
		for(int i=0;i<readline.size();i++) {
			sss.insert({{nc,i},readline[i]});
		}
	}
	ifs.close();
}














/*

void Filter::totalene_eterm() {
	etot_.assign(ENE::SITEPAIR+1,0.0);
	for(auto &es:enes_res) {
		for(auto &e:es) {
			if(fabs(e.ephipsi-BigValue)>0.00001) etot_[ENE::PHIPSI] += e.ephipsi;
			if(fabs(e.elocal-BigValue)>0.00001) etot_[ENE::LOCAL] += e.elocal;
			if(fabs(e.escconf-BigValue)>0.00001) etot_[ENE::SCCONF] += e.escconf;
			for(auto &p:e.esitepair) etot_[ENE::SITEPAIR] += p.second;
			for(auto &p:e.escpacking) etot_[ENE::SCPACKING] += p.second;
			for(auto &p:e.elocalbbhb) etot_[ENE::LOCALBBHB] += p.second;
		}
	}
	for(auto &es:enes_atom) {
		for(auto &e1:es) {
			for(auto &e2:e1) {
				Filter::Ene_Atom & e=e2.second;
				for(auto &p:e.ebond) etot_[ENE::BOND] += p.second;
				for(auto &p:e.eang) etot_[ENE::ANG] += p.second;
				for(auto &p:e.eimpdih) etot_[ENE::IMPDIH] += p.second;
				for(auto &p:e.esteric) etot_[ENE::STERIC] += p.second;
			}
		}
	}
	etot_[ENE::BOND] /= 2.0;
	etot_[ENE::ANG] /= 3.0;
	etot_[ENE::IMPDIH] /= 4.0;
	etot_[ENE::STERIC] /= 2.0;
	etot_[ENE::SCCONF] /= 1.0;
	etot_[ENE::PHIPSI] /= 1.0;
	etot_[ENE::LOCAL] /= 1.0;
	etot_[ENE::SCPACKING] /= 2.0;
	etot_[ENE::LOCALBBHB] /= 2.0;
	etot_[ENE::SITEPAIR] /= 2.0;
	for(int i=0;i<10;i++) {
		etot_[ENE::TOTAL] += etot_[i+1];
	}
}

bool Filter::equalatoms(const std::vector<atominres> &air1, const std::vector<atominres> &air2) {
	if(air1.size()!=air2.size()) return false;
	for(int i=0;i<air1.size();i++) {
		bool finded{false};
		for(int j=0;j<air2.size();j++) {
			if(air1[i]!=air2[j]) continue;
			finded=true;
			break;
		}
		if(!finded) return false;
	}
	return true;
}

void Filter::island_space(const std::vector<std::vector<NSPgeometry::XYZ>> &crds,
		double discut, std::vector<std::vector<int>> &ild) {
	if(crds.empty()) return;
	ild.push_back({0});
	double dc2=discut*discut;
	for(int i=1;i<crds.size();i++) {
		std::vector<int> rs;
		std::set<int> rs1;
		for(int j=0;j<crds[i].size();j++) {
			for(int m=0;m<ild.size();m++) {
				if(rs1.find(m)!=rs1.end()) continue;
				for(int n=0;n<ild[m].size();n++) {
					for(int a=0;a<crds[ild[m][n]].size();a++) {
						NSPgeometry::XYZ c=crds[ild[m][n]][a];
						if((c-crds[i][j]).squarednorm()>dc2) continue;
						rs.push_back(m);
						rs1.insert(m);
						break;
					}
					if(rs1.find(m)!=rs1.end()) break;
				}
			}
			bool isall{true};
			for(int k=0;k<ild.size();k++) {
				if(rs1.find(k)!=rs1.end()) continue;
				isall=false;
				break;
			}
			if(isall) break;
		}
		if(rs.empty()) {
			ild.push_back({i});
		} else {
			ild[rs[0]].push_back(i);
			for(int j=1;j<rs.size();j++) {
				for(int k:ild[rs[j]]) ild[rs[0]].push_back(k);
				ild[rs[j]].clear();
			}
			std::vector<std::vector<int>> ild1;
			for(int j=0;j<ild.size();j++) {
				if(ild[j].empty()) continue;
				ild1.push_back(ild[j]);
			}
			ild=ild1;
		}
	}
}

void Filter::island_sequence(const std::vector<std::vector<int>>&hrs, int discut, std::vector<std::vector<int>>&ild) {
	if(hrs.empty()) return;
	ild.push_back({0});
	for(int i=1;i<hrs.size();i++) {
		std::vector<int> rs;
		std::set<int> rs1;
		for(int j=0;j<hrs[i].size();j+=2) {
			int cid0=hrs[i][j];
			int rid0=hrs[i][j+1];
			for(int m=0;m<ild.size();m++) {
				if(rs1.find(m)!=rs1.end()) continue;
				for(int n=0;n<ild[m].size();n++) {
					for(int a=0;a<hrs[ild[m][n]].size();a+=2) {
						int cid1=hrs[ild[m][n]][a];
						int rid1=hrs[ild[m][n]][a+1];
						if(cid0!=cid1) continue;
						if(rid0-rid1<-discut || rid0-rid1>discut) continue;
						rs.push_back(m);
						rs1.insert(m);
						break;
					}
					if(rs1.find(m)!=rs1.end()) break;
				}
			}
			bool isall{true};
			for(int k=0;k<ild.size();k++) {
				if(rs1.find(k)!=rs1.end()) continue;
				isall=false;
				break;
			}
			if(isall) break;
		}
		if(rs.empty()) {
			ild.push_back({i});
		} else {
			ild[rs[0]].push_back(i);
			for(int j=1;j<rs.size();j++) {
				for(int k:ild[rs[j]]) ild[rs[0]].push_back(k);
				ild[rs[j]].clear();
			}
			std::vector<std::vector<int>> ild1;
			for(int j=0;j<ild.size();j++) {
				if(ild[j].empty()) continue;
				ild1.push_back(ild[j]);
			}
			ild=ild1;
		}
	}
}*/

/*void Filter::island_atom(std::string label) {
	double mincutoff = (label=="bond"?highregion_.pbond.first:(
			label=="ang"?highregion_.pang.first:(
					label=="impdih"?highregion_.pimpdih.first:highregion_.psteric.first)));
	std::vector<std::vector<atominres>> & highene = (
			label=="bond"?highregion_.ubond:(
					label=="ang"?highregion_.uang:(
							label=="impdih"?highregion_.uimpdih:highregion_.usteric)));
	int & ne= (label=="bond"?highregion_.nbond:(
			label=="ang"?highregion_.nang:(
					label=="impdih"?highregion_.nimpdih:highregion_.nsteric)));
	for(auto &p1:enes_atom) {
		for(auto &p2:p1) {
			for(auto &p3:p2) {
				std::vector<std::pair<std::vector<atominres>,double>> & es=(
						label=="bond"?p3.second.ebond:(
								label=="ang"?p3.second.eang:(
										label=="impdih"?p3.second.eimpdih:p3.second.esteric)));
				for(auto &p4:es) {
					ne++;
					if(p4.second<mincutoff) continue;
					std::vector<atominres> va=p4.first;
					va.push_back({{p3.second.chainid,p3.second.resid},p3.second.atomname});
					bool finded{false};
					for(int i=0;i<highene.size();i++) {
						if(!equalatoms(va,highene[i])) continue;
						finded=true;
						break;
					}
					if(!finded) highene.push_back(va);
				}
			}
		}
	}
	int nn= (label=="ang"?3:(label=="impdih"?4:2));
	ne /= nn;

	if(highene.empty()) return;
	std::vector<std::vector<int>> & ild = (label=="bond"?highregion_.ibond:(
			label=="ang"?highregion_.iang:(
					label=="impdih"?highregion_.iimpdih:highregion_.isteric)));
	double discut = (label=="bond"?highregion_.pbond.second:(
			label=="ang"?highregion_.pang.second:(
					label=="impdih"?highregion_.pimpdih.second:highregion_.psteric.second)));
	std::vector<std::vector<NSPgeometry::XYZ>> ps;
	for(int i=0;i<highene.size();i++) {
		ps.push_back(std::vector<NSPgeometry::XYZ>());
		for(int j=0;j<highene[i].size();j++) {
			ps.back().push_back(chs_[highene[i][j].first[0]][highene[i][j].first[1]].rds().at(highene[i][j].second).crd);
		}
	}
	island_space(ps,discut,ild);*/
	/*
	std::vector<NSPgeometry::XYZ> cens;
	for(int i=0;i<highene.size();i++) {
		int sz=high[0].size();
		NSPgeometry::XYZ o(0.0,0.0,0.0);
		for(int j=0;j<sz;j++) {
			o = o + chs_[highene[i][j].first[0]][highene[i][j].first[1]].rds().at([highene[i][j].second).crd;
		}
		o=o/(double)sz;
		cens.push_back(o);
	}
	double dc2=discut*discut;
	if(cens.empty()) return;
	ild.push_back({0});
	for(int i=1;i<cens.size();i++) {
		std::vector<int> rs;
		for(int j=0;j<ild.size();j++) {
			for(int k=0;k<ild[j].size();k++) {
				NSPgeometry::XYZ c=cens[ild[j][k]];
				if((c-cens[i]).squarednorm()>dc2) continue;
				rs.push_back(j);
				break;
			}
		}
		if(rs.empty()) {
			ild.push_back({i});
		} else {
			ild[rs[0]].push_back(i);
			for(int j=1;j<rs.size();j++) {
				for(int k:ild[rs[j]]) ild[rs[0]].push_back(k);
				ild[rs[j]].clear();
			}
			std::vector<std::vector<int>> ild1;
			for(int j=0;j<ild.size();j++) {
				if(ild[j].empty()) continue;
				ild1.push_back(ild[j]);
			}
			ild=ild1;
		}
	}*/
//}
/*
void Filter::island_respair(std::string label) {
	double mincutoff = (label=="sitepair"?highregion_.psitepair.first:(
			label=="scpacking"?highregion_.pscpacking.first:highregion_.plocalbbhb.first));
	std::vector<std::vector<int>> & highene = (label=="sitepair"?highregion_.usitepair:(
			label=="scpacking"?highregion_.uscpacking:highregion_.ulocalbbhb));
	int & ne= (label=="sitepair"?highregion_.nsitepair:(
			label=="scpacking"?highregion_.nscpacking:highregion_.nlocalbbhb));
	for(auto &p1:enes_res) {
		for(Ene_Res &p2:p1) {
			std::map<std::pair<int,int>,double> & es = (label=="sitepair"?p2.esitepair:(
					label=="scpacking"?p2.escpacking:p2.elocalbbhb));
			for(auto &p3:es) {
				ne++;
				if(p3.second<mincutoff) continue;
				std::vector<int> airs{p2.chainid,p2.resid,p3.first.first,p3.first.second};
				bool finded{false};
				for(int i=0;i<highene.size();i++) {
					if(airs[0]==highene[i][0] && airs[1]==highene[i][1] && airs[2]==highene[i][2] && airs[3]==highene[i][3])
						finded=true;
					if(airs[0]==highene[i][2] && airs[1]==highene[i][3] && airs[2]==highene[i][0] && airs[3]==highene[i][1])
						finded=true;
					if(finded) break;
				}
				if(!finded) highene.push_back(airs);
			}
		}
	}
	ne /= 2;

	if(highene.empty()) return;
	std::vector<std::vector<int>> & ild = (label=="localbbhb"?highregion_.ilocalbbhb:(
			label=="scpacking"?highregion_.iscpacking:highregion_.isitepair));
	double discut=(label=="localbbhb"?highregion_.plocalbbhb.second:(
			label=="scpacking"?highregion_.pscpacking.second:highregion_.psitepair.second));
	if(label=="localbbhb") {
		island_sequence(highene,discut,ild);
	} else {
		std::vector<std::vector<NSPgeometry::XYZ>> cs;
		std::set<std::string> mcs{"N","CA","C","O"};
		for(int i=0;i<highene.size();i++) {
			cs.push_back(std::vector<NSPgeometry::XYZ>());
			for(int j=0;j<highene[i].size();j+=2) {
				int cid=highene[i][j];
				int rid=highene[i][j+1];
				if(label=="scpacking") {
					auto & cs1=chs_[cid][rid].rds();
					for(auto &p:cs1) {
						if(mcs.find(p.first)!=mcs.end()) continue;
						cs.back().push_back(p.second.crd);
					}
				} else {//sitepair
					for(auto &s:mcs) {
						cs.back().push_back(chs_[cid][rid].rds().at(s).crd);
					}
				}
			}
		}
		island_space(cs,discut,ild);
	}
}*/
/*
void Filter::island_single(std::string label) {
	double mincutoff = (label=="scconf"?highregion_.pscconf.first:(
			label=="phipsi"?highregion_.pphipsi.first:highregion_.plocal.first));
	std::vector<std::vector<int>> & highene = (label=="scconf"?highregion_.uscconf:(
			label=="phipsi"?highregion_.uphipsi:highregion_.ulocal));
	int & ne= (label=="scconf"?highregion_.nscconf:(
			label=="phipsi"?highregion_.nphipsi:highregion_.nlocal));
	for(auto &p1:enes_res) {
		for(auto &p2:p1) {
			double ee=(label=="scconf"?p2.escconf:(label=="phipsi"?p2.ephipsi:p2.elocal));
			if(fabs(ee-BigValue)<0.0001) continue;
			ne++;
			if(mincutoff>ee) continue;
			highene.push_back({p2.chainid,p2.resid});
		}
	}
	std::vector<std::vector<int>> & ild = (label=="scconf"?highregion_.iscconf:(
			label=="phipsi"?highregion_.iphipsi:highregion_.ilocal));
	int discut = (label=="scconf"?highregion_.pscconf.second:(
			label=="phipsi"?highregion_.pphipsi.second:highregion_.plocal.second));
	if(highene.empty()) return;
	island_sequence(highene,discut,ild);*/
	/*
	ild.push_back({0});
	for(int i=1;i<highene.size();i++) {
		int cid0=highene[i][0];
		int rid0=highene[i][1];
		std::vector<int> rs;
		for(int j=0;j<ild.size();j++) {
			for(int k=0;k<ild[j].size();k++) {
				int cid1=highene[ild[j][k]][0];
				int rid1=highene[ild[j][k]][1];
				if(cid0!=cid1) continue;
				if(cid0-cid1<-discut || cid0-cid1>discut) continue;
				rs.push_back(j);
				break;
			}
		}
		if(rs.empty()) {
			ild.push_back({i});
		} else {
			ild[rs[0]].push_back(i);
			for(int j=1;j<rs.size();j++) {
				for(int k:ild[rs[j]]) ild[rs[0]].push_back(k);
				ild[rs[j]].clear();
			}
			std::vector<std::vector<int>> ild1;
			for(int j=0;j<ild.size();j++) {
				if(ild[j].empty()) continue;
				ild1.push_back(ild[j]);
			}
			ild=ild1;
		}
	}*/
//}
/*
void Filter::island_bond(double e, double d) {
	highregion_.ubond.clear();
	highregion_.ibond.clear();
	highregion_.pbond={e,d};
	island_atom(eneterm[ENE::BOND]);
}

void Filter::island_ang(double e, double d) {
	highregion_.uang.clear();
	highregion_.iang.clear();
	highregion_.pang={e,d};
	island_atom(eneterm[ENE::ANG]);
}

void Filter::island_impdih(double e, double d) {
	highregion_.uimpdih.clear();
	highregion_.iimpdih.clear();
	highregion_.pimpdih={e,d};
	island_atom(eneterm[ENE::IMPDIH]);
}

void Filter::island_steric(double e, double d) {
	highregion_.usteric.clear();
	highregion_.isteric.clear();
	highregion_.psteric={e,d};
	island_atom(eneterm[ENE::STERIC]);
}

void Filter::island_scconf(double e, int d) {
	highregion_.uscconf.clear();
	highregion_.iscconf.clear();
	highregion_.pscconf={e,d};
	island_single(eneterm[ENE::SCCONF]);
}

void Filter::island_phipsi(double e, int d) {
	highregion_.uphipsi.clear();
	highregion_.iphipsi.clear();
	highregion_.pphipsi={e,d};
	island_single(eneterm[ENE::PHIPSI]);
}

void Filter::island_local(double e, int d) {
	highregion_.ulocal.clear();
	highregion_.ilocal.clear();
	highregion_.plocal={e,d};
	island_single(eneterm[ENE::LOCAL]);
}

void Filter::island_localbbhb(double e, int d) {
	highregion_.ulocalbbhb.clear();
	highregion_.ilocalbbhb.clear();
	highregion_.plocalbbhb={e,d};
	island_respair(eneterm[ENE::LOCALBBHB]);
}

void Filter::island_scpacking(double e, double d) {
	highregion_.uscpacking.clear();
	highregion_.iscpacking.clear();
	highregion_.pscpacking={e,d};
	island_respair(eneterm[ENE::SCPACKING]);
}

void Filter::island_sitepair(double e, double d) {
	highregion_.usitepair.clear();
	highregion_.isitepair.clear();
	highregion_.psitepair={e,d};
	island_respair(eneterm[ENE::SITEPAIR]);
}*/
/*
void Filter::eneisland(std::map<std::string,std::pair<double,double>> &lb1, std::map<std::string,std::pair<double,int>> &lb2) {
	island_bond(lb1.at(eneterm[ENE::BOND]).first, lb1.at(eneterm[ENE::BOND]).second);
	island_ang(lb1.at(eneterm[ENE::ANG]).first, lb1.at(eneterm[ENE::ANG]).second);
	island_impdih(lb1.at(eneterm[ENE::IMPDIH]).first, lb1.at(eneterm[ENE::IMPDIH]).second);
	island_steric(lb1.at(eneterm[ENE::STERIC]).first, lb1.at(eneterm[ENE::STERIC]).second);
	island_scconf(lb2.at(eneterm[ENE::SCCONF]).first, lb2.at(eneterm[ENE::SCCONF]).second);
	island_phipsi(lb2.at(eneterm[ENE::PHIPSI]).first, lb2.at(eneterm[ENE::PHIPSI]).second);
	island_local(lb2.at(eneterm[ENE::LOCAL]).first, lb2.at(eneterm[ENE::LOCAL]).second);
	island_localbbhb(lb2.at(eneterm[ENE::LOCALBBHB]).first, lb2.at(eneterm[ENE::LOCALBBHB]).second);
	island_scpacking(lb1.at(eneterm[ENE::SCPACKING]).first, lb1.at(eneterm[ENE::SCPACKING]).second);
	island_sitepair(lb1.at(eneterm[ENE::SITEPAIR]).first, lb1.at(eneterm[ENE::SITEPAIR]).second);*/

	/*
	eneunitclear();
	island_atom("bond");
	island_atom("ang");
	island_atom("impdih");
	island_atom("steric");
	island_respair("scpacking");
	island_respair("localbbhb");
	island_respair("sitepair");
	island_single("scconf");
	island_single("phipsi");
	island_single("local");*/
//}
/*
void Filter::eneisland1() {
	std::map<std::string,std::pair<double,double>> lb1{{"bond",highregion_.pbond},{"ang",highregion_.pang},
			{"impdih",highregion_.pimpdih},{"steric",highregion_.psteric},{"scpacking",highregion_.pscpacking},
			{"sitepair",highregion_.psitepair}};
	std::map<std::string,std::pair<double,int>> lb2{{"scconf",highregion_.pscconf},{"phipsi",highregion_.pphipsi},
			{"local",highregion_.plocal},{"localbbhb",highregion_.plocalbbhb}};
	eneisland(lb1,lb2);
}*/

/*
void Filter::eneunitclear() {
	highregion_.ubond.clear(); highregion_.ibond.clear();
	highregion_.uang.clear(); highregion_.iang.clear();
	highregion_.uimpdih.clear(); highregion_.iimpdih.clear();
	highregion_.usteric.clear(); highregion_.isteric.clear();
	highregion_.uscconf.clear(); highregion_.iscconf.clear();
	highregion_.uphipsi.clear(); highregion_.iphipsi.clear();
	highregion_.ulocal.clear(); highregion_.ilocal.clear();
	highregion_.ulocalbbhb.clear(); highregion_.ilocalbbhb.clear();
	highregion_.uscpacking.clear(); highregion_.iscpacking.clear();
	highregion_.usitepair.clear(); highregion_.isitepair.clear();
}
*/
/*
Filter::EneInKind Filter::geteneinkind() {
	EneInKind eik;
	for(auto &es:enes_res) {
		for(Ene_Res &er:es) {
			std::string resname=chs_[er.chainid][er.resid].resname();
			if(fabs(er.escconf-BigValue)>0.0001) eik.escconf[resname].push_back(er.escconf);
			if(fabs(er.ephipsi-BigValue)>0.0001) eik.ephipsi[resname].push_back(er.ephipsi);
			if(fabs(er.elocal-BigValue)>0.0001) eik.elocal[resname].push_back(er.elocal);
			for(auto &em:er.esitepair) {
				if(em.first.first<er.chainid || (em.first.first==er.chainid && em.first.second<=er.resid)) continue;
				std::string rn1=chs_[em.first.first][em.first.second].resname();
				std::pair<std::string,std::string> p{resname,rn1};
				if(resname>rn1) p={rn1,resname};
				eik.esitepair[p].push_back(em.second);
			}
			for(auto &em:er.escpacking) {
				if(em.first.first<er.chainid || (em.first.first==er.chainid && em.first.second<=er.resid)) continue;
				std::string rn1=chs_[em.first.first][em.first.second].resname();
				std::pair<std::string,std::string> p{resname,rn1};
				if(resname>rn1) p={rn1,resname};
				eik.escpacking[p].push_back(em.second);
			}
			for(auto &em:er.elocalbbhb) {
				if(em.first.first<er.chainid || (em.first.first==er.chainid && em.first.second<=er.resid)) continue;
				std::string rn1=chs_[em.first.first][em.first.second].resname();
				std::pair<std::string,std::string> p{resname,rn1};
				if(resname>rn1) p={rn1,resname};
				eik.elocalbbhb[p].push_back(em.second);
			}
		}
	}
	for(auto &es:enes_atom) {
		for(auto &emp:es) {
			for(auto &ep:emp) {
				Ene_Atom & ea = ep.second;
				AtomName an0{chs_[ea.chainid][ea.resid].resname(),ea.atomname};
				for(auto &em:ea.ebond) {
					int cid=em.first[0].first[0], rid=em.first[0].first[1];
					std::string aid=em.first[0].second;
					if(cid<ea.chainid) continue;
					if(cid==ea.chainid && rid<ea.resid) continue;
					if(cid==ea.chainid && rid==ea.resid && aid<=ea.atomname) continue;
					AtomName an1{chs_[cid][rid].resname(),aid};
					std::pair<AtomName,AtomName> p{an0,an1};
					if(an1<an0) p={an1,an0};
					eik.ebond[p].push_back(em.second);
				}
				for(auto &em:ea.esteric) {
					int cid=em.first[0].first[0], rid=em.first[0].first[1];
					std::string aid=em.first[0].second;
					if(cid<ea.chainid) continue;
					if(cid==ea.chainid && rid<ea.resid) continue;
					if(cid==ea.chainid && rid==ea.resid && aid<=ea.atomname) continue;
					AtomName an1{chs_[cid][rid].resname(),aid};
					std::pair<AtomName,AtomName> p{an0,an1};
					if(an1<an0) p={an1,an0};
					eik.esteric[p].push_back(em.second);
				}
				for(auto &em:ea.eang) {
					int cid=em.first[0].first[0], rid=em.first[0].first[1];
					std::string aid=em.first[0].second;
					if(cid<ea.chainid) continue;
					if(cid==ea.chainid && rid<ea.resid) continue;
					if(cid==ea.chainid && rid==ea.resid && aid<=ea.atomname) continue;
					AtomName an1{chs_[cid][rid].resname(),aid};

					cid=em.first[1].first[0];
					rid=em.first[1].first[1];
					aid=em.first[1].second;
					if(cid<ea.chainid) continue;
					if(cid==ea.chainid && rid<ea.resid) continue;
					if(cid==ea.chainid && rid==ea.resid && aid<=ea.atomname) continue;
					AtomName an2{chs_[cid][rid].resname(),aid};

					std::vector<AtomName> ans{an0,an1,an2};
					std::sort(ans.begin(),ans.end());
					std::pair<AtomName,std::pair<AtomName,AtomName>> p{ans[0],{ans[1],ans[2]}};
					eik.eang[p].push_back(em.second);
				}
				for(auto &em:ea.eimpdih) {
					int cid=em.first[0].first[0], rid=em.first[0].first[1];
					std::string aid=em.first[0].second;
					if(cid<ea.chainid) continue;
					if(cid==ea.chainid && rid<ea.resid) continue;
					if(cid==ea.chainid && rid==ea.resid && aid<=ea.atomname) continue;
					AtomName an1{chs_[cid][rid].resname(),aid};

					cid=em.first[1].first[0];
					rid=em.first[1].first[1];
					aid=em.first[1].second;
					if(cid<ea.chainid) continue;
					if(cid==ea.chainid && rid<ea.resid) continue;
					if(cid==ea.chainid && rid==ea.resid && aid<=ea.atomname) continue;
					AtomName an2{chs_[cid][rid].resname(),aid};

					cid=em.first[2].first[0];
					rid=em.first[2].first[1];
					aid=em.first[1].second;
					if(cid<ea.chainid) continue;
					if(cid==ea.chainid && rid<ea.resid) continue;
					if(cid==ea.chainid && rid==ea.resid && aid<=ea.atomname) continue;
					AtomName an3{chs_[cid][rid].resname(),aid};

					std::vector<AtomName> ans{an0,an1,an2,an3};
					std::sort(ans.begin(),ans.end());
					std::pair<std::pair<AtomName,AtomName>,std::pair<AtomName,AtomName>> p{{ans[0],ans[1]},{ans[2],ans[3]}};
					eik.eimpdih[p].push_back(em.second);
				}
			}
		}
	}
	return eik;
}

void Filter::geteneunitpar() {
	std::string parfile=NSPdataio::datafilename("IslandCutOffMC95");
	std::set<std::string> integers{Filter::eneterm[Filter::ENE::SCCONF],
		Filter::eneterm[Filter::ENE::PHIPSI],Filter::eneterm[Filter::ENE::LOCAL],
		Filter::eneterm[Filter::ENE::LOCALBBHB]};
	std::string readline;
	std::stringstream ss;
	std::ifstream ifs(parfile);
	while(std::getline(ifs,readline)) {
		std::string str;
		int ix{-1};
		double dx{-1.0};
		double e;
		int sz;
		ss << readline;
		ss >> str;
		if(integers.find(readline)==integers.end()) ss >>dx;
		else ss >>ix;
		ss >>e >>sz;
		ss.clear();
		if(str=="bond") {
			highregion_.pbond = {e,dx};
			highregion_.szbond = sz;
		}
		if(str=="ang") {
			highregion_.pang = {e,dx};
			highregion_.szang = sz;
		}
		if(str=="impdih") {
			highregion_.pimpdih = {e,dx};
			highregion_.szimpdih = sz;
		}
		if(str=="scconf") {
			highregion_.pscconf = {e,dx};
			highregion_.szscconf = sz;
		}
		if(str=="scpacking") {
			highregion_.pscpacking = {e,dx};
			highregion_.szscpacking = sz;
		}
		if(str=="steric") {
			highregion_.psteric = {e,dx};
			highregion_.szsteric = sz;
		}
		if(str=="local") {
			highregion_.plocal = {e,dx};
			highregion_.szlocal = sz;
		}
		if(str=="phipsi") {
			highregion_.pphipsi = {e,dx};
			highregion_.szphipsi = sz;
		}
		if(str=="localbbhb") {
			highregion_.plocalbbhb = {e,dx};
			highregion_.szlocalbbhb = sz;
		}
		if(str=="sitepair") {
			highregion_.psitepair = {e,dx};
			highregion_.szsitepair = sz;
		}
	}
	ifs.close();
}

std::string Filter::heainlocalregion() {
	for(auto &vi:highregion_.ibond) {
		if(vi.size()>highregion_.szbond) return "bond";
	}
	for(auto &vi:highregion_.iang) {
		if(vi.size()>highregion_.szang) return "ang";
	}
	for(auto &vi:highregion_.iimpdih) {
		if(vi.size()>highregion_.szimpdih) return "impdih";
	}
	for(auto &vi:highregion_.isteric) {
		if(vi.size()>highregion_.szsteric) return "steric";
	}
	for(auto &vi:highregion_.iphipsi) {
		if(vi.size()>highregion_.szphipsi) return "phipsi";
	}
	for(auto &vi:highregion_.ilocal) {
		if(vi.size()>highregion_.szlocal) return "local";
	}
	for(auto &vi:highregion_.ilocalbbhb) {
		if(vi.size()>highregion_.szlocalbbhb) return "localbbhb";
	}
	for(auto &vi:highregion_.isitepair) {
		if(vi.size()>highregion_.szsitepair) return "sitepair";
	}
	return "nohea";
}*/
/*
std::string Filter::heainlocalregion() {
	for(auto &vi:highregion_.ibond) {
		if(vi.size()>highregion_.szbond) return "bond";
	}
	for(auto &vi:highregion_.iang) {
		if(vi.size()>highregion_.szang) return "ang";
	}
	for(auto &vi:highregion_.iimpdih) {
		if(vi.size()>highregion_.szimpdih) return "impdih";
	}
	for(auto &vi:highregion_.isteric) {
		if(vi.size()>highregion_.szsteric) return "steric";
	}
	for(auto &vi:highregion_.iphipsi) {
		if(vi.size()>highregion_.szphipsi) return "phipsi";
	}
	for(auto &vi:highregion_.ilocal) {
		if(vi.size()>highregion_.szlocal) return "local";
	}
	for(auto &vi:highregion_.ilocalbbhb) {
		if(vi.size()>highregion_.szlocalbbhb) return "localbbhb";
	}
	for(auto &vi:highregion_.isitepair) {
		if(vi.size()>highregion_.szsitepair) return "sitepair";
	}
	return "nohea";
}
*/









/*
int Filter::exposed(std::vector<NSPgeometry::XYZ>&cs, double rad) {
	int ne=0;
	double r2=rad*rad;
	for(auto &c:cs) {
		bool exp{true};
		for(auto &ch:chs_) {
			for(auto &res:ch) {
				for(auto &rd:res.rds()) {
					NSPgeometry::XYZ r=rd.second.crd;
					if((c-r).squarednorm()<r2) {
						exp=false;
						ne++;
						break;
					}
				}
				if(!exp) break;
			}
			if(!exp) break;
		}
	}
	return ne;
}

void Filter::findinnerpolar(int pmax, double rad, std::vector<std::pair<std::vector<int>,std::string>>&pas) {
	std::vector<NSPgeometry::XYZ> aas=AA_Atom::sphere256();
	for(auto &a:aas) a=a*rad;
	for(auto &hbs:hbps_) {
		for(auto &hs:hbs) {
			for(auto &h:hs.hbs) {
				if(!h.second.empty()) continue;
				std::vector<NSPgeometry::XYZ> cs=aas;
				NSPgeometry::XYZ cen=chs_[hs.chainid][hs.resid].rds().at(h.first).crd;
				for(auto &c:cs) c=c+cen;
				if(exposed(cs,rad)>pmax) pas.push_back({{hs.chainid,hs.resid},h.first});
			}
		}
	}
}

std::vector<double> Filter::gethighene() {
	std::vector<double> es(6,0.0); // tot, steric, phipsi, local, localbbhb, sitepair
	for(int i=0;i<chs_.size();i++) {
		for(int j=0;j<chs_[i].size();j++) {
			for(auto &pr:enes_atom[i][j]) {
				std::vector<std::pair<std::vector<atominres>,double>> &est=pr.second.esteric;
				for(auto &p:est) {
					if(p.second<highregion_.psteric.first) continue;
					es[1] += p.second;
				}
			}
			if(enes_res[i][j].ephipsi>highregion_.pphipsi.first) {
				es[2] += enes_res[i][j].ephipsi;
			}
			if(enes_res[i][j].elocal>highregion_.plocal.first) {
				es[3] += enes_res[i][j].elocal;
			}

			std::map<std::pair<int,int>,double> &elo=enes_res[i][j].elocalbbhb;
			for(auto &p:elo) {
				if(p.second<highregion_.plocalbbhb.first) continue;
				es[4] += p.second;
			}
			std::map<std::pair<int,int>,double> &esi=enes_res[i][j].esitepair;
			for(auto &p:esi) {
				if(p.second<highregion_.psitepair.first) continue;
				es[5] += p.second;
			}
		}
	}
	es[1] /= 2.0;
	es[4] /= 2.0;
	es[5] /= 2.0;
	es[0] = es[1] + es[2] + es[3] + es[4] + es[5];
	return es;
}

void Filter::printislandatom(std::vector<atominres> &ar, std::string lb, std::ostream &os) {
	for(auto &a:ar) os <<'\t' <<a.first[0] <<'\t' <<a.first[1] <<'\t' <<a.second;
	Ene_Atom & ea = enes_atom[ar[0].first[0]][ar[0].first[1]].at(ar[0].second);
	std::vector<std::pair<std::vector<atominres>,double>> &eterm=(
			lb=="bond"?ea.ebond:(lb=="ang"?ea.eang:(lb=="impdih"?ea.eimpdih:ea.esteric)));
	std::set<atominres> ar1;
	for(auto &a:ar) ar1.insert(a);
	for(int i=0;i<eterm.size();i++) {
		bool isall{true};
		for(int j=0;j<eterm[i].first.size();j++) {
			if(ar1.find(eterm[i].first[j])==ar1.end()) {
				isall=false;
				break;
			}
		}
		if(isall) {
			os <<'\t' <<eterm[i].second <<std::endl;
			break;
		}
	}
}

void Filter::printislandres(std::vector<int> &rs, std::string lb, std::ostream &os) {
	for(auto &r:rs) os <<'\t' <<r;
	Ene_Res & er=enes_res[rs[0]][rs[1]];
	if(rs.size()==2) {
		os <<'\t' <<(lb=="phipsi"?er.ephipsi:(lb=="local"?er.elocal:er.escconf)) <<std::endl;
	} else {
		std::map<std::pair<int,int>,double> &ep=(lb=="sitepair"?er.esitepair:(lb=="localbbhb"?er.elocalbbhb:er.escpacking));
		for(auto &p:ep) {
			if(p.first.first==rs[2] && p.first.second==rs[3]) {
				os <<'\t' <<p.second <<std::endl;
				break;
			}
		}
	}
}

void Filter::printisland(std::ostream &os) {
	os <<"bond:\t" <<highregion_.ibond.size() <<std::endl;
	for(int i=0;i<highregion_.ibond.size();i++) {
		os <<"island " <<i <<"   " <<highregion_.ibond[i].size() <<std::endl;
		for(int j=0;j<highregion_.ibond[i].size();j++) {
			int ix=highregion_.ibond[i][j];
			printislandatom(highregion_.ubond[ix],"bond",os);
		}
	}
	os <<std::endl;
	os <<"ang:\t" <<highregion_.iang.size() <<std::endl;
	for(int i=0;i<highregion_.iang.size();i++) {
		os <<"island " <<i <<"   " <<highregion_.iang[i].size() <<std::endl;
		for(int j=0;j<highregion_.iang[i].size();j++) {
			int ix=highregion_.iang[i][j];
			printislandatom(highregion_.uang[ix],"ang",os);
		}
	}
	os <<std::endl;
	os <<"impdih:\t" <<highregion_.iimpdih.size() <<std::endl;
	for(int i=0;i<highregion_.iimpdih.size();i++) {
		os <<"island " <<i <<"   " <<highregion_.iimpdih[i].size() <<std::endl;
		for(int j=0;j<highregion_.iimpdih[i].size();j++) {
			int ix=highregion_.iimpdih[i][j];
			printislandatom(highregion_.uimpdih[ix],"impdih",os);
		}
	}
	os <<std::endl;
	os <<"steric:\t" <<highregion_.isteric.size() <<std::endl;
	for(int i=0;i<highregion_.isteric.size();i++) {
		os <<"island " <<i <<"   " <<highregion_.isteric[i].size() <<std::endl;
		for(int j=0;j<highregion_.isteric[i].size();j++) {
			int ix=highregion_.isteric[i][j];
			printislandatom(highregion_.usteric[ix],"steric",os);
		}
	}
	os <<std::endl;
	os <<"phipsi:\t" <<highregion_.iphipsi.size() <<std::endl;
	for(int i=0;i<highregion_.iphipsi.size();i++) {
		os <<"island " <<i <<"   " <<highregion_.iphipsi[i].size() <<std::endl;
		for(int j=0;j<highregion_.iphipsi[i].size();j++) {
			int ix=highregion_.iphipsi[i][j];
			printislandres(highregion_.uphipsi[ix],"phipsi",os);
		}
	}
	os <<std::endl;
	os <<"local:\t" <<highregion_.ilocal.size() <<std::endl;
	for(int i=0;i<highregion_.ilocal.size();i++) {
		os <<"island " <<i <<"   " <<highregion_.ilocal[i].size() <<std::endl;
		for(int j=0;j<highregion_.ilocal[i].size();j++) {
			int ix=highregion_.ilocal[i][j];
			printislandres(highregion_.ulocal[ix],"local",os);
		}
	}
	os <<std::endl;
	os <<"localbbhb:\t" <<highregion_.ilocalbbhb.size() <<std::endl;
	for(int i=0;i<highregion_.ilocalbbhb.size();i++) {
		os <<"island " <<i <<"   " <<highregion_.ilocalbbhb[i].size() <<std::endl;
		for(int j=0;j<highregion_.ilocalbbhb[i].size();j++) {
			int ix=highregion_.ilocalbbhb[i][j];
			printislandres(highregion_.ulocalbbhb[ix],"localbbhb",os);
		}
	}
	os <<std::endl;
	os <<"sitepair:\t" <<highregion_.isitepair.size() <<std::endl;
	for(int i=0;i<highregion_.isitepair.size();i++) {
		os <<"island " <<i <<"   " <<highregion_.isitepair[i].size() <<std::endl;
		for(int j=0;j<highregion_.isitepair[i].size();j++) {
			int ix=highregion_.isitepair[i][j];
			printislandres(highregion_.usitepair[ix],"sitepair",os);
		}
	}
	os <<std::endl;
}
*/





