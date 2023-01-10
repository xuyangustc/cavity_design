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

std::vector<std::string> getenergynames() {
	return {"total","bond","ang","impdih","steric","scconf","phipsi","local","scpacking","localbbhb","sitepair"};
}

std::vector<std::string> Filter::eneterm=getenergynames();

void Filter::assign_ene(std::string label,
		const std::vector<std::pair<std::vector<atominres>,double>> &es) {
	int size=es[0].first.size();
	for(auto &e:es) {
		std::vector<int> vi;
		std::vector<int> vj;
		std::vector<std::string> vs;
		auto &v=e.first;
		for(int a=0;a<size;a++) {
			vi.push_back(v[a].first[0]);
			vj.push_back(v[a].first[1]);
			vs.push_back(v[a].second);
		}
		for(int a=0;a<size;a++) {
			if(enes_atom[vi[a]][vj[a]][vs[a]].atomname.empty()) {
				enes_atom[vi[a]][vj[a]][vs[a]].chainid=vi[a];
				enes_atom[vi[a]][vj[a]][vs[a]].resid=vj[a];
				enes_atom[vi[a]][vj[a]][vs[a]].atomname=vs[a];
			}
			std::vector<std::pair<std::vector<atominres>,double>> & esub=(
					label=="bond"?enes_atom[vi[a]][vj[a]][vs[a]].ebond:(
							label=="ang"?enes_atom[vi[a]][vj[a]][vs[a]].eang:(
									label=="impdih"?enes_atom[vi[a]][vj[a]][vs[a]].eimpdih:
											enes_atom[vi[a]][vj[a]][vs[a]].esteric)));
			esub.push_back(std::pair<std::vector<atominres>,double>());
			atominres air={{vi[a],vj[a]},vs[a]};
			for(int b=0;b<size;b++) {
				if(air==e.first[b]) continue;
				esub.back().first.push_back(e.first[b]);
			}
			esub.back().second=e.second;
		}
	}
}

void Filter::assign_ene(std::string label,
		const std::vector<std::pair<std::vector<std::vector<int>>,double>> &es) {
	//int size=es[0].first.size();
	int size=2;
	for(auto &e:es) {
		std::vector<int> vi;
		std::vector<int> vj;
		auto &v=e.first;
		for(int a=0;a<size;a++) {
			vi.push_back(v[a][0]);
			vj.push_back(v[a][1]);
		}
		for(int a=0;a<size;a++) {
			if(enes_res[vi[a]][vj[a]].chainid!=vi[a] || enes_res[vi[a]][vj[a]].resid!=vj[a]) {
				enes_res[vi[a]][vj[a]].chainid=vi[a];
				enes_res[vi[a]][vj[a]].resid=vj[a];
			}
			std::map<std::pair<int,int>,double> & esub=(
					label=="sitepair"?enes_res[vi[a]][vj[a]].esitepair:(
							label=="scpacking"?enes_res[vi[a]][vj[a]].escpacking:
									enes_res[vi[a]][vj[a]].elocalbbhb));
			for(int b=0;b<size;b++) {
				if(a==b) continue;
				esub.insert({{vi[b],vj[b]},e.second});
			}
		}
	}
}

void Filter::assign_ene(std::string label,
		const std::vector<std::pair<std::vector<int>,double>> &es) {
	for(auto &e:es) {
		int i0=e.first[0];
		int j0=e.first[1];
		(label=="phipsi"?enes_res[i0][j0].ephipsi:(
				label=="local"?enes_res[i0][j0].elocal:
						enes_res[i0][j0].escconf)) = e.second;
	}
}

void Filter::totalene_site() {
	for(int i=0;i<enes_res.size();i++) {
		for(int j=0;j<enes_res[i].size();j++) {
			if(fabs(enes_res[i][j].escconf-BigValue)>0.00001) enes_res[i][j].etot += enes_res[i][j].escconf;
			if(fabs(enes_res[i][j].ephipsi-BigValue)>0.00001) enes_res[i][j].etot += enes_res[i][j].ephipsi;
			if(fabs(enes_res[i][j].elocal-BigValue)>0.00001) enes_res[i][j].etot += enes_res[i][j].elocal;
			for(auto &p:enes_res[i][j].escpacking) {
				enes_res[i][j].etot += p.second / 2.0;
				enes_res[p.first.first][p.first.second].etot += p.second / 2.0;
			}
			for(auto &p:enes_res[i][j].elocalbbhb) {
				enes_res[i][j].etot += p.second / 2.0;
				enes_res[p.first.first][p.first.second].etot += p.second / 2.0;
			}
			for(auto &p:enes_res[i][j].esitepair) {
				enes_res[i][j].etot += p.second / 2.0;
				enes_res[p.first.first][p.first.second].etot += p.second / 2.0;
			}
		}
	}
	for(int i=0;i<enes_atom.size();i++) {
		for(int j=0;j<enes_atom[i].size();j++) {
			for(auto &p1:enes_atom[i][j]) {
				for(auto &p2:p1.second.ebond) {
					enes_res[i][j].etot += p2.second / 2.0;
					enes_res[p2.first[0].first[0]][p2.first[0].first[1]].etot += p2.second / 2.0;
				}
				for(auto &p2:p1.second.eang) {
					enes_res[i][j].etot += p2.second / 3.0;
					enes_res[p2.first[0].first[0]][p2.first[0].first[1]].etot += p2.second / 3.0;
					enes_res[p2.first[1].first[0]][p2.first[1].first[1]].etot += p2.second / 3.0;
				}
				for(auto &p2:p1.second.eimpdih) {
					enes_res[i][j].etot += p2.second / 4.0;
					enes_res[p2.first[0].first[0]][p2.first[0].first[1]].etot += p2.second / 4.0;
					enes_res[p2.first[1].first[0]][p2.first[1].first[1]].etot += p2.second / 4.0;
					enes_res[p2.first[2].first[0]][p2.first[2].first[1]].etot += p2.second / 4.0;
				}
				for(auto &p2:p1.second.esteric) {
					enes_res[i][j].etot += p2.second / 2.0;
					enes_res[p2.first[0].first[0]][p2.first[0].first[1]].etot += p2.second / 2.0;
				}
			}
		}
	}
}

void Filter::getenes() {
	std::string defaultffparfile=NSPdataio::datafilename("DefaultFF.par");
	NSPdataio::ControlFile cf;
	cf.readfile(defaultffparfile);
	std::vector<std::string> ffcontrolines=cf.getcontrolines("ForceField");
	std::string ffcontrolname{"control_ff"};
	NSPsd::defineforcefieldcontrol(ffcontrolname,ffcontrolines);
	std::vector<std::vector<NSPproteinrep::BackBoneSite>> bsss(chs_.size());
	for(int i=0;i<chs_.size();i++) {
		for(int j=0;j<chs_[i].size();j++) {
			bsss[i].push_back(chs_[i][j].getbackbonesite());
		}
	}
	NSPsd::ForceField ff=NSPsd::make_forcefield_allatom(bsss, ffcontrolname);
	std::string allfixed, mcfixed;
	NSPsd::ActiveSelections acts(&ff,allfixed,mcfixed);
	std::vector<std::vector<NSPproteinrep::FullSite>> fsss(chs_.size());
	for(int i=0;i<chs_.size();i++) {
		for(int j=0;j<chs_[i].size();j++) {
			fsss[i].push_back(chs_[i][j].getfullsite());
		}
	}
	std::vector<std::pair<std::vector<int>,std::string>> atomseq;
	std::vector<double> cs=NSPsd::GenChain::getcrd(fsss,atomseq);
	for(double & c:cs) c=c*A2NM;
	NSPsd::NeighborList nbl(cs,ff);
	std::vector<std::pair<std::vector<std::pair<std::vector<int>,std::string>>,double>> esubbond;
	std::vector<std::pair<std::vector<std::pair<std::vector<int>,std::string>>,double>> esubang;
	std::vector<std::pair<std::vector<std::pair<std::vector<int>,std::string>>,double>> esubimpdih;
	std::vector<std::pair<std::vector<std::pair<std::vector<int>,std::string>>,double>> esubsteric;
	std::vector<std::pair<std::vector<int>,double>> esubphipsi;
	std::vector<std::pair<std::vector<int>,double>> esublocal;
	std::vector<std::pair<std::vector<int>,double>> esubscconf;
	std::vector<std::pair<std::vector<std::vector<int>>,double>> esubsitepair;
	std::vector<std::pair<std::vector<std::vector<int>>,double>> esubscpacking;
	std::vector<std::pair<std::vector<std::vector<int>>,double>> esublocalbbhb;
	std::vector<double> pots;
	std::vector<int> sizes;
	for(int i=0;i<chs_.size();i++) sizes.push_back(chs_[i].size());
	ff.forces_subentry(cs, nbl, &pots, acts, atomseq, esubbond, esubang, esubimpdih, esubsteric,
			esubphipsi, esublocal, esubscconf, esubsitepair, esubscpacking, esublocalbbhb, sizes);
	enes_res.resize(chs_.size());
	enes_atom.resize(chs_.size());
	for(int i=0;i<chs_.size();i++) {
		enes_res[i].resize(chs_[i].size());
		enes_atom[i].resize(chs_[i].size());
	}
	assign_ene("bond",esubbond);
	assign_ene("ang",esubang);
	assign_ene("impdih",esubimpdih);
	assign_ene("steric",esubsteric);
	assign_ene("phipsi",esubphipsi);
	assign_ene("local",esublocal);
	assign_ene("scconf",esubscconf);
	assign_ene("sitepair",esubsitepair);
	assign_ene("scpacking",esubscpacking);
	assign_ene("localbbhb",esublocalbbhb);
	totalene_site();
	totalene_eterm();


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
	std::cout <<"SCConf:\t" <<pots[NSPsd::ForceField::ESCCONF] <<"   " <<ps[NSPsd::ForceField::ESCCONF] <<std::endl;
	std::cout <<"SCPacking:\t" <<pots[NSPsd::ForceField::ESCPACKING] <<"   " <<ps[NSPsd::ForceField::ESCPACKING] <<std::endl;
	std::cout <<"LocalBBHB:\t" <<pots[NSPsd::ForceField::ELOCALHB] <<"   " <<ps[NSPsd::ForceField::ELOCALHB] <<std::endl;
	std::cout <<std::endl;*/
}

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
}

void Filter::island_atom(std::string label) {
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
	island_space(ps,discut,ild);
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
}

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
}

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
	island_sequence(highene,discut,ild);
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
}

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
}

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
	island_sitepair(lb1.at(eneterm[ENE::SITEPAIR]).first, lb1.at(eneterm[ENE::SITEPAIR]).second);

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
}

void Filter::eneisland1() {
	std::map<std::string,std::pair<double,double>> lb1{{"bond",highregion_.pbond},{"ang",highregion_.pang},
			{"impdih",highregion_.pimpdih},{"steric",highregion_.psteric},{"scpacking",highregion_.pscpacking},
			{"sitepair",highregion_.psitepair}};
	std::map<std::string,std::pair<double,int>> lb2{{"scconf",highregion_.pscconf},{"phipsi",highregion_.pphipsi},
			{"local",highregion_.plocal},{"localbbhb",highregion_.plocalbbhb}};
	eneisland(lb1,lb2);
}

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
}
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








void Filter::gethbps() {
	std::set<std::pair<std::string,std::string>> drs, acs;
	std::set<std::pair<char,std::string>> drs1=AA_Atom::donors();
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
	//check hb in the same residue
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
					double d3=d1+d2;
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
								double d1=rads.at({hbps_[i0][j0].resname,h0.first})+rads.at({hbps_[i1][j1].resname,h1.first});
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

void Filter::mainchainHBfilter(std::vector<std::vector<int>>&helixs, std::vector<std::vector<int>> &hbh,
		std::vector<std::vector<int>>&strands, std::vector<std::vector<int>> &hbs) {
	for(auto &p:helixs) {
		hbh.push_back(hbinhelix(p[0],p[1],p[2]));
	}
	for(auto &p:strands) {
		hbs.push_back(hbinstrand(p[0],p[1],p[2],p[3],p[4],p[5]));
	}
}

std::vector<std::vector<double>> Filter::mainchainHBfilter(std::vector<std::vector<int>>&helixs,
		std::vector<std::vector<int>>&strands) {
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






