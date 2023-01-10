/*
 * polaratompacking.cpp
 *
 *  Created on: Apr 25, 2019
 *      Author: xuyang
 */
#include "allatom/polaratompacking.h"
using namespace NSPallatom;
using namespace NSPgeometry;

XYZ PolarAtomGroup::center() {
	XYZ cen;
	if(inMC) {
		cen = (crds.find("N")==crds.end())?crds.at("O"):crds.at("N");
	} else {
		if(resname=="SER") {
			cen = crds.at("OG");
		} else if(resname=="THR") {
			cen = crds.at("OG1");
		} else if(resname=="ASP") {
			cen = (crds.at("OD1")+crds.at("OD2"))/2.0;
		} else if(resname=="GLU") {
			cen = (crds.at("OE1")+crds.at("OE2"))/2.0;
		} else if(resname=="ASN") {
			cen = (crds.at("OD1")+crds.at("ND2"))/2.0;
		} else if(resname=="GLN") {
			cen = (crds.at("OE1")+crds.at("NE2"))/2.0;
		} else if(resname=="LYS") {
			cen = crds.at("NZ");
		} else if(resname=="ARG") {
			cen = (crds.at("NE")+crds.at("NH1")+crds.at("NH2"))/3.0;
		} else if(resname=="HIS") {
			cen = (crds.at("ND1")+crds.at("NE2"))/2.0;
		} else if(resname=="TYR") {
			cen = crds.at("OH");
		} else if(resname=="TRP") {
			cen = crds.at("NE1");
		} else {
			std::cout <<"Can Not Find PolarAtomGroup Center:   " <<resname <<std::endl;
			exit(1);
		}
	}
	return cen;
}

XYZ PolarAtomGroup::center_all() {
	XYZ cen(0.0,0.0,0.0);
	for(auto &c:crds) cen = cen + c.second;
	cen = cen / (double)(crds.size());
	return cen;
}

LocalFrame PolarAtomGroup::buildlocalframe() {
	XYZ origin=center(), px, pxy;
	if(inMC) {
		px = crds.at("C") - origin;
		pxy = crds.at("CA") - origin;
	} else {
		if(resname=="SER") {
			px = crds.at("CB") - origin;
			pxy = crds.at("CA") - origin;
		} else if(resname=="THR") {
			px = crds.at("CB") - origin;
			pxy = crds.at("CA") - origin;
		} else if(resname=="ASP") {
			px = crds.at("OD1") - origin;
			pxy = crds.at("CG") - origin;
		} else if(resname=="GLU") {
			px = crds.at("OE1") - origin;
			pxy = crds.at("CD") - origin;
		} else if(resname=="ASN") {
			px = crds.at("OD1") - origin;
			pxy = crds.at("CG") - origin;
		} else if(resname=="GLN") {
			px = crds.at("OE1") - origin;
			pxy = crds.at("CD") - origin;
		} else if(resname=="LYS") {
			px = crds.at("CE") - origin;
			pxy = crds.at("CD") - origin;
		} else if(resname=="ARG") {
			px = crds.at("NH1") - origin;
			pxy = crds.at("NE") - origin;
		} else if(resname=="HIS") {
			px = crds.at("ND1") - origin;
			pxy = crds.at("CG") - origin;
		} else if(resname=="TYR") {
			px = crds.at("CZ") - origin;
			pxy = crds.at("CE1") - origin;
		} else if(resname=="TRP") {
			px = crds.at("CE2") - origin;
			pxy = crds.at("CD1") - origin;
		} else {
			std::cout <<"Can Not Find PolarAtomGroup:   " <<resname <<std::endl;
			exit(1);
		}
	}
	return make_localframe(origin,px,pxy);
}

void PolarAtomGroup::global2local(const LocalFrame &lf) {
	for(auto &p:crds) {
		p.second = lf.global2localcrd(p.second);
	}
}

void PolarAtomGroup::findunits() {
	if(crds.find("N")!=crds.end()) {
		units.insert({"N","C","CA"});
		polaratoms.insert("N");
	} else if(crds.find("O")!=crds.end()) {
		units.insert({"O","C","CA"});
		polaratoms.insert("O");
	} else {
		std::map<std::string,std::map<std::string,std::set<std::vector<std::string>>>> scus=AA_Atom::sidechainHBunits();
		auto & sc = scus.at(resname);
		for(auto &p:sc) {
			polaratoms.insert(p.first);
			for(auto &v:p.second) {
				units.insert({p.first,v[0],v[1]});
			}
		}
	}
}

std::vector<HBUnit> PolarAtomGroup::extractHBunits() {
	if(units.empty()) findunits();
	std::vector<HBUnit> us;
	for(auto &u:units) {
		us.push_back(HBUnit(crds.at(u[0]),crds.at(u[1]),crds.at(u[2]),resname,u[0],u[1],u[2],0,0));
	}
	return us;
}

void PolarAtomGroup::translate(const XYZ &lth) {
	for(auto &p:crds) p.second = p.second + lth;
}

void PolarAtomGroup::rotate(const Rotation &rot) {
	for(auto &p:crds) rot.apply(&(p.second));
}

void PolarAtomGroupPair::extracthbunits(
		std::vector<std::shared_ptr<HBUnit>> &u1s, std::vector<std::shared_ptr<HBUnit>> &u2s) {
	;
}

void PolarAtomGroupPair::extracthbpairs(std::vector<HBUnitPair> &hbps, bool adjusthbp) {
	;
}

void PolarAtomPacking::extractpagsinSC(std::vector<std::shared_ptr<PolarAtomGroup>> &gs) {
	std::map<std::string,std::set<std::string>> resas;
	std::map<std::string,std::map<std::string,std::set<std::vector<std::string>>>> scus=AA_Atom::sidechainHBunits();
	for(auto &p:scus) {
		std::set<std::string> ss;
		for(auto &p1:p.second) {
			ss.insert(p1.first);
			for(auto &v:p1.second) {
				for(auto &s:v) {
					ss.insert(s);
				}
			}
		}
		resas.insert({p.first,ss});
	}

	for(int i=0;i<chains_.size();i++) {
		for(int j=0;j<chains_[i].size();j++) {
			std::string rn = chains_[i][j].resname();
			if(resas.find(rn)==resas.end()) continue;
			std::map<std::string,XYZ> cs;
			std::set<std::string> &as = resas.at(rn);
			for(auto &a:as) cs.insert({a,chains_[i][j].rds().at(a).crd});
			gs.push_back(std::shared_ptr<PolarAtomGroup>(new PolarAtomGroup(cs,rn,i,j)));
		}
	}
}

void PolarAtomPacking::extractpagsinMC(std::vector<std::shared_ptr<PolarAtomGroup>> &gs) {
	for(int i=0;i<chains_.size();i++) {
		for(int j=0;j<chains_[i].size();j++) {
			std::string rn = chains_[i][j].resname();
			auto & rs = chains_[i][j].rds();
			if(j!=0) { // N, C, CA
				std::map<std::string,XYZ> cs;
				cs.insert({"N",rs.at("N").crd});
				cs.insert({"C",chains_[i][j-1].rds().at("C").crd});
				cs.insert({"CA",rs.at("CA").crd});
				gs.push_back(std::shared_ptr<PolarAtomGroup>(new PolarAtomGroup(cs,rn,i,j)));
				gs.back()->inMC = true;
			}
			if(j!=chains_[i].size()-1) { // O, C, CA
				std::map<std::string,XYZ> cs;
				cs.insert({"O",rs.at("O").crd});
				cs.insert({"C",rs.at("C").crd});
				cs.insert({"CA",rs.at("CA").crd});
				gs.push_back(std::shared_ptr<PolarAtomGroup>(new PolarAtomGroup(cs,rn,i,j)));
				gs.back()->inMC = true;
			}
		}
	}
}

void PolarAtomPacking::findgpps() {
	std::vector<std::shared_ptr<PolarAtomGroup>> gs;
	extractpagsinMC(gs);
	//std::cout <<gs.size() <<std::endl;
	extractpagsinSC(gs);
	//std::cout <<gs.size() <<std::endl;

	for(auto &g:gs) g->findunits();

	//std::map<std::string,std::set<std::string>> resas;
	std::map<std::string,std::map<std::string,std::set<std::vector<std::string>>>> scus=AA_Atom::sidechainHBunits();

	double dis2 = dcutoff_ * dcutoff_;

	for(int i=0;i<gs.size();i++) {
		for(int j=i+1;j<gs.size();j++) {
			if(gs[i]->crds.find("N")!=gs[i]->crds.end() && gs[j]->crds.find("N")!=gs[j]->crds.end()) continue;
			if(gs[i]->crds.find("O")!=gs[i]->crds.end() && gs[j]->crds.find("O")!=gs[j]->crds.end()) continue;
			if(gs[i]->chainseq==gs[j]->chainseq) {
				if((gs[i]->crds.find("N")!=gs[i]->crds.end() &&
						gs[j]->crds.find("O")!=gs[j]->crds.end()) || (
								gs[i]->crds.find("O")!=gs[i]->crds.end() &&
								gs[j]->crds.find("N")!=gs[j]->crds.end())) {
					int delta = fabs(gs[i]->resseq-gs[j]->resseq);
					if(delta<2) continue;
				}
			}
			bool contanct{false};
			for(auto &a1:gs[i]->polaratoms) {
				//std::cout <<a1 <<'\t' <<gs[i]->resname <<'\t' <<gs[i]->chainseq <<'\t' <<gs[i]->resseq <<'\t' <<gs[i]->crds.size() <<std::endl;
				XYZ c1 = gs[i]->crds.at(a1);
				for(auto &a2:gs[j]->polaratoms) {
					//std::cout <<a2 <<'\t' <<gs[j]->resname <<'\t' <<gs[j]->chainseq <<'\t' <<gs[j]->resseq <<'\t' <<gs[j]->crds.size() <<std::endl;
					XYZ c2 = gs[j]->crds.at(a2);
					double d2 = (c1-c2).squarednorm();
					if(d2>dis2) continue;
					contanct = true;
					break;
				}
			}
			if(!contanct) continue;
			gpairs_.push_back(PolarAtomGroupPair(gs[i],gs[j]));
		}
	}
}

void PolarAtomPacking::extracthbunitsinSC(std::vector<std::vector<std::shared_ptr<HBUnit>>> &us) {
	//int nt=0;
	std::map<std::string,std::map<std::string,std::set<std::vector<std::string>>>> scus=AA_Atom::sidechainHBunits();
	/*int nt=0;
	for(auto &p1:scus) {
		for(auto &p2:p1.second) {
			nt+=p2.second.size();
		}
	}
	std::cout <<"Number of HBUnit: " <<nt <<std::endl;*/

	//std::vector<std::vector<HBUnit>> us;
	for(int i=0;i<chains_.size();i++) {
		for(int j=0;j<chains_[i].size();j++) {
			std::vector<std::shared_ptr<HBUnit>> us1;
			std::string rn = chains_[i][j].resname();
			if(scus.find(rn)==scus.end()) continue;
			for(auto &p1:scus.at(rn)) {
				std::string p0 = p1.first;
				if(chains_[i][j].rds().find(p0)==chains_[i][j].rds().end()) continue;
				XYZ cp0 = chains_[i][j].rds().at(p0).crd;
				for(auto &p2:p1.second) {
					std::string c1 = p2[0];
					std::string c2 = p2[1];
					if(chains_[i][j].rds().find(c1)==chains_[i][j].rds().end()) continue;
					if(chains_[i][j].rds().find(c2)==chains_[i][j].rds().end()) continue;
					XYZ cc1 = chains_[i][j].rds().at(c1).crd;
					XYZ cc2 = chains_[i][j].rds().at(c2).crd;
					us1.push_back(std::shared_ptr < HBUnit > (new HBUnit(cp0,cc1,cc2,rn,p0,c1,c2,i,j)));
					//nt++;
				}
			}
			us.push_back(us1);
		}
	}
	//std::cout <<"SC Unit: " <<nt <<std::endl;

	//return us;
}

void PolarAtomPacking::extracthbunitsinMC(std::vector<std::vector<std::shared_ptr<HBUnit>>> &us) {
	//std::vector<std::vector<HBUnit>> us;

	for(int i=0;i<chains_.size();i++) {
		us.push_back(std::vector<std::shared_ptr<HBUnit>>());
		us.back().push_back(std::shared_ptr < HBUnit > (new HBUnit(chains_[i][0].rds().at("O").crd,
				chains_[i][0].rds().at("C").crd,
				chains_[i][0].rds().at("CA").crd,
				chains_[i][0].resname(),"O","C","CA",i,0)));
		int j=1;
		for(;j<chains_[i].size()-1;j++) {
			us.push_back(std::vector<std::shared_ptr<HBUnit>>());
			us.back().push_back(std::shared_ptr < HBUnit > (new HBUnit(chains_[i][j].rds().at("N").crd,
					chains_[i][j-1].rds().at("C").crd,
					chains_[i][j].rds().at("CA").crd,
					chains_[i][j].resname(),"N","C","CA",i,j)));
			us.back().push_back(std::shared_ptr < HBUnit > (new HBUnit(chains_[i][j].rds().at("O").crd,
					chains_[i][j].rds().at("C").crd,
					chains_[i][j].rds().at("CA").crd,
					chains_[i][j].resname(),"O","C","CA",i,j)));
		}
		us.push_back(std::vector<std::shared_ptr<HBUnit>>());
		us.back().push_back(std::shared_ptr < HBUnit > (new HBUnit(chains_[i][j].rds().at("N").crd,
				chains_[i][j-1].rds().at("C").crd,
				chains_[i][j].rds().at("CA").crd,
				chains_[i][j].resname(),"N","C","CA",i,j)));
	}

	//return us;
}

std::vector<std::vector<std::shared_ptr<HBUnit>>> PolarAtomPacking::findhbps(bool adjusthbp) {
	double dis2 = dcutoff_ * dcutoff_;

	std::vector<std::vector<std::shared_ptr<HBUnit>>> us;
	extracthbunitsinMC(us);
	//std::cout <<us.size() <<std::endl;
	extracthbunitsinSC(us);
	//std::cout <<us.size() <<std::endl;

	if(adjusthbp) {
		for(auto &vu:us) {
			for(auto &u:vu) {
				if(u->np0=="N") u->resname = "MC";
				if(u->np0=="O") u->resname = "MC";
			}
		}
		for(auto &vu:us) {
			for(auto &u:vu) {
				if(u->resname=="ASP" && u->np0=="OD1") u->np0="OD";
				else if(u->resname=="ASP" && u->np0=="OD2") u->np0="OD";
				else if(u->resname=="GLU" && u->np0=="OE1") u->np0="OE";
				else if(u->resname=="GLU" && u->np0=="OE2") u->np0="OE";
				else if(u->resname=="ARG" && u->np0=="NH1") u->np0="NH";
				else if(u->resname=="ARG" && u->np0=="NH2") u->np0="NH";
				else if(u->resname=="TYR" && u->nc2=="CE1") u->nc2="CE";
				else if(u->resname=="TYR" && u->nc2=="CE2") u->nc2="CE";
			}
		}
	}

	for(int i=0;i<us.size();i++) {
		for(int j=i+1;j<us.size();j++) {
			for(int i1=0;i1<us[i].size();i1++) {
				std::shared_ptr<HBUnit> & u1 = us[i][i1];
				for(int j1=0;j1<us[j].size();j1++) {
					std::shared_ptr<HBUnit> & u2 = us[j][j1];
					if(u1->np0=="N" && u2->np0=="N") continue;
					if(u1->np0=="O" && u2->np0=="O") continue;
					if((u1->resseq+1==u2->resseq || u1->resseq==u2->resseq+1) && u1->chainseq==u2->chainseq) {
						if(u1->np0=="N" && u2->np0=="O") continue;
						if(u1->np0=="O" && u2->np0=="N") continue;
					}
					XYZ pi = us[i][i1]->p0;
					XYZ pj = us[j][j1]->p0;
					double d2 = (pi-pj).squarednorm();
					if(d2>dis2) continue;
					hbpairs_.push_back(HBUnitPair(us[i][i1],us[j][j1]));
				}
			}
		}
	}

	return us;


	//int nmchp=0;
/*	std::vector<std::vector<HBUnit>> mcus = gethbunitsinMC();
	//for(auto &vi:mcus) nmchp+=vi.size();
	//std::cout <<"unit mc : " <<nmchp <<std::endl;
	//nmchp=0;
	for(int i=0;i<mcus.size();i++) {
		for(int j=i+1;j<mcus.size();j++) {
			if(j==i+1 && (mcus[i].size()!=1 || mcus[j].size()!=1)) continue;
			for(int i1=0;i1<mcus[i].size();i1++) {
				for(int j1=0;j1<mcus[j].size();j1++) {
					if(mcus[i][i1].np0==mcus[j][j1].np0) continue;
					XYZ pi = mcus[i][i1].p0;
					XYZ pj = mcus[j][j1].p0;
					double d2 = (pi-pj).squarednorm();
					//nmchp++;
					//std::cout <<i <<'\t' <<mcus[i][i1].np0 <<'\t' <<j <<'\t' <<mcus[j][j1].np0 <<std::endl;
					if(d2>dis2) continue;
					hbps_.push_back(HBUnitPair(mcus[i][i1],mcus[j][j1]));
				}
			}
		}
	}*/
	//std::cout <<"Unit Pair in MC: " <<nmchp <<"   " <<hbps_.size() <<std::endl;


/*	std::vector<std::vector<HBUnit>> scus = gethbunitsinSC();
	for(int i=0;i<scus.size();i++) {

		std::vector<XYZ> cs1;
		for(HBUnit &hu:scus[i]) cs1.push_back(hu.p0);

		for(int i1=0;i1<mcus.size();i1++) {
			for(int j1=0;j1<mcus[i1].size();j1++) {
				XYZ o1 = mcus[i1][j1][0].p0;
				XYZ n1 = mcus[i1][j1][1].p0;
				for(int j=0;j<cs1.size();j++) {
					double d2 = (n1-cs1[j]).squarednorm();
					if(d2<dis2) {
						hbps_.push_back(HBUnitPair(mcus[i1][j1][1],scus[i][j]));
					}
					d2 = (o1-cs1[j]).squarednorm();
					if(d2<dis2) {
						hbps_.push_back(HBUnitPair(mcus[i1][j1][0],scus[i][j]));
					}
				}
			}
		}

		for(int j=i+1;j<scus.size();j++) {
			std::vector<XYZ> cs2;
			for(HBUnit &hu:scus[j]) cs2.push_back(hu.p0);

			for(int i1=0;i1<cs1.size();i1++) {
				for(int j1=0;j1<cs2.size();j1++) {
					double d2 = (cs1[i1]-cs2[j1]).squarednorm();
					if(d2<dis2) {*/
						/*HBUnitPair hp(scus[i][i1],scus[j][j1]);
						if(scus[i][i1].resname > scus[j][j1].resname) {
							hp.u1 = scus[j][j1];
							hp.u2 = scus[i][i1];
						} else if(scus[i][i1].resname == scus[j][j1].resname) {
							if(scus[i][i1].np0 > scus[j][j1].np0) {
								hp.u1 = scus[j][j1];
								hp.u2 = scus[i][i1];
							} else if(scus[i][i1].np0 == scus[j][j1].np0) {
								if(scus[i][i1].nc1 > scus[j][j1].nc1) {
									hp.u1 = scus[j][j1];
									hp.u2 = scus[i][i1];
								} else if(scus[i][i1].nc1 == scus[j][j1].nc1) {
									if(scus[i][i1].nc2 > scus[j][j1].nc2) {
										hp.u1 = scus[j][j1];
										hp.u2 = scus[i][i1];
									}
								}
							}
						}
						hbps_.push_back(hp);*/
/*						hbps_.push_back(HBUnitPair(scus[i][i1],scus[j][j1]));
					}
				}
			}
		}

	}*/
}

