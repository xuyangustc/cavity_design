/*
 * grouppacking.cpp
 *
 *  Created on: Jun 17, 2019
 *      Author: xuyang
 */
#include "dstl/randomengine.h"
#include "allatom/grouppacking.h"

using namespace NSPgeometry;
using namespace NSPallatom;

bool GroupPair::clash(const std::vector<std::string>&n1, const std::vector<XYZ>&cs1,
		const std::vector<std::string>&n2, const std::vector<XYZ>&cs2) {
	std::map<std::vector<std::vector<std::string>>,double> apc=AA_Atom::atompairclashdis();
	for(int i=0;i<cs1.size();i++) {
		std::vector<std::string> vs1{n1[0],n1[i+1]};
		if(vs1[0].substr(0,2)=="MC") vs1[0]="MC";
		for(int j=0;j<cs2.size();j++) {
			std::vector<std::string> vs2{n2[0],n2[j+1]};
			if(vs2[0].substr(0,2)=="MC") vs2[0]="MC";
			double d1 = apc.at({vs1,vs2});
			double d2 = d1 * d1;
			double d22 = (cs1[i]-cs2[j]).squarednorm();
			if(d2>d22) return true;
		}
	}
	return false;
}

void GroupPair::rerank() {
	std::vector<std::string> v1=g1->makelabel();
	std::vector<std::string> v2=g2->makelabel();
	if(v2<v1) changeorder();
}

void GroupPair::rerank_samegroup(int i, int j) {
	std::vector<std::string> v1=g1->makelabel();
	std::vector<std::string> v2=g2->makelabel();
	if(v1!=v2) return;
	XYZ c1=g1->crds_v[i].second-g2->crds_v[j].second;
	XYZ c2=g1->crds_v[j].second-g2->crds_v[i].second;
	if(c2.squarednorm()<c1.squarednorm()) changeorder();
}


void ChemGroup::rotate(const Rotation &rt) {
	for(auto &p:crds) {
		rt.apply(&(p.second));
	}
}
void ChemGroup::translate(const XYZ &md) {
	for(auto &p:crds) {
		p.second = p.second + md;
	}
}

bool GroupPair::clash() {
	std::map<std::vector<std::vector<std::string>>,double> apc=AA_Atom::atompairclashdis();
	for(auto &p1:g1->crds) {
		for(auto &p2:g2->crds) {
			//std::cout <<g1->resname <<p1.first <<g2->resname <<p2.first <<std::endl;
			std::vector<std::string> vs1{g1->resname,p1.first};
			if(vs1[0].substr(0,2)=="MC") vs1[0]="MC";
			std::vector<std::string> vs2{g2->resname,p2.first};
			if(vs2[0].substr(0,2)=="MC") vs2[0]="MC";
			double d1 = apc.at({vs1,vs2});
			double d2 = d1 * d1;
			double d22 = (p1.second-p2.second).squarednorm();
			if(d2>d22) return true;
		}
	}
	return false;
}

std::vector<double> GroupPair::atom3par5() {
	double dis=sqrt((g1->crds_v[1].second-g2->crds_v[1].second).squarednorm());
	double ang1=angle(g1->crds_v[0].second,g1->crds_v[1].second,g2->crds_v[1].second);
	double dih1=torsion(g1->crds_v[0].second,g1->crds_v[2].second,g1->crds_v[1].second,g2->crds_v[1].second);
	double ang2=angle(g2->crds_v[0].second,g2->crds_v[1].second,g1->crds_v[1].second);
	double dih2=torsion(g2->crds_v[0].second,g2->crds_v[2].second,g2->crds_v[1].second,g1->crds_v[1].second);
	return {dis,ang1,dih1,ang2,dih2};
}

std::vector<double> GroupPair::atom3par6() {
	LocalFrame lf1=make_localframe(g1->crds_v[1].second,g1->crds_v[0].second,g1->crds_v[2].second);
	LocalFrame lf2=make_localframe(g2->crds_v[1].second,g2->crds_v[0].second,g2->crds_v[2].second);
	XYZ c1=lf2.global2localcrd(g1->crds_v[1].second);
	XYZ c2=lf1.global2localcrd(g2->crds_v[1].second);
	return {c1.x_,c1.y_,c1.z_,c2.x_,c2.y_,c2.z_};
}

double GroupPair::minatomdis() {
	double d2=1000000.0;
	for(auto &p1:g1->crds_v) {
		for(auto &p2:g2->crds_v) {
			double d=(p1.second-p2.second).squarednorm();
			if(d<d2) d2=d;
		}
	}
	return sqrt(d2);
}









XYZ ChemGroup::RandomRotateWithFixedCenter(std::vector<XYZ>&cs) {
	XYZ cen1(0.0,0.0,0.0);
	for(XYZ &c:cs) {
		cen1 = cen1 + c;
	}
	cen1 = cen1/(double)(cs.size());
	auto &rng01=NSPdstl::RandomEngine<>::getinstance().realrng(0.0,1.0);
	QuaternionCrd qc(rng01,0);
	Rotation rt(qc,cen1);
	for(XYZ &c:cs) {
		rt.apply(&c);
	}
	return cen1;
}

void ChemGroup::findcenter() {
	std::set<std::string> lbs;
	for(auto &p:crds) lbs.insert(p.first);
	std::map<std::string,std::vector<std::set<std::string>>> gpns = AA_Atom::groups();
	std::map<std::string,std::vector<std::set<std::string>>> gpns_cen = AA_Atom::groups_center();
	std::set<std::string> cens;
	for(int i=0;i<gpns.at(resname).size();i++) {
		if(gpns.at(resname)[i]==lbs) {
			cens = gpns_cen.at(resname)[i];
			break;
		}
	}

	center = XYZ(0.0,0.0,0.0);
	for(auto &s:cens) {
		center = center + crds.at(s);
	}
	center = center / (double)(cens.size());
}

XYZ ChemGroup::findcenter(int n) {
	std::set<std::string> lbs;
	for(auto &p:decoys[n]) lbs.insert(p.first);
	std::map<std::string,std::vector<std::set<std::string>>> gpns = AA_Atom::groups();
	std::map<std::string,std::vector<std::set<std::string>>> gpns_cen = AA_Atom::groups_center();
	std::set<std::string> cens;
	for(int i=0;i<gpns.at(resname).size();i++) {
		if(gpns.at(resname)[i]==lbs) {
			cens = gpns_cen.at(resname)[i];
			break;
		}
	}

	XYZ cn(0.0,0.0,0.0);
	for(auto &s:cens) {
		cn = cn + decoys[n].at(s);
	}
	cn = cn / (double)(cens.size());

	return cn;
}

void ChemGroup::finddisatom() {
	std::map<std::string,std::vector<std::vector<std::string>>> aags = AA_Atom::groups(2);
	for(auto &v:aags.at(resname)) {
		bool isthis{true};
		for(auto &s:v) {
			if(crds.find(s)==crds.end()) {
				isthis=false;
				break;
			}
		}
		if(!isthis) continue;
		std::set<std::string> ss;
		for(auto &s:v) ss.insert(s);
		for(int i=0;i<crds_v.size();i++) {
			if(ss.find(crds_v[i].first)==ss.end()) continue;
			disatom.insert(i);
		}
		break;
	}
}

void GroupPacking::addgroup(Residue &res, std::vector<std::string> &ans, std::string resname, int cid, int rid) {
	std::shared_ptr<ChemGroup> cg = std::shared_ptr < ChemGroup > (new ChemGroup(cid,rid));
	cg->resname = resname;
	bool isall{true};
	for(auto &an:ans) {
		if(res.rds().find(an)==res.rds().end()) {
			isall = false;
			break;
		}
		cg->crds.insert({an,res.rds().at(an).crd});
	}
	//cg->findcenter();


	if(isall) {
		for(auto &an:ans) {
			cg->crds_v.push_back({an,res.rds().at(an).crd});
		}
		//cg->makecenter();
		//cg->findcrdatom();
		cg->finddisatom();
		cgs_.push_back(cg);
	}
}

void GroupPacking::findgroup() {
	//std::cout <<1 <<std::endl;
	std::map<std::string,std::vector<std::set<std::string>>> gpns = AA_Atom::groups();
	//std::cout <<2 <<std::endl;
	std::map<std::string,std::vector<std::vector<std::string>>> gpns_v = AA_Atom::groups_v();
	std::map<std::string,char> a31=AA_Atom::aa31();
	//std::cout <<3 <<std::endl;
	//std::map<std::string,std::vector<std::set<std::string>>> gpns_cen = AA_Atom::groups_center();
	std::set<std::string> mcas{"CA","C","N","O"};

	/*for(auto &p:gpns) {
		std::cout <<p.first <<'\t' <<p.second.size() <<std::endl;
		for(auto &v:p.second) {
			std::cout <<v.size() <<'\t';
		}
		std::cout <<std::endl;
	}*/

	for(int i=0;i<chains_.size();i++) {
		for(int j=0;j<chains_[i].size();j++) {
			//std::cout <<i <<' ' <<j <<' ' <<cgs_.size() <<std::endl;
			std::string rn = chains_[i][j].resname();
			if(gpns_v.find(rn)==gpns_v.end()) {
				if(a31.find(rn)==a31.end()) std::cout <<"Can Not Recognize Residue: " <<rn <<std::endl;
			} else {
				std::vector<std::vector<std::string>> &gs = gpns_v.at(rn);
				//std::vector<std::set<std::string>> &gs_cen = gpns_cen.at(rn);
				for(int k=0;k<gs.size();k++) {
					addgroup(chains_[i][j],gs[k],rn,i,j);
				}
			}
			/*if(gpns.find("MC")!=gpns.end()) {
				std::vector<std::set<std::string>> &gsm = gpns.at("MC");
				//std::vector<std::set<std::string>> &gsm_cen = gpns_cen.at("MC");
				for(int k=0;k<gsm.size();k++) {
					if(gsm[k].find("N")!=gsm[k].end() && gsm[k].find("O")!=gsm[k].end()) {
						if(j!=0) {
							std::shared_ptr<ChemGroup> cg = std::shared_ptr < ChemGroup > (new ChemGroup(i,j));
							cg->resname = "MC";
							for(auto &an:gsm[k]) {
								if(an=="C" || an=="O") {
									cg->crds.insert({an,chains_[i][j-1].rds().at(an).crd});
								} else cg->crds.insert({an,chains_[i][j].rds().at(an).crd});
							}
							cg->findcenter();
							cgs_.push_back(cg);
						}
					} else {
						addgroup(chains_[i][j],gsm[k],"MC",i,j);
					}
				}
			}
			if(gpns.find("MCE")!=gpns.end()) {
				std::vector<std::set<std::string>> &gse = gpns.at("MC");
				//std::vector<std::set<std::string>> &gse_cen = gpns_cen.at("MC");
				for(int k=0;k<gse.size();k++) {
					if(j==0 && gse[k].find("N")!=gse[k].end()) addgroup(chains_[i][j],gse[k],"MCE",i,j);
					if(j==chains_[i].size()-1 && gse[k].find("O")!=gse[k].end()) addgroup(chains_[i][j],gse[k],"MCE",i,j);
				}
			}*/
			if(j!=0) {
				if(gpns_v.find("MCCIS")==gpns_v.end() && gpns_v.find("MCTRANS")==gpns_v.end()) continue;
				double t=chains_[i][j].omega(chains_[i][j-1]);
				std::shared_ptr<ChemGroup> cg = std::shared_ptr < ChemGroup > (new ChemGroup(i,j));
				if(t>-90.0 && t<90.0) {
					if(gpns_v.find("MCCIS")!=gpns_v.end()) {
						cg->resname = "MCCIS";
					} else continue;
				} else {
					if(gpns_v.find("MCTRANS")!=gpns_v.end()) {
						cg->resname = "MCTRANS";
					} else continue;
				}
				//std::set<std::string> ans1{"CA","C","O"};
				//std::set<std::string> ans2{"N","CA"};
				//cg->crds.insert({"CA",chains_[i][j-1].rds().at("CA").crd});
				cg->crds.insert({"C",chains_[i][j-1].rds().at("C").crd});
				cg->crds.insert({"O",chains_[i][j-1].rds().at("O").crd});
				cg->crds.insert({"N",chains_[i][j].rds().at("N").crd});
				//cg->cacrd=chains_[i][j].rds().at("CA").crd;

				//cg->crds_v.push_back({"CA",chains_[i][j-1].rds().at("CA").crd});
				cg->crds_v.push_back({"N",chains_[i][j].rds().at("N").crd});
				cg->crds_v.push_back({"C",chains_[i][j-1].rds().at("C").crd});
				cg->crds_v.push_back({"O",chains_[i][j-1].rds().at("O").crd});
				//cg->crds_v.push_back({"CA",chains_[i][j].rds().at("CA").crd});

				//cg->findcrdatom();
				cg->finddisatom();
				cg->chainseq = i;
				cg->resseq = j;
				cgs_.push_back(cg);
			}

			/*for(auto &inf:chains_[i][j].rds()) {
				if(mcas.find(inf.first)!=mcas.end()) continue;
				for(auto &g:gs) {
					if(g[0].find(inf.first)==g[0].end()) continue;
					std::shared_ptr<ChemGroup> cg = std::shared_ptr < ChemGroup > (new ChemGroup(i,j));
					cg->resname = rn;
					for(auto &an:g[0]) {
						cg->crds.insert({an,chains_[i][j].rds().at(an).crd});
					}
					XYZ cn(0.0,0.0,0.0);
					double na = (double)(g[1].size());
					for(auto &an:g[1]) {
						cn = cn + chains_[i][j].rds().at(an).crd;
					}
					cg->center = cn / na;
					cgs_.push_back(cg);
				}
			}*/
			/*if(j!=0) {
				std::shared_ptr<ChemGroup> cg = std::shared_ptr < ChemGroup > (new ChemGroup(i,j));
				cg->resname = "MC";
				cg->crds.insert({"C",chains_[i][j-1].rds().at("C").crd});
				cg->crds.insert({"N",chains_[i][j].rds().at("N").crd});
				cg->crds.insert({"CA",chains_[i][j].rds().at("CA").crd});
				//cg->center = (cg->crds.at("O")+cg->crds.at("N")) / 2.0;
				cg->center = cg->crds.at("N");
				cgs_.push_back(cg);
			}
			std::shared_ptr<ChemGroup> cg = std::shared_ptr < ChemGroup > (new ChemGroup(i,j));
			cg->resname = "MC";
			cg->crds.insert({"CA",chains_[i][j].rds().at("CA").crd});
			cg->crds.insert({"C",chains_[i][j].rds().at("C").crd});
			cg->crds.insert({"O",chains_[i][j].rds().at("O").crd});
			cg->center = cg->crds.at("O");
			cgs_.push_back(cg);*/
		}
	}
}
/*
void GroupPacking::findgrouppair() {
	if(cgs_.empty()) findgroup();
	//double dc2 = dcutoff_*dcutoff_;
	std::set<std::string> mcas_noca{"C","N","O"};
	for(int i=0;i<cgs_.size();i++) {
		for(int j=i+1;j<cgs_.size();j++) {
			if(cgs_[i]->crds.find("C")!=cgs_[i]->crds.end() && cgs_[j]->crds.find("C")!=cgs_[j]->crds.end()) {
				continue;
			} else if(cgs_[i]->crds.find("C")==cgs_[i]->crds.end() && cgs_[j]->crds.find("C")==cgs_[j]->crds.end()) {
				if(cgs_[i]->chainseq==cgs_[j]->chainseq && cgs_[i]->resseq==cgs_[j]->resseq) {
					continue;
				} else if(cgs_[i]->chainseq==cgs_[j]->chainseq && cgs_[i]->resseq+1==cgs_[j]->resseq) {
					//GroupPair gp(cgs_[i],cgs_[j]);
					//if(gp.clash()) continue;
					//gps_scsc_1.push_back(gp);
					continue;
				} else if(cgs_[i]->chainseq==cgs_[j]->chainseq && cgs_[i]->resseq==cgs_[j]->resseq+1) {
					//GroupPair gp(cgs_[j],cgs_[i]);
					//if(gp.clash()) continue;
					//gps_scsc_1.push_back(gp);
					continue;
				} else {
					GroupPair gp(cgs_[i],cgs_[j]);
					gp.findmindis();
					//if(gp.dis>dcutoff_) continue;
					//if(gp.mindis<5.0 && gp.clash()) continue;
					if(gp.mindis>dcutoff_) continue;
					if(gp.mindis<4.0 && gp.clash()) continue;
					gp.rerank();
					gps_scsc_no1.push_back(gp);
				}
			} else {
				bool ms1,sm1,ms;
				ms1=sm1=ms=false;
				int a=i, b=j;
				if(cgs_[i]->crds.find("O")==cgs_[i]->crds.end() && cgs_[j]->crds.find("O")==cgs_[j]->crds.end()) {
					if(cgs_[j]->crds.find("C")!=cgs_[j]->crds.end()) {
						a=j;
						b=i;
					}
					if(cgs_[a]->chainseq==cgs_[b]->chainseq && cgs_[a]->resseq==cgs_[b]->resseq) continue;
					else if(cgs_[a]->chainseq==cgs_[b]->chainseq && cgs_[a]->resseq+1==cgs_[b]->resseq) ms1=true;
					else if(cgs_[a]->chainseq==cgs_[b]->chainseq && cgs_[a]->resseq==cgs_[b]->resseq+1) sm1=true;
					else ms=true;
				} else {
					if(cgs_[j]->crds.find("O")!=cgs_[j]->crds.end()) {
						a=j;
						b=i;
					}
					if(cgs_[a]->chainseq==cgs_[b]->chainseq && cgs_[a]->resseq==cgs_[b]->resseq) ms1=true;
					else if(cgs_[a]->chainseq==cgs_[b]->chainseq && cgs_[a]->resseq==cgs_[b]->resseq+1) sm1=true;
					else ms=true;
				}
				//if(ms1) gps_mcsc_1.push_back(GroupPair(cgs_[a],cgs_[b]));
				//if(sm1) gps_scmc_1.push_back(GroupPair(cgs_[b],cgs_[a]));
				if(ms) {
					GroupPair gp(cgs_[a],cgs_[b]);
					if(gp.dis>dcutoff_) continue;
					if(gp.dis<5.0 && gp.clash()) continue;
					gps_mcsc.push_back(gp);
				}

			}
			//if(cgs_[i]->resname=="MC" && cgs_[j]->resname=="MC") continue;
			//std::cout <<cgs_[i]->resseq <<' ' <<cgs_[j]->resseq <<"    ";
			//if(cgs_[i]->chainseq==cgs_[j]->chainseq) {
			//	int itv = fabs(cgs_[i]->resseq-cgs_[j]->resseq);
				//int itv = cgs_[i]->resseq==cgs_[j]->resseq;
				//std::cout <<itv <<std::endl;
			//	if(itv<=seqitv_) continue;
			//}
			//GroupPair gp(cgs_[i],cgs_[j]);
			//if(gp.dis<dcutoff_) gps_.push_back(GroupPair(cgs_[i],cgs_[j]));
		}
	}
}
*/


bool GroupPacking::gpfilter(GroupPair &gp, bool removeclash) {
	gp.findmindis();
	if(gp.mindis>dcutoff_) return false;
	if(removeclash) if(gp.mindis<4.0 && gp.clash()) return false;
	gp.rerank();
	return true;
}
void GroupPacking::findgrouppair(bool removemcmc, bool removeclash, bool samedouble) {
	if(cgs_.empty()) findgroup();
	//std::cout <<cgs_.size() <<std::endl;
	//double dc2 = dcutoff_*dcutoff_;
	for(int i=0;i<cgs_.size();i++) {
		for(int j=i+1;j<cgs_.size();j++) {
			if(removemcmc) {
				if(cgs_[i]->crds.find("C")!=cgs_[i]->crds.end()
						&& cgs_[j]->crds.find("C")!=cgs_[j]->crds.end()) continue;
			}
			if(!GroupPair::farinteraction(*(cgs_[i]), *(cgs_[j]), nbondmin_)) continue;
			GroupPair gp(cgs_[i],cgs_[j]);
			if(gpfilter(gp,removeclash)) {
				gps_all.push_back(gp);
				if(samedouble) {
					if(gp.samegroup()) {
						GroupPair gp1(cgs_[j],cgs_[i]);
						gp1.findmindis();
						gps_all.push_back(gp1);
					}
				}
			}
		}
	}
}








void ChemGroup::MakeDecoysWithFixedCenter(int num) {
	auto &rng01=NSPdstl::RandomEngine<>::getinstance().realrng(0.0,1.0);
	for(int i=0;i<num;i++) {
		QuaternionCrd qc(rng01,0);
		Rotation rt(qc,center);
		std::map<std::string,XYZ> cs = crds;
		for(auto &c:cs) {
			rt.apply(&(c.second));
		}
		decoys.push_back(cs);
	}
}

void ChemGroup::findcrdatom() {
	std::map<std::string,std::vector<std::vector<std::string>>> aags = AA_Atom::groups_crd();
	for(auto &v:aags.at(resname)) {
		bool isthis{true};
		for(auto &s:v) {
			if(crds.find(s)==crds.end()) {
				isthis=false;
				break;
			}
		}
		if(!isthis) continue;
		crd_atoms=v;
		break;
	}
}

std::vector<LocalFrame> ChemGroup::makelocalframe(const std::map<std::string,XYZ>&cs) {
	std::vector<LocalFrame> lfs;
	for(int i=0;i<crd_atoms.size();i+=3) {
		lfs.push_back(make_localframe(cs.at(crd_atoms[i]),cs.at(crd_atoms[i+1]),cs.at(crd_atoms[i+2])));
	}
	return lfs;
}

std::vector<LocalFrame> ChemGroup::makelocalframe() {
	return makelocalframe(crds);
}

std::vector<LocalFrame> ChemGroup::makelocalframe(int n) {
	return makelocalframe(decoys[n]);
}

std::vector<XYZ> GroupPair::crd4localframe() {
	if(g1->crd_atoms.empty()) g1->findcrdatom();
	std::vector<LocalFrame> lfs1=g1->makelocalframe();
	if(g2->crd_atoms.empty()) g2->findcrdatom();
	std::vector<LocalFrame> lfs2=g2->makelocalframe();
	std::vector<XYZ> cs;
	for(int i=0;i<lfs1.size();i++) {
		for(int j=0;j<lfs2.size();j++) {
			cs.push_back(lfs1[i].global2localcrd(g2->crds.at(g2->crd_atoms[0])));
			cs.push_back(lfs2[j].global2localcrd(g1->crds.at(g1->crd_atoms[0])));
		}
	}
	return cs;
}

std::vector<XYZ> GroupPair::crd4localframe(int n) {
	if(g1->crd_atoms.empty()) g1->findcrdatom();
	std::vector<LocalFrame> lfs1=g1->makelocalframe(n);
	if(g2->crd_atoms.empty()) g2->findcrdatom();
	std::vector<LocalFrame> lfs2=g2->makelocalframe(n);
	std::vector<XYZ> cs;
	for(int i=0;i<lfs1.size();i++) {
		for(int j=0;j<lfs2.size();j++) {
			cs.push_back(lfs1[i].global2localcrd(g2->decoys[n].at(g2->crd_atoms[0])));
			cs.push_back(lfs2[j].global2localcrd(g1->decoys[n].at(g1->crd_atoms[0])));
		}
	}
	return cs;
}




int GroupPair::nbondin2resin1chain(int r1, std::string res1, std::string atom1,
		int r2, std::string res2, std::string atom2) {
	if(r1==r2) return 0;
	if(r2<r1) {
		int tempi=r2;
		r2=r1;
		r1=tempi;
		std::string temp=res2;
		res2=res1;
		res1=temp;
		temp=atom2;
		atom2=atom1;
		atom1=temp;
	}
	std::map<std::string,std::map<std::string,std::vector<int>>> cads=AA_Atom::CADis();
	//std::cout <<res1 <<' ' <<atom1 <<std::endl;
	int c1=cads.at(res1).at(atom1)[1];
	//std::cout <<res2 <<' ' <<atom2 <<std::endl;
	int c2=cads.at(res2).at(atom2)[0];
	return (r2-r1)*3+c1+c2;
}
bool GroupPair::farinteraction(const ChemGroup&cg1, const ChemGroup&cg2, int nbond) {
	std::string rn1=cg1.resname;
	std::string rn2=cg2.resname;
	//std::cout <<rn1 <<' ' <<rn2 <<std::endl;
	if(rn1.substr(0,2)=="MC") rn1="GLY";
	if(rn2.substr(0,2)=="MC") rn2="GLY";
	//std::cout <<"st\t" <<rn1 <<' ' <<rn2 <<std::endl;
	for(int m=0;m<cg1.crds_v.size();m++) {
		for(int n=0;n<cg2.crds_v.size();n++) {
			//std::cout <<rn1 <<' ' <<rn2 <<std::endl;
			//int nb=GroupPair::nbondin2resin1chain(cg1.resseqs[m],rn1,cg1.crds_v[m].first,
			//		cg2.resseqs[n],rn2,cg2.crds_v[n].first);
			int nb=GroupPair::nbondin2resin1chain(cg1.resseq,rn1,cg1.crds_v[m].first,
					cg2.resseq,rn2,cg2.crds_v[n].first);
			if(nb<nbond) return false;
		}
	}
	return true;
}



void Atom3::findatom3() {
	std::map<std::string,std::vector<std::vector<std::string>>> a3=AA_Atom::Atom3();
	for(int i=0;i<chains_.size();i++) {
		/*if(chains_[i].size()!=1) {
			std::pair<std::string,XYZ> p0{"N",chains_[i][0].rds().at("N").crd};
			std::pair<std::string,XYZ> p1{"CA",chains_[i][0].rds().at("CA").crd};
			std::pair<std::string,XYZ> p2{"C",chains_[i][0].rds().at("C").crd};

			std::shared_ptr<ChemGroup> cg = std::shared_ptr < ChemGroup > (new ChemGroup(i,0));
			cg->crds.insert(p0);
			cg->crds.insert(p1);
			cg->crds.insert(p2);

			cg->crds_v.push_back(p0);
			cg->crds_v.push_back(p1);
			cg->crds_v.push_back(p2);

			cg->resseqs.push_back(0);
			cg->resseqs.push_back(0);
			cg->resseqs.push_back(0);

			cg->resname = "MCN";
			cg->makecenter();
			cgs_.push_back(cg);
		}*/
		for(int j=1;j<chains_[i].size();j++) {
			std::pair<std::string,XYZ> p0{"N",chains_[i][j].rds().at("N").crd};
			std::pair<std::string,XYZ> p1{"C",chains_[i][j-1].rds().at("C").crd};
			std::pair<std::string,XYZ> p2{"O",chains_[i][j-1].rds().at("O").crd};

			std::shared_ptr<ChemGroup> cg = std::shared_ptr < ChemGroup > (new ChemGroup(i,j));
			cg->crds.insert(p0);
			cg->crds.insert(p1);
			cg->crds.insert(p2);

			cg->crds_v.push_back(p0);
			cg->crds_v.push_back(p1);
			cg->crds_v.push_back(p2);

			cg->resseqs.push_back(j);
			cg->resseqs.push_back(j-1);
			cg->resseqs.push_back(j-1);

			double t=chains_[i][j].omega(chains_[i][j-1]);
			if(t>-90.0 && t<90.0) {
				if(a3.find("MCCIS")!=a3.end()) {
					cg->resname = "MCCIS";
				}
			} else {
				if(a3.find("MCTRANS")!=a3.end()) {
					cg->resname = "MCTRANS";
				}
			}
			cg->makecenter();
			cgs_.push_back(cg);
		}
		/*if(chains_[i].size()!=1) {
			std::pair<std::string,XYZ> p0{"O",chains_[i].back().rds().at("O").crd};
			std::pair<std::string,XYZ> p1{"C",chains_[i].back().rds().at("C").crd};
			std::pair<std::string,XYZ> p2{"CA",chains_[i].back().rds().at("CA").crd};

			std::shared_ptr<ChemGroup> cg = std::shared_ptr < ChemGroup > (new ChemGroup(i,chains_[i].size()-1));
			cg->crds.insert(p0);
			cg->crds.insert(p1);
			cg->crds.insert(p2);

			cg->crds_v.push_back(p0);
			cg->crds_v.push_back(p1);
			cg->crds_v.push_back(p2);

			cg->resseqs.push_back(chains_[i].size()-1);
			cg->resseqs.push_back(chains_[i].size()-1);
			cg->resseqs.push_back(chains_[i].size()-1);

			cg->resname = "MCC";
			cg->makecenter();
			cgs_.push_back(cg);
		}*/
	}
	std::set<std::string> mcs{"N","C","O"};
	for(int i=0;i<chains_.size();i++) {
		for(int j=0;j<chains_[i].size();j++) {
			std::string rn = chains_[i][j].resname();
			if(a3.find(rn)==a3.end()) {
				std::cout <<"Can Not Find Chemical Group: " <<rn <<std::endl;
				continue;
			}
			for(int k=0;k<a3.at(rn).size();k++) {
				const std::vector<std::string> &an=a3.at(rn)[k];
				bool cninmc{false};
				if(j==0 || j==chains_[i].size()-1) {
					for(const auto &s:an) {
						if(mcs.find(s)==mcs.end()) continue;
						cninmc=true;
						break;
					}
				}
				const auto &rds=chains_[i][j].rds();
				bool isall = (rds.find(an[0])!=rds.end() &&
						rds.find(an[1])!=rds.end() && rds.find(an[2])!=rds.end());
				if(!isall) continue;
				std::shared_ptr<ChemGroup> cg = std::shared_ptr < ChemGroup > (new ChemGroup(i,j));
				std::pair<std::string,XYZ> p0{an[0],rds.at(an[0]).crd};
				std::pair<std::string,XYZ> p1{an[1],rds.at(an[1]).crd};
				std::pair<std::string,XYZ> p2{an[2],rds.at(an[2]).crd};
				cg->crds.insert(p0);
				cg->crds.insert(p1);
				cg->crds.insert(p2);
				cg->crds_v.push_back(p0);
				cg->crds_v.push_back(p1);
				cg->crds_v.push_back(p2);
				cg->resseqs.push_back(j);
				cg->resseqs.push_back(j);
				cg->resseqs.push_back(j);
				cg->resname = rn;
				cg->makecenter();
				cg->aa20cn=cninmc;
				cgs_.push_back(cg);
			}
		}
	}
}

void Atom3::cover2NM() {
	for(auto &cg:cgs_) {
		for(auto &p:cg->crds) {
			p.second.x_ *= 0.1;
			p.second.y_ *= 0.1;
			p.second.z_ *= 0.1;
		}
		for(auto &p:cg->crds_v) {
			p.second.x_ *= 0.1;
			p.second.y_ *= 0.1;
			p.second.z_ *= 0.1;
		}
	}
}

void Atom3::findatom3pair(bool removeclash) {
	//if(cgs_.empty()) findatom3();
	double md2=MaxDistance_*MaxDistance_;
	int nmainchain=5;
	for(int i=0;i<cgs_.size();i++) {
		if(cgs_[i]->aa20cn) continue;
		for(int j=i+1;j<cgs_.size();j++) {
			if(cgs_[j]->aa20cn) continue;
			//2 groups can not be in the same AA
			if(cgs_[i]->chainseq==cgs_[j]->chainseq&&cgs_[i]->resseq==cgs_[j]->resseq) continue;
			//make group pair
			double dis2=(cgs_[i]->center-cgs_[j]->center).squarednorm();
			if(dis2>md2) continue;
			double dis=sqrt(dis2);
			GroupPair gp(cgs_[i],cgs_[j]);
			if(removeclash) {
				if(dis<5.0 && gp.clash()) continue;
			}
			gp.dis=dis;
			//number of bond in 2 group must be less than MinAtomInterval_
			int nbond=1000;
			if(cgs_[i]->chainseq==cgs_[j]->chainseq) {
				int nm=fabs(cgs_[i]->resseq==cgs_[j]->resseq);
				if(nm<nmainchain) {
					std::string rni=cgs_[i]->resname;
					std::string rnj=cgs_[j]->resname;
					//if(cgs_[i]->resseqs[0]!=cgs_[i]->resseqs[2]) rni="GLY";
					//if(cgs_[j]->resseqs[0]!=cgs_[j]->resseqs[2]) rnj="GLY";
					if(cgs_[i]->resname.substr(0,2)=="MC") rni="GLY";
					if(cgs_[j]->resname.substr(0,2)=="MC") rnj="GLY";
					for(int m=0;m<cgs_[i]->crds_v.size();m++) {
						for(int n=0;n<cgs_[j]->crds_v.size();n++) {
							int nb=GroupPair::nbondin2resin1chain(cgs_[i]->resseqs[m],rni,
									cgs_[i]->crds_v[m].first,cgs_[j]->resseqs[n],
									rnj,cgs_[j]->crds_v[n].first);
							if(nb<nbond) nbond=nb;
							//if(nbond<MinAtomInterval_) break;
						}
						//if(nbond<MinAtomInterval_) break;
					}
					//if(nbond<MinAtomInterval_) continue;
				}
			}
			gp.nbond=nbond;
			gp.rerank();
			gp.rerank_samegroup(0,2);

			if(nbond<MinAtomInterval_) continue;
			/*if(nbond<MinAtomInterval_) {
				//if(nbond>1) gps_inbond.push_back(gp);
				continue;
			} else if(cgs_[i]->inmainchain() && cgs_[j]->inmainchain()) {
				gps_mcmc.push_back(gp);
			} else {
				gps_sc.push_back(gp);
			}*/
			gps_.push_back(gp);
		}
	}
}

std::vector<GroupPair> Atom3::gpsinatomdis(double dis) {
	std::vector<GroupPair> gp;
	for(GroupPair &p:gps_) {
		if(p.minatomdis()<dis) gp.push_back(p);
	}
	return gp;
}











double TREEf::GetPoint(double query) {
	static std::map<double,double,Comp_Double> mp;
#pragma omp critical
	{
	if(mp.empty()) {
		double half = BinSize / 2.0;
		for(double d=half;d<TreeMaxVal;d+=BinSize) {
			mp.insert({d,d});
		}
		for(double d=-half;d>-TreeMaxVal;d-=BinSize) {
			mp.insert({d,d});
		}
	}
	}
	if(query>TreeMaxVal || query<-TreeMaxVal) {
		std::cout <<"Value is not in Range!!!" <<std::endl;
		exit(1);
	}
	if(mp.find(query)==mp.end()) {
		std::cout <<"Query Not Found     " <<query <<std::endl;
		if(query<0.0) query += 0.00000001;
		else query -= 0.00000001;
	}
	return mp.at(query);
}
double TREEf::GetPoint2(double query) {
	static std::map<double,double,Comp_Double> mp;
#pragma omp critical
	{
	if(mp.empty()) {
		double half = BinSize2 / 2.0;
		for(double d=half;d<TreeMaxVal;d+=BinSize2) {
			mp.insert({d,d});
		}
		for(double d=-half;d>-TreeMaxVal;d-=BinSize2) {
			mp.insert({d,d});
		}
	}
	}
	if(query>TreeMaxVal || query<-TreeMaxVal) {
		std::cout <<"Value is not in Range!!!" <<std::endl;
		exit(1);
	}
	if(mp.find(query)==mp.end()) {
		std::cout <<"Query Not Found     " <<query <<std::endl;
		if(query<0.0) query += 0.00000001;
		else query -= 0.00000001;
	}
	return mp.at(query);
}
std::vector<double> TREEf::MakeKey(const std::vector<double> &query) {
	std::vector<double> key;
	for(double d:query) {
		key.push_back(GetPoint(d));
	}
	return key;
}
std::vector<double> TREEf::MakeKey2(const std::vector<double> &query) {
	std::vector<double> key;
	for(double d:query) {
		key.push_back(GetPoint2(d));
	}
	return key;
}
void TREEf::BuildTree() {
	for(long i=0;i<dts.size();i++) {
		std::vector<double> key=MakeKey(dts[i]);
		std::vector<double> key2=MakeKey2(dts[i]);
		auto iter2=tree.find(key2);
		if(iter2==tree.end()) {
			BRANCH bh;
			bh.insert({key,{i}});
			tree.insert({key2,bh});
		} else {
			auto iter=iter2->second.find(key);
			if(iter==iter2->second.end()) {
				iter2->second.insert({key,{i}});
			} else {
				iter->second.insert(i);
			}
		}
	}
}
void TREEf::BuildTree_omp() {
#pragma omp parallel for
	for(long i=0;i<dts.size();i++) {
		std::vector<double> key=MakeKey(dts[i]);
		std::vector<double> key2=MakeKey2(dts[i]);
#pragma omp critical
		{
			auto iter2=tree.find(key2);
			if(iter2==tree.end()) {
				BRANCH bh;
				bh.insert({key,{i}});
				tree.insert({key2,bh});
			} else {
				auto iter=iter2->second.find(key);
				if(iter==iter2->second.end()) {
					iter2->second.insert({key,{i}});
				} else {
					iter->second.insert(i);
				}
			}
		}
	}
}
std::vector<std::vector<double>> TREEf::GetBins(
		const std::vector<double> &query, const std::vector<double> &ranges) {
	std::vector<std::vector<double>> vals;
	for(int i=0;i<query.size();i++) {
		double low = query[i]-ranges[i];
		double high = query[i]+ranges[i];
		low = GetPoint(low);
		high = GetPoint(high);
		std::vector<double> vd;
		for(double d=low;d<high+0.000001;d+=BinSize) {
			vd.push_back(d);
		}
		//assert(vd.size()>0);
		vals.push_back(vd);
	}

	std::vector<std::vector<double>> bins(1);
	for(int i=0;i<vals.size();i++) {
		bins[0].push_back(vals[i][0]);
	}
	for(int i=0;i<vals.size();i++) {
		int nz = bins.size();
		for(int j=1;j<vals[i].size();j++) {
			for(int k=0;k<nz;k++) {
				std::vector<double> vd = bins[k];
				vd[i] = vals[i][j];
				bins.push_back(vd);
			}
		}
	}

	//int nsize=1;
	//for(int i=0;i<vals.size();i++) nsize *= vals[i].size();

	return bins;
}
std::set<long> TREEf::findneighbor(const std::vector<double> &query, const std::vector<double>&ranges) {
	std::set<long> si;
	std::vector<std::vector<double>> stages = GetBins(query,ranges);
	/*for(auto &s0:stages) {
		auto iter2=tree.find(s0);
		if(iter2==tree.end()) continue;
		for(auto i:iter->second) {
			si.insert(i);
		}
	}*/
	return si;
}
double TREEf::f10_old(double val, double rad) {
	if(val>rad) return 0.0;
	double hr = rad / 2.0;
	if(val<hr) return 1.0;
	double r34 = rad * 3.0 / 4.0;
	double a = 8.0 / rad / rad;
	double delta = r34 - val;
	double v = a*delta*delta;
	if(val<r34) v=1.0-v;
	return v;
}
double TREEf::f10(double val, double rad) {
	if(val>rad) return 0.0;
	double r1 = 0.1;
	if(val<r1) return 1.0;
	double r2 = (rad + r1) / 2.0;
	double r21 = r2 - r1;
	double a = 0.5 / r21 / r21;
	if(val<r2) {
		double delta = val-r1;
		return 1.0-a*delta*delta;
	}
	double delta = rad - val;
	return a*delta*delta;
}
double TREEf::score1(std::set<long> &si, const std::vector<double> &qy, double rad) {
	double hr = rad/2.0;
	double sc=0.0;
	for(long n:si) {
		double s1=1.0;
		for(int i=0;i<Ndim;i++) {
			double delta = fabs(qy[i]-dts[n][i]);
			if(delta>rad) {
				s1=0.0;
				break;
			} else s1 *= f10(delta,rad);
		}

		double w = 1.0;
		if(weights.size()==dts.size()) w= weights[n];
		sc += s1 * w;
	}
	return sc;
}
double TREEf::score(const std::vector<double> &query, const std::vector<double>&ranges, double rad) {
	std::set<long> sl=findneighbor(query,ranges);
	return score1(sl,query,rad);
}


bool NSPallatom::samechemgroup(const ChemGroup&c1, const ChemGroup&c2) {
	if(c1.resname!=c2.resname) return false;
	if(c1.chainseq!=c2.chainseq) return false;
	if(c1.resseq!=c2.resseq) return false;
	if(c1.crds_v.size()!=c2.crds_v.size()) return false;
	for(int i=0;i<c1.crds_v.size();i++) {
		XYZ c=c1.crds_v[i].second-c2.crds_v[i].second;
		if(c.squarednorm()>0.00001) return false;
	}
	return true;
}











