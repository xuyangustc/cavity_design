/*
 * pdbreader_xy.cpp
 *
 *  Created on: Nov 20, 2018
 *      Author: xuyang
 */
#include "allatom/pdbreader_xy.h"
#include "sd/sidechainff.h"
#include "allatom/basic_info.h"
#include "fullsite/fullsite.h"
using namespace NSPgeometry;
using namespace NSPallatom;
using namespace NSPproteinrep;

XYZ Residue::cb() {
	if(records_.find("CB")==records_.end()) {
		NSPgeometry::XYZ c=records_.at("C").crd;
		NSPgeometry::XYZ ca=records_.at("CA").crd;
		NSPgeometry::XYZ n=records_.at("N").crd;
		double theta=109.5*3.14159265/180.0;
		double t=120*3.14159265/180.0;
		double b0=1.5;
		return NSPgeometry::InternaltoXYZ(ca,c,n,b0,theta,t);
	} else return records_.at("CB").crd;
}

std::vector<std::vector<double>> NSPallatom::getphipsi(const std::vector<Residue>&chain) {
	std::vector<std::vector<double>> dihs;
	for(int i=1;i<chain.size()-1;i++) {
		double phi = chain[i].phi(chain[i-1].rds().at("C").crd);
		double psi = chain[i].psi(chain[i+1].rds().at("N").crd);
		dihs.push_back({phi,psi});
	}
	return dihs;
}

void NSPallatom::crd2h2opdb(std::ostream &os,std::vector<NSPgeometry::XYZ>&crds, int ast, int rst) {
	NSPproteinrep::PdbRecord pr;
	pr.label = "HETATM";
	pr.namesymbol = "O";
	pr.atomname = "O";
	pr.residuename = "HOH";
	pr.chainid = 'A';
	pr.elementname[1] = 'O';

	for(int i=0;i<crds.size();i++) {
		auto p = pr;
		p.atomid = ast++;
		p.residueid = rst++;
		p.x = crds[i].x_;
		p.y = crds[i].y_;
		p.z = crds[i].z_;
		os <<p.toString() <<std::endl;
	}
}

std::string NSPallatom::getchainid() {
	std::string chid;
	for(int i=65;i<=90;i++) chid.push_back(i);
	for(int i=97;i<=122;i++) chid.push_back(i);
	for(int i=48;i<=57;i++) chid.push_back(i);
	for(int i=33;i<=47;i++) chid.push_back(i);
	for(int i=58;i<=64;i++) chid.push_back(i);
	for(int i=91;i<=96;i++) chid.push_back(i);
	for(int i=123;i<=126;i++) chid.push_back(i);
	return chid;
}

int NSPallatom::residue2pdb(std::ostream &os, const std::vector<std::vector<Residue>>& outpep, bool onlyheavyatom) {
	int aid=0;
	std::string chid=getchainid();
	char cid='A';
	for(int i=0;i<outpep.size();i++) {
		//if(i<26) cid = 'A' + i;
		//else if(i<52) cid = 'a' + i - 26;
		//else cid = '0' + i - 52;
		cid = chid[i];
		for(int j=0;j<outpep[i].size();j++) {
			std::vector<NSPproteinrep::PdbRecord> prs=outpep[i][j].getpdbrecords(onlyheavyatom);
			for(NSPproteinrep::PdbRecord &p:prs) {
				p.residueid = j;
				p.atomid = aid++;
				p.chainid = cid;
				if(aid>99999) {
					aid = -9999;
				}
				os <<p.toString() <<std::endl;
			}
		}
	}
	return aid;
}

void NSPallatom::residue2pdb(std::string outfile, const std::vector<std::vector<Residue>>& outpep, bool onlyheavyatom) {
	std::ofstream ofs(outfile);
	residue2pdb(ofs,outpep,onlyheavyatom);
	ofs.close();
}

void NSPallatom::residue2pdb(std::string outfile, const std::vector<std::vector<Residue>>& outpep,
		std::vector<std::vector<int>> toprint) {
	std::vector<std::vector<std::vector<int>>> lbs;
	for(auto &vi:toprint) {
		bool added{false};
		for(auto &vvi:lbs) {
			std::vector<int> &vii=vvi.back();
			if(vii[0]==vi[0] && vii[1]==vi[1]-1) {
				vvi.push_back(vi);
				added=true;
				break;
			}
		}
		if(!added) {
			lbs.push_back(std::vector<std::vector<int>>());
			lbs.back().push_back(vi);
		}
		bool isall{false};
		while(!isall) {
			isall = true;
			for(int i=0;i<lbs.size();i++) {
				for(int j=i+1;j<lbs.size();j++) {
					if(lbs[i][0][0]==lbs[j].back()[0] && lbs[i][0][1]-1==lbs[j].back()[1]) {
						for(auto &vii:lbs[i]) lbs[j].push_back(vii);
						lbs[i].clear();
						isall = false;
					} else if(lbs[j][0][0]==lbs[i].back()[0] && lbs[j][0][1]-1==lbs[i].back()[1]) {
						for(auto &vii:lbs[j]) lbs[i].push_back(vii);
						lbs[j].clear();
						isall = false;
					}
					if(!isall) {
						std::vector<std::vector<std::vector<int>>> lbs1;
						for(auto &vvi:lbs) {
							if(vvi.empty()) continue;
							lbs1.push_back(vvi);
						}
						lbs = lbs1;
						break;
					}
				}
				if(!isall) break;
			}
		}
	}
	std::vector<std::vector<Residue>> ress;
	for(auto &vvi:lbs) {
		ress.push_back(std::vector<Residue>());
		for(auto &vi:vvi) {
			ress.back().push_back(outpep[vi[0]][vi[1]]);
		}
	}
	residue2pdb(outfile,ress);
}


void Residue::init(const std::map<std::string,PdbRecord> &prs, int cid, int rid) {
	chainid_ = prs.at("CA").chainid;
	resid_ = prs.at("CA").residueid;
	resname_ = prs.at("CA").residuename;
	insertionid_ = prs.at("CA").insertionid;
	chainid_new = cid;
	resid_new = rid;
	for(auto &p:prs) records_.insert({p.first,AtomDataInPDB(p.second)});
	std::set<std::string> atoms;//std::cout <<'\t' <<resname_ <<std::endl;
	auto aa20=AA_Atom::aa31();
	if(aa20.find(resname_)!=aa20.end()) {//std::cout <<resname_ <<std::endl;
		auto &vsc=NSPsd::VSCType::getVSCType(resname_);
		std::vector<std::string> scatoms=vsc.atomnames;
		for(std::string & s:scatoms) atoms.insert(s);
	}
	atoms.insert("N");
	atoms.insert("CA");
	atoms.insert("C");
	atoms.insert("O");
	for(auto &a:atoms) {
		if(records_.find(a)==records_.end()) atoms_lack.insert(a);
	}
	argnh1nh2();
}

void Peptide::init(const std::vector<std::map<std::string,PdbRecord>>&prs, bool removeMCCAclash) {
	int cid=0;
	int rid=0;
	chains_new.resize(cid+1);
	XYZ ccrd(prs[0].at("N").x,prs[0].at("N").y,prs[0].at("N").z);
	XYZ ncrd;
	char cid_prev = prs[0].at("N").chainid;
	char cid_now;
	double cutoff=2.0*2.0;
	for(const auto &ps:prs) {
		ncrd=XYZ(ps.at("N").x,ps.at("N").y,ps.at("N").z);
		cid_now = ps.at("N").chainid;
		if((ncrd-ccrd).squarednorm()>cutoff || cid_now!=cid_prev) {
			cid++;
			rid=0;
			chains_new.resize(cid+1);
		}
		Residue restemp(ps,cid,rid++);
		if(removeMCCAclash) {
			if(chains_new.back().size()!=0) {
				XYZ ca1=chains_new.back().back().rds().at("CA").crd;
				XYZ ca2=restemp.rds().at("CA").crd;
				double d2=(ca1-ca2).squarednorm();
				if(d2<4.0) {
					rid--;
					continue;
				}
			}
		}
		chains_new.back().push_back(restemp);
		ccrd=chains_new.back().back().rds().at("C").crd;
		cid_prev=cid_now;
	}
	char cc=chains_new[0][0].chid();
	chains_.push_back(std::vector<Residue>());
	for(const auto &ch:chains_new) {
		char cc1=ch[0].chid();
		if(cc1!=cc) {
			cc=cc1;
			chains_.push_back(std::vector<Residue>());
		}
		if(removeMCCAclash) {
			for(const Residue &r:ch){
				if(chains_.back().size()!=0) {
					XYZ ca1=chains_.back().back().rds().at("CA").crd;
					XYZ ca2=r.rds().at("CA").crd;
					double d2=(ca1-ca2).squarednorm();
					if(d2<4.0) {
						continue;
					}
				}
			}
		} else {
			for(const Residue &r:ch) chains_.back().push_back(r);
		}
	}
}

std::vector<std::vector<Residue>> Peptide::continuechain(const std::vector<Residue> &ress) {
	std::vector<std::vector<Residue>> cch;
	std::vector<Residue> temp{ress[0]};
	double cutoff=2.0*2.0;
	for(int i=1;i<ress.size();i++) {
		XYZ ccrd=ress[i-1].rds().at("C").crd;
		XYZ ncrd=ress[i].rds().at("N").crd;
		if((ccrd-ncrd).squarednorm()>cutoff) {
			cch.push_back(temp);
			temp.clear();
		}
		temp.push_back(ress[i]);
	}
	cch.push_back(temp);
	return cch;
}

void PdbReader_xy::init(std::vector<std::vector<PdbRecord>>&prs, bool removeMCCAclash) {
	std::vector<std::map<std::string,PdbRecord>> ress;
	for(const auto &p:prs) {
		if(p[0].residuename=="HOH") {
			for(const auto &p1:p) h2os_.push_back(p1);
			continue;
		}
		bool n,ca,c,o;
		n=ca=c=o=false;
		for(const PdbRecord &p1:p) {
			if(p1.atomname=="N") n=true;
			else if(p1.atomname=="CA") ca=true;
			else if(p1.atomname=="C") c=true;
			else if(p1.atomname=="O") o=true;
		}
		if(!n||!ca||!c||!o) {
			if(p[0].label=="HETATM") ligands_.push_back(p);
			continue;
		}
		std::map<std::string,PdbRecord> rs;
		for(const auto &p1:p) {
			rs.insert({p1.atomname,p1});
		}
		ress.push_back(rs);
	}
	if(ress.empty()) return ;
	Peptide pep1;
	pep1.init(ress,removeMCCAclash);
	peps_.push_back(pep1);
}
void PdbReader_xy::init(std::string pdbfile, bool removeMCCAclash) {
	std::vector<std::vector<std::string>> lines;
	std::ifstream ifs(pdbfile);
	std::string line2;
	std::vector<std::string> line1;
	while(std::getline(ifs,line2)) {
		std::string subl = line2.substr(0,6);
		if(subl=="ENDMDL") {
			if(!line1.empty()) lines.push_back(line1);
			line1.clear();
			continue;
		}
		if(subl!="ATOM  " && subl!="HETATM") continue;
		line1.push_back(line2);
	}
	if(!line1.empty()) lines.push_back(line1);
	ifs.close();
	for(int i=0;i<lines.size();i++) {
		std::vector<std::vector<PdbRecord>> prs;
		std::vector<PdbRecord> ps;//int nl=0;
		for(std::string &line:lines[i]) {
			if(ps.empty()) {
				ps.push_back(PdbRecord(line));
				continue;
			}
			PdbRecord p(line);
			PdbRecord p2=ps.back();
			if(p2.chainid!=p.chainid || p2.residueid!=p.residueid || p2.insertionid!=p.insertionid) {
				prs.push_back(ps);
				ps.clear();
			}
			ps.push_back(p);
		}
		if(!ps.empty()) prs.push_back(ps);
		if(!prs.empty()) init(prs,removeMCCAclash);
	}
}

Residue::Residue(const NSPproteinrep::BackBoneSite &bbs) {
	chainid_ = bbs.chainid;
	resid_ = bbs.resid;
	resname_ = bbs.resname;
	resid_new = bbs.resseq;
	records_.insert({"N",AtomDataInPDB('N')});
	records_.at("N").crd.x_=bbs.data_[4];
	records_.at("N").crd.y_=bbs.data_[5];
	records_.at("N").crd.z_=bbs.data_[6];
	records_.insert({"CA",AtomDataInPDB('C')});
	records_.at("CA").crd.x_=bbs.data_[7];
	records_.at("CA").crd.y_=bbs.data_[8];
	records_.at("CA").crd.z_=bbs.data_[9];
	records_.insert({"C",AtomDataInPDB('C')});
	records_.at("C").crd.x_=bbs.data_[10];
	records_.at("C").crd.y_=bbs.data_[11];
	records_.at("C").crd.z_=bbs.data_[12];
	records_.insert({"O",AtomDataInPDB('O')});
	records_.at("O").crd.x_=bbs.data_[13];
	records_.at("O").crd.y_=bbs.data_[14];
	records_.at("O").crd.z_=bbs.data_[15];
}

Residue::Residue(const NSPproteinrep::FullSite &fs) {
	chainid_ = fs.chainid();
	resid_ = fs.resid();
	resname_ = fs.resname();
	resid_new = fs.resseq();
	std::map<std::string,XYZ> cs=fs.getcrds();
	for(auto &c:cs) {
		records_.insert({c.first,AtomDataInPDB(c.first[0])});
		records_.at(c.first).crd=c.second;
	}
}

Residue::Residue(const NSPproteinrep::BackBoneSite &bbs, std::string aaname) {
	chainid_ = bbs.chainid;
	resid_ = bbs.resid;
	resname_ = aaname;
	resid_new = bbs.resseq;
	BackBoneSite bs=bbs;
	bs.resname=aaname;
	FullSite fs=make_fullsite(bs);
	std::map<std::string,NSPgeometry::XYZ> cs=fs.getcrds();
	for(auto &c:cs) {
		records_.insert({c.first,AtomDataInPDB(c.first[0])});
		records_.at(c.first).crd=c.second;
	}
}

void Residue::changeaa(std::string newname) {
	resname_ = newname;
	BackBoneSite bs=getbackbonesite();
	bs.resname=newname;
	records_.clear();
	FullSite fs=make_fullsite(bs);
	std::map<std::string,NSPgeometry::XYZ> cs=fs.getcrds();
	for(auto &c:cs) {
		records_.insert({c.first,AtomDataInPDB()});
		records_.at(c.first).crd=c.second;
		records_.at(c.first).element.push_back(' ');
		records_.at(c.first).element.push_back(c.first[0]);
	}
	argnh1nh2();
}

NSPproteinrep::FullSite Residue::getfullsite() const {
	NSPproteinrep::FullSite fs;
	std::map<std::string,XYZ> cs;
	for(const auto &c:records_) cs.insert({c.first,c.second.crd});
	fs.changecrd(cs);
	fs.resname()=resname_;
	fs.chainid()=chainid_;
	fs.resid()=resid_;
	fs.resseq()=resid_new;
	fs.insertionid()=insertionid_;
	return fs;
}

NSPproteinrep::BackBoneSite Residue::getbackbonesite() const {
	NSPproteinrep::BackBoneSite bbs;
	bbs.chainid = chainid_;
	bbs.resid = resid_;
	bbs.resname = resname_;
	bbs.resseq = resid_new;
	bbs.data_[4] = records_.at("N").crd.x_;
	bbs.data_[5] = records_.at("N").crd.y_;
	bbs.data_[6] = records_.at("N").crd.z_;
	bbs.data_[7] = records_.at("CA").crd.x_;
	bbs.data_[8] = records_.at("CA").crd.y_;
	bbs.data_[9] = records_.at("CA").crd.z_;
	bbs.data_[10] = records_.at("C").crd.x_;
	bbs.data_[11] = records_.at("C").crd.y_;
	bbs.data_[12] = records_.at("C").crd.z_;
	bbs.data_[13] = records_.at("O").crd.x_;
	bbs.data_[14] = records_.at("O").crd.y_;
	bbs.data_[15] = records_.at("O").crd.z_;
	return bbs;
}

std::vector<NSPproteinrep::PdbRecord> Residue::getpdbrecords(bool onlyheavyatom) const {
	std::set<std::string> al=AA_Atom::mainchain();
	auto asc = AA_Atom::sidechain3();
	if(asc.find(resname_)!=asc.end()) {
		std::set<std::string> al1=AA_Atom::sidechain3().at(resname_);
		for(auto &a:al1) al.insert(a);
	}
	std::vector<PdbRecord> prs;
	char cid{'A'};
	for(int i=1;i<chainid_new;i++) cid = cid + 1;
	for(auto &r:records_) {
		if(onlyheavyatom) {
			if(al.find(r.first)==al.end()) continue;
		}
		PdbRecord p;
		p.label = "ATOM";
		p.namesymbol = r.first.substr(0,1);
		p.namemodifier = r.first.substr(1);
		p.atomname = r.first;
		p.residuename = resname_;
		p.chainid = cid;
		p.occupation = r.second.occupancy;
		p.bfactor = r.second.tempfac;
		p.elementname[0] = r.second.element[0];
		p.elementname[1] = r.second.element[1];
		p.x = r.second.crd.x_;
		p.y = r.second.crd.y_;
		p.z = r.second.crd.z_;
		prs.push_back(p);
	}
	return prs;
}

std::vector<NSPproteinrep::PdbRecord> Residue::getpdbrecords(int &atomidst, int residst) const {
	std::vector<PdbRecord> prs;
	char cid{'A'};
	//for(int i=1;i<chainid_new;i++) cid = cid + 1;
	for(auto &r:records_) {
		PdbRecord p;
		p.label = "ATOM";
		p.atomid = atomidst++;
		p.namesymbol = r.first.substr(0,1);
		p.namemodifier = r.first.substr(1);
		p.atomname = r.first;
		p.residuename = resname_;
		p.chainid = cid;
		p.residueid = residst;
		p.occupation = r.second.occupancy;
		p.bfactor = r.second.tempfac;
		p.elementname[0] = r.second.element[0];
		p.elementname[1] = r.second.element[1];
		p.x = r.second.crd.x_;
		p.y = r.second.crd.y_;
		p.z = r.second.crd.z_;
		prs.push_back(p);
	}
	return prs;
}

double Residue::phi(const XYZ &prevc) const {
	double t = NSPgeometry::torsion(prevc, records_.at("N").crd,
			records_.at("CA").crd, records_.at("C").crd);
	double rad = 180.0 / 3.14159265358979323846;
	t *= rad;
	while (t > 180.0)
		t -= 360.0;
	while (t < -180.0)
		t += 360.0;
	return t;
}

double Residue::psi(const XYZ &behn) const {
	double t = NSPgeometry::torsion(records_.at("N").crd,
			records_.at("CA").crd, records_.at("C").crd,behn);
	double rad = 180.0 / 3.14159265358979323846;
	t *= rad;
	while (t > 180.0)
		t -= 360.0;
	while (t < -180.0)
		t += 360.0;
	return t;
}
/*
void PdbReader_xy::printpdb_new(std::ostream &os) {
	const std::vector<std::vector<Residue>> & ress=pep_.chs_new();
	int aid=0;
	for(int i=0;i<ress.size();i++) {
		for(int j=0;j<ress[i].size();j++) {
			std::vector<PdbRecord> prs=ress[i][j].getpdbrecords();
			for(PdbRecord &p:prs) {
				p.residueid = j;
				p.atomid = aid++;
				os <<p.toString() <<std::endl;
			}
		}
	}
}*/

void Residue::mse2met() {
	if(resname_!="MSE") return;
	resname_="MET";
	//if(records_.find("SE")==records_.end()) return;
	if(records_.find("SE")!=records_.end()) {
		records_.insert({"SD",records_.at("SE")});
		records_.at("SD").element="S";
		records_.erase("SE");
	}
	std::set<std::string> atoms;//std::cout <<'\t' <<resname_ <<std::endl;
	auto aa20=AA_Atom::aa31();
	auto &vsc=NSPsd::VSCType::getVSCType("MET");
	std::vector<std::string> scatoms=vsc.atomnames;
	for(std::string & s:scatoms) atoms.insert(s);
	atoms.insert("N");
	atoms.insert("CA");
	atoms.insert("C");
	atoms.insert("O");
	for(auto &a:atoms) {
		if(records_.find(a)==records_.end()) atoms_lack.insert(a);
	}
}

void Peptide::mse2met() {
	for(auto &ch:chains_) {
		for(Residue &r:ch) r.mse2met();
	}
	for(auto &ch:chains_new) {
		for(Residue &r:ch) r.mse2met();
	}
}


void Residue::argnh1nh2() {
	if(resname_!="ARG") return;
	if(records_.find("CD") == records_.end()) return;
	if(records_.find("NH1") == records_.end()) return;
	if(records_.find("NH2") == records_.end()) return;

	XYZ c0 = records_.at("CD").crd;
	XYZ c1 = records_.at("NH1").crd;
	XYZ c2 = records_.at("NH2").crd;

	double d21 = (c1-c0).squarednorm();
	double d22 = (c2-c0).squarednorm();

	if(d21<d22) return;

	records_.at("NH1").crd = c2;
	records_.at("NH2").crd = c1;

}

double Residue::omega(const Residue &prevRes) {
	double t = NSPgeometry::torsion(prevRes.rds().at("CA").crd, prevRes.rds().at("C").crd,
			records_.at("N").crd, records_.at("CA").crd);
	double rad = 180.0 / 3.14159265358979323846;
	t *= rad;
	while (t > 180.0)
		t -= 360.0;
	while (t < -180.0)
		t += 360.0;
	return t;
}








