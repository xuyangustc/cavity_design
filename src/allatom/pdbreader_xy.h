/*
 * pdbreader_xy.h
 *
 *  Created on: Nov 18, 2018
 *      Author: xuyang
 */

#ifndef ALLATOM_PDBREADER_XY_H_
#define ALLATOM_PDBREADER_XY_H_

#include "geometry/xyz.h"
#include "proteinrep/pdbrecord.h"
#include "backbone/backbonesite.h"
#include "fullsite/fullsite.h"
#include <map>
#include <set>
#include <vector>
using namespace NSPgeometry;

namespace NSPallatom {

struct AtomDataInPDB {
	XYZ crd;
	double occupancy{0.0};
	double tempfac{0.0};
	std::string element;
	std::string charge;
	AtomDataInPDB() {;}
	AtomDataInPDB(const NSPproteinrep::PdbRecord &p):occupancy(p.occupation), tempfac(p.bfactor) {
		crd.x_ = p.x;
		crd.y_ = p.y;
		crd.z_ = p.z;
		element.push_back(p.elementname[0]);
		element.push_back(p.elementname[1]);
	}
	AtomDataInPDB(const XYZ &c, char e) {
		element.push_back(' ');
		element.push_back(e);
		crd.x_ = c.x_;
		crd.y_ = c.y_;
		crd.z_ = c.z_;
	}
	AtomDataInPDB(char e) {
		element.push_back(' ');
		element.push_back(e);
	}
};

class Residue {
public:
	Residue() {;}
	void init(const std::map<std::string,NSPproteinrep::PdbRecord> &prs, int cid, int rid);
	Residue(const std::map<std::string,NSPproteinrep::PdbRecord> &prs, int cid, int rid) {
		init(prs,cid,rid);
	}
	Residue(const NSPproteinrep::BackBoneSite &bbs);
	Residue(const NSPproteinrep::BackBoneSite &bbs, std::string aaname);
	Residue(const NSPproteinrep::FullSite &fs);
	//char& chid() {return chainid_;}
	const char& chid() const {return chainid_;}
	const int & rid() const {return resid_;}
	const int & chid_new() const {return chainid_new;}
	const int & rid_new() const {return resid_new;}
	const std::string & resname() const {return resname_;}
	const std::map<std::string,AtomDataInPDB> & rds() const {return records_;}
	const std::set<std::string> & lack() const {return atoms_lack;}
	NSPproteinrep::BackBoneSite getbackbonesite() const;
	NSPproteinrep::FullSite getfullsite() const;
	std::vector<NSPproteinrep::PdbRecord> getpdbrecords(bool onlyheavyatom=false) const;
	std::vector<NSPproteinrep::PdbRecord> getpdbrecords(int &atomidst, int residst) const;
	char& chid() {return chainid_;}
	int& rid() {return resid_;};
	int& chid_new() {return chainid_new;}
	int& rid_new() {return resid_new;}
	std::map<std::string,AtomDataInPDB> & rds() {return records_;}
	std::string &resname() {return resname_;}
	double phi(const XYZ &prevc) const;
	double psi(const XYZ &behn) const;
	void changeaa(std::string newname);
	void mse2met();
	void argnh1nh2();
	XYZ cb();
	double omega(const Residue &prevRes);//return dihedral CA-C-N-CA, judge peptide-plane is cis or trans
	char& ssid() {return ss_;}
	char& insertionid() {return insertionid_;}
private:
	char chainid_{'X'};
	int resid_{-10000};
	char insertionid_{' '};
	std::string resname_;
	std::map<std::string,AtomDataInPDB> records_;
	std::set<std::string> atoms_lack;
	int chainid_new{-10000};//start from 0
	int resid_new{-10000};//start from 0
	char ss_{'-'};
};

//size of return value is 2 less than chain
std::vector<std::vector<double>> getphipsi(const std::vector<Residue>&chain);

void crd2h2opdb(std::ostream &os,std::vector<NSPgeometry::XYZ>&crds, int ast, int rst);
std::string getchainid();
int residue2pdb(std::ostream &os, const std::vector<std::vector<Residue>>& outpep, bool onlyheavyatom=false);
void residue2pdb(std::string outfile, const std::vector<std::vector<Residue>>& outpep, bool onlyheavyatom=false);
void residue2pdb(std::string outfile, const std::vector<std::vector<Residue>>& outpep,
		std::vector<std::vector<int>> toprint); //chainid, resid

class Peptide {
public:
	static std::vector<std::vector<Residue>> continuechain(const std::vector<Residue> &ress);
	void init(const std::vector<std::map<std::string,NSPproteinrep::PdbRecord>>&prs, bool removeMCCAclash);
	const std::vector<std::vector<Residue>> & chs() const {return chains_;}
	const std::vector<std::vector<Residue>> & chs_new() const {return chains_new;}
	void mse2met();
private:
	std::vector<std::vector<Residue>> chains_;
	std::vector<std::vector<Residue>> chains_new;
};

class PdbReader_xy {
public:
	void init(std::string pdbfile, bool removeMCCAclash=false);//only residue with all of "C", "CA", "N", "O" can add to Peptide
	void init(std::vector<std::vector<NSPproteinrep::PdbRecord>>&prs, bool removeMCCAclash);
	//Peptide & pep() {return pep_;}
	//void printpdb_new(std::ostream &os);
	void mse2met() {
		for(Peptide&p:peps_) p.mse2met();
	}
	std::vector<NSPproteinrep::PdbRecord> waters() {
		return h2os_;
	}
	std::vector<std::vector<NSPproteinrep::PdbRecord>> ligands() {
		return ligands_;
	}
	std::vector<std::vector<Residue>> chains() {
		std::vector<std::vector<Residue>> rts;
		for(Peptide &p:peps_) for(const auto &ch:p.chs()) rts.push_back(ch);
		return rts;
	}
	std::vector<std::vector<Residue>> chains_new() {
		std::vector<std::vector<Residue>> rts;
		for(Peptide &p:peps_) for(const auto &ch:p.chs_new()) rts.push_back(ch);
		return rts;
	}
private:
	//Peptide pep_;
	std::vector<Peptide> peps_;
	std::vector<std::vector<NSPproteinrep::PdbRecord>> ligands_;
	std::vector<NSPproteinrep::PdbRecord> h2os_;
};


}


#endif /* ALLATOM_PDBREADER_XY_H_ */
