/*
 * pdbreader_xy.h
 *
 *  Created on: Aug 13, 2018
 *      Author: xuyang
 */

#ifndef ALLATOM_PDBREADER_XY_OLD_H_
#define ALLATOM_PDBREADER_XY_OLD_H_

//#define Square_Length 2.0
#define Water_Radius 1.4
#define Atom_Radius 1.5

#include "proteinrep/pdbrecord.h"
#include "fullsite/fullsite.h"
#include "sd/stochasticdynamics.h"
using namespace NSPgeometry;

namespace NSPallatom {

/*
 * Hydrogen Bond Criterion:
 * 1.	distance between donor and accepter is less than sum of their radius;
 * 2.	SER OG or THR OG1 could not couple with O or N in the same residue.
 */

struct AtomPair {
	static std::vector<std::vector<std::string>> crdmainchain;
	static std::map<std::string,std::vector<std::vector<std::string>>> crdsidechain;
	static std::set<std::string> residues;
	static std::set<std::string> getresidues();
	static std::set<std::string> mainchains;
	static std::set<std::string> getmainchains();
	static std::map<std::string,std::set<std::string>> sidechains;
	static std::map<std::string,std::set<std::string>> getsidechains();
	typedef std::pair<std::string,std::string> AtomName;
	static std::set<AtomName> atomnames;
	static std::set<AtomName> getatomnames();
	static std::map<AtomName,std::pair<AtomName,AtomName>> neighbors;
	static std::map<AtomName,std::pair<AtomName,AtomName>> getneighbors();
	static std::map<AtomName,double> radius;
	static std::map<AtomName,double> getradius();
	static std::set<AtomName> donors;
	static std::set<AtomName> getdonors();
	static std::set<AtomName> accepters;
	static std::set<AtomName> getaccepters();

	std::string resname;
	std::string cen;
	std::string fi;
	std::string se;
	std::string envresname;
	std::string envatomname;
	std::vector<XYZ> crds; //first is cen, second is fi, third is se, others are env
};

class Residue {
public:
	void init(const std::vector<NSPproteinrep::PdbRecord> &prs);
	int & chainseq() {return chainseq_;}
	int & resseq() {return resseq_;}
	std::map<std::string,NSPgeometry::XYZ> & crds() {return crds_;}
	std::string & resname() {return resname_;}
	const int & chainseq() const {return chainseq_;}
	const int & resseq() const {return resseq_;}
	const std::map<std::string,NSPgeometry::XYZ> & crds() const {return crds_;}
	const std::string & resname() const {return resname_;}
	NSPproteinrep::FullSite Tofullsite();
	std::map<std::pair<std::string,std::string>,std::vector<XYZ>> getneighbor();
	std::vector<XYZ> covert20(std::string resn="");
	char & chid() {return chainid_;}
	int & resid() {return resid_;}
	char & ssid() {return ss_;}
	void printrecords(int &atominitnum, std::ostream &os) const;
private:
	char chainid_{'A'};
	int resid_{-1000};
	std::string resname_{"XXX"};
	char insertionid_{' '};
	char ss_{'X'};

	int chainseq_{-1};
	int resseq_{-1};
	//bool complete_{true};
	std::map<std::string,NSPgeometry::XYZ> crds_;
};

class PdbReader_xy {
public:
	typedef std::vector<Residue> Chain;
	void init(std::string pdbfile, std::string out20="");
	//first is cen, second is fi, third is se
	std::map<std::pair<std::string,std::string>,std::vector<XYZ>> neighbors();
	std::map<std::pair<std::string,std::string>,std::vector<XYZ>> peps();
	int print(std::ostream &os);
	XYZ center_pep();
	double maxrad_pep() {
		XYZ center=center_pep();
		return maxrad_pep(center);
	}
	double maxrad_pep(XYZ &center);
	static NSPproteinrep::PdbRecord getarecord(std::string resname, std::string atomname, const XYZ &c);
	std::string & pdbid() {return pdbid_;}
	std::vector<Chain> &getchain() {return protein_;}
	void addss(std::string pdbfile);
	std::vector<std::vector<NSPproteinrep::PdbRecord>> &getligands() {return ligands_;}
	struct OneAtom {
		int chainid;
		int resid;
		std::string atomname;
		OneAtom(int c,int r,std::string a):chainid(c),resid(r),atomname(a) {;}
	};
	std::vector<std::pair<OneAtom,OneAtom>> hbpairs() const;
	std::vector<std::string> seqs();
/*
	void getpar(std::vector<NSPgeometry::XYZ>&fixeds,std::vector<double>&rads,
			NSPgeometry::XYZ &center,double &maxdis,double wr,double ar) const;

	//fixeds contain protein_ and ligands_, rads is coresponding to fixeds
	//center and maxdis are calculated just from protein_
	void getpar_grid(std::vector<NSPgeometry::XYZ>&fixeds,std::vector<double>&rads,
			NSPgeometry::XYZ &center,double &maxdis) const ;

	std::vector<Chain> protein() {return protein_;}
	Chain pep() {
		Chain ch;
		for(Chain &c:protein_) {
			for(Residue &r:c) ch.push_back(r);
		}
		return ch;
	}
	std::map<std::pair<std::string,std::string>,std::vector<XYZ>> crds();*/
private:
	std::vector<Chain> protein_; //if assemble has "N", "CA", "C", "O"
	std::set<std::pair<int,int>> beyond20s_;
	std::vector<std::vector<NSPproteinrep::PdbRecord>> ligands_;
	std::vector<NSPproteinrep::PdbRecord> waters_;
	std::string pdbid_;
};

/*struct Square {
	double x;
	double y;
	double z;
	Square(double a, double b, double c):x(a),y(b),z(c) {;}
	//Square(const XYZ &c) {
	//	x=c.x_;
	//	y=c.y_;
	//	z=c.z_;
	//}
	//double length{0.0};
	//std::vector<std::pair<std::pair<std::string,std::string>,std::vector<XYZ>>> nghs;
	//std::vector<XYZ> lgds;
	//std::vector<XYZ> wats;
};*/
typedef XYZ Square;
/*
bool operator < (const Square &r1, const Square &r2) {
	double halflth=Square_Length/2.0;
	if(r1.x_<r2.x_-halflth) return true;
	else if(r1.x_>r2.x_+halflth) return false;
	else if(r1.y_<r2.y_-halflth) return true;
	else if(r1.y_>r2.y_+halflth) return false;
	else if(r1.z_<r2.z_-halflth) return true;
	else return false;
}*/
/*
bool operator < (const Range &r1, const Range &r2) {
	if(r1.x1<r2.x1) return true;
	else if(r1.x1>r2.x1) return false;
	else if(r1.y1<r2.y1) return true;
	else if(r1.y1>r2.y1) return false;
	else if(r1.z1<r2.z1) return true;
	else return false;
}


bool range_comp(Range r1, Range r2) {
	double x=r1.x1-r2.x1;
	if(x<0.0) {
		if(x>-PRECISION)
	}
	if(r1.x1<r2.x1) return true;
	else if(r1.x1>r2.x1) return false;
	else if(r1.y1<r2.y1) return true;
	else if(r1.y1>r2.y1) return false;
	else if(r1.z1<r2.z1) return true;
	else return false;
}

typedef std::map<Range,std::vector<XYZ>,range_comp> space_part;*/

struct Grid {
	static std::vector<NSPgeometry::XYZ> getwater(
			const NSPgeometry::XYZ &cen, double rad, double wtrad=Water_Radius);
	static std::vector<Square> getrange(const NSPgeometry::XYZ &cen, double rad);
};

class MixModel {
public:
	void init(PdbReader_xy &p);
	MixModel(PdbReader_xy &pr, double r):rad_add(r) {
		init(pr);
	}
	MixModel(std::string fn, double r):rad_add(r) {
		PdbReader_xy pr;
		pr.init(fn);
		init(pr);
	}
	MixModel() {;}
	void rot_random();
	void microenvironment(std::string outpath, bool calatom, bool calwater);
	typedef std::pair<AtomPair::AtomName,AtomPair::AtomName> APS;
	void microenvironment(std::map<APS,std::vector<std::vector<XYZ>>> &atompairs, bool calatom, bool calwater);
	void statistic(int ntime, std::string outpath, bool calatom, bool calwater);
	std::map<APS,std::vector<std::vector<XYZ>>> statistic(int ntime, bool calatom, bool calwater);
	std::string & pdbid() {return pdbid_;}
	static void print(std::string fn, std::map<APS,std::vector<std::vector<XYZ>>> &atompairs, bool fixcrd=true);
	XYZ & center() { return center_; }
	std::vector<XYZ> & waters() { return waters_; }
	double & radadd() {return rad_add;}
	double & radwater() {return rad_water;}
	double & radsquare() {return rad_square;}
	struct Content {
		//nghs is same with nghs_
		std::map<std::pair<std::string,std::string>,std::vector<XYZ>> nghs;
		std::vector<XYZ> wats;
		std::vector<int> w_label;//1-clash, 0-no clash
	};
	void partion(bool calwater, double buff_rad);
	void nearregion(std::vector<XYZ> &atoms, const XYZ &cen, double max_rad, int n_add, double atomrad, double buffrad);
	void watermenv(std::map<APS,std::vector<std::vector<XYZ>>> &atompairs, double max_rad, double buffrad);
	bool getrt(const XYZ &crd, const XYZ &axis, Rotation &rt);
	bool getrt(const XYZ &crd, const XYZ &axis, const XYZ &plane, Rotation &rt);
	void rot(std::map<MixModel::APS,std::vector<std::vector<NSPgeometry::XYZ>>> &atompairs);
	void waterisclash(const XYZ & cen, int n_add, double dis);
private:
	double rad_add{0.0};
	NSPgeometry::XYZ center_;
	double rad_water{0.0};
	std::vector<NSPgeometry::XYZ> waters_;
	double rad_square{0.0};
	std::map<Square,Content> ranges_;
	void add_square(const XYZ &crd, XYZ &ori);
	//3 XYZ is a unit
	std::map<std::pair<std::string,std::string>,std::vector<XYZ>> nghs_;
	//1 XYZ is a unit
	std::map<std::pair<std::string,std::string>,std::vector<XYZ>> crds_;
	std::string pdbid_{""};
};

class SDWater {
public:
	SDWater(NSPdataio::ParameterSet &pset);
	void add_water();
	double stericene(const XYZ &c1,const XYZ &c2,std::vector<XYZ>*fs);
	std::vector<double> forces(const std::vector<double> &cs, const std::set<std::pair<int,int>>&nbl, double f_cen, double discut);
	std::set<std::pair<int,int>> getneighbor(double rcut2);
	void runstep(int nstep, bool isopt, int itvl, std::string path);
	void move2center();
	void opt();
	void printtraj(std::string fn);
	void printdis(std::string fn);
	std::vector<double> diswater();
	double maxrad() {return maxrad_;}
	void printene(std::string fn1, std::string fn2);
private:
	NSPdataio::ParameterSet pset_;
	PdbReader_xy pr_;
	std::vector<NSPgeometry::XYZ> fixeds_;
	std::vector<double> rads_;
	NSPgeometry::XYZ center_;
	double maxrad_;
	std::vector<NSPgeometry::XYZ> waters_;
	std::vector<bool> forceoff_;
	std::vector<double> masses_;
	std::vector<double> potenes_;
	std::shared_ptr<NSPsd::StochasticDynamics> sd_ { nullptr };
	std::shared_ptr<NSPsd::StochasticDynamics::State> state_ { nullptr };
	std::shared_ptr<NSPsd::StochasticDynamics::State> buffstate_ { nullptr };
};

typedef NSPdataio::TypedMultiInstanceControls<SDWater> SDWaterControls;
void definesdwatercontrol(std::string name,const std::vector<std::string> &controllines);
void readcontrols_sdwater(const std::string &filename,std::string name);

struct BLOSUM62 {
	static std::map<std::pair<char,char>,int> matrix;
	static std::map<std::pair<char,char>,int> getmatrix();
};
}


#endif /* ALLATOM_PDBREADER_XY_OLD_H_ */
