/*
 * filter.h
 *
 *  Created on: Jan 14, 2019
 *      Author: xuyang
 */

#ifndef ALLATOM_FILTER_H_
#define ALLATOM_FILTER_H_

#include "pdbreader_xy.h"
#include "basic_info.h"
#include "sd/forcefield.h"

namespace NSPallatom {
class Filter000000000 {
public:
	/*
	 * return proportion of hydrogen bond formation in the given betapairs
	 */
	static std::vector<double> betasheets(const std::vector<std::vector<Residue>> &chains,
			const std::vector<std::pair<int,int>,std::pair<int,int>>& betapairs);
	/*
	 * find residues without hydrogen bond and embeded inner of protein
	 * return the sequence number of residue
	 */
	static std::vector<std::pair<int,int>> mainchainHB(const std::vector<std::vector<Residue>> &chs);
	/*
	 * find residues without hydrogen bond and embeded inner of protein
	 * return the sequence number of residue
	 */
	static std::vector<std::pair<int,int>> allatomHB(const std::vector<std::vector<Residue>> &chs);
private:
};

struct StrideRecord {
	StrideRecord(std::string line);
	std::string resname;
	char chainid;
	int resnumber;
	int resorder;
	char sscode;
	std::string ssname;
	double phi;
	double psi;
	double sasa;
};
std::vector<StrideRecord> PDB2Stride(std::string pdbfile);
//loop, helix, strand,,,,,,,   chainid, st, en
std::map<std::string,std::vector<int>> getssfromstride(const std::vector<StrideRecord> &srs);

class Filter {
public:
	typedef std::pair<std::vector<int>,std::string> atominres; //chainid, resid, atomname
	enum ENE {
		TOTAL, BOND, ANG, IMPDIH, STERIC, SCCONF, PHIPSI, LOCAL, SCPACKING, LOCALBBHB, SITEPAIR
	};
	static std::vector<std::string> eneterm;

	std::vector<std::vector<Residue>> & chs() {return chs_;}

	struct Ene_Res {
		int chainid{-1};
		int resid{-1};
		double ephipsi{BigValue}, elocal{BigValue}, escconf{BigValue}, etot{0.0};
		std::map<std::pair<int,int>,double> esitepair, escpacking, elocalbbhb;
	};
	struct Ene_Atom {
		int chainid{-1};
		int resid{-1};
		std::string atomname;
		std::vector<std::pair<std::vector<atominres>,double>> ebond, eang, eimpdih, esteric;
	};
	void getenes();
	std::vector<double> &totalene() {return etot_;}
	std::vector<std::vector<Ene_Res>> & e_res() {return enes_res;}
	std::vector<std::vector<std::map<std::string,Ene_Atom>>> & e_atom() {return enes_atom;}
	typedef std::pair<std::string,std::string> AtomName;
	struct EneInKind {
		std::map<std::pair<AtomName,AtomName>,std::vector<double>> ebond,esteric;
		std::map<std::pair<AtomName,std::pair<AtomName,AtomName>>,std::vector<double>> eang;
		std::map<std::pair<std::pair<AtomName,AtomName>,std::pair<AtomName,AtomName>>,std::vector<double>> eimpdih;
		std::map<std::string,std::vector<double>> escconf,ephipsi,elocal;
		std::map<std::pair<std::string,std::string>,std::vector<double>> escpacking,elocalbbhb,esitepair;
	};
	EneInKind geteneinkind();

	struct EneUnit{
		//total number
		int nbond{0},nang{0},nimpdih{0},nsteric{0},nscconf{0},nphipsi{0},nlocal{0},nscpacking{0},nsitepair{0},nlocalbbhb{0};
		//energy unit
		std::vector<std::vector<atominres>> ubond,uang,uimpdih,usteric;
		std::vector<std::vector<int>> uscconf,uphipsi,ulocal; //size==2
		std::vector<std::vector<int>> uscpacking,usitepair,ulocalbbhb; //size==4
		//energy island
		std::vector<std::vector<int>> ibond,iang,iimpdih,isteric,iscconf,iphipsi,ilocal,iscpacking,isitepair,ilocalbbhb;
		//energy, interval
		std::pair<double,int> pscconf;
		std::pair<double,int> pphipsi;
		std::pair<double,int> plocal;
		std::pair<double,int> plocalbbhb; //the second is distance in sequence
		std::pair<double,double> pbond;
		std::pair<double,double> pang;
		std::pair<double,double> pimpdih;
		std::pair<double,double> psteric;
		std::pair<double,double> pscpacking;
		std::pair<double,double> psitepair;
		//maxisland size cutoff, if more than this, the deleted
		int szbond, szang, szimpdih, szsteric, szscconf, szscpacking, szphipsi, szlocal, szlocalbbhb, szsitepair;
		//view as 1 island if distance smaller than this par
		//int disscconf, disphipsi, dislocal, dislocalbbhb;
		//double dissteric, disscpacking, dissitepair;
		//allowable max energy [1] & island energy [0]   ->    cutoff
		std::map<std::string,std::vector<double>> cutscconf, cutphipsi, cutlocal;
		std::map<std::vector<std::string>,std::vector<double>> cutlocalbbhb, cutscpacking, cutdissitepair, cutsteric;
	};
	void geteneunitpar();
	void island_atom(std::string label);
	void island_respair(std::string label);
	void island_single(std::string label);
	void island_bond(double e, double d);
	void island_ang(double e, double d);
	void island_impdih(double e, double d);
	void island_steric(double e, double d);
	void island_scconf(double e, int d);
	void island_phipsi(double e, int d);
	void island_local(double e, int d);
	void island_localbbhb(double e, int d);
	void island_scpacking(double e, double d);
	void island_sitepair(double e, double d);
	void eneisland(std::map<std::string,std::pair<double,double>> &lb1, std::map<std::string,std::pair<double,int>> &lb2);
	void eneisland1();
	EneUnit & geteneunit() {return highregion_;}
	std::vector<double> gethighene();
	std::string heainlocalregion();
	void printisland(std::ostream &os);

	struct HB {
		int chainid{-1};
		int resid{-1};
		std::string resname;
		std::map<std::string,std::set<atominres>> hbs;
	};
	void gethbps();
	std::vector<std::vector<HB>> & hbps() {return hbps_;}
	void mainchainHBfilter(std::vector<std::vector<int>>&helixs, std::vector<std::vector<int>> &hbh,
			std::vector<std::vector<int>>&strands, std::vector<std::vector<int>> &hbs); //start, end+1
	std::vector<std::vector<double>> mainchainHBfilter(std::vector<std::vector<int>>&helixs, std::vector<std::vector<int>>&strands);
	void findinnerpolar(int pmax, double rad, std::vector<atominres>&pas);

	static std::vector<std::vector<int>> looplthfilter(std::vector<std::vector<int>>&acts, std::string pdbfile); // acts, chainid, st, en
private:
	std::vector<std::vector<Residue>> chs_;
	//std::map<int,std::map<int,int>> fixed_; // chain id, start, end
	std::vector<std::vector<Ene_Res>> enes_res;
	std::vector<std::vector<std::map<std::string,Ene_Atom>>> enes_atom;
	std::vector<double> etot_;
	void assign_ene(std::string label,
			const std::vector<std::pair<std::vector<std::pair<std::vector<int>,std::string>>,double>> &es);
	void assign_ene(std::string label,
			const std::vector<std::pair<std::vector<std::vector<int>>,double>> &es);
	void assign_ene(std::string label,
			const std::vector<std::pair<std::vector<int>,double>> &es);
	void totalene_site();
	void totalene_eterm();

	std::vector<double> enecutoff_;
	EneUnit highregion_;
	std::vector<std::vector<int>> island_;
	bool equalatoms(const std::vector<atominres> &air1, const std::vector<atominres> &air2);
	void island_space(const std::vector<std::vector<NSPgeometry::XYZ>> &crds,
			double discut, std::vector<std::vector<int>> &ild);
	void island_sequence(const std::vector<std::vector<int>>&hrs, int discut, std::vector<std::vector<int>>&ild);
	std::vector<std::vector<int>> eneisland();
	void printislandatom(std::vector<atominres> &ar, std::string lb, std::ostream &os);
	void printislandres(std::vector<int> &rs, std::string lb, std::ostream &os);

	//EneInKind eik_;


	std::vector<std::vector<HB>> hbps_;
	std::vector<int> hbinhelix(int cid, int st, int en);
	std::vector<int> hbinstrand(int c0, int s0, int e0, int c1, int s1, int e1);
	int exposed(std::vector<NSPgeometry::XYZ>&cs, double rad);
};
}


#endif /* ALLATOM_FILTER_H_ */
