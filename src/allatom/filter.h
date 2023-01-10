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

typedef std::vector<int> SingleResidue; //size==2, chainid, residueid
typedef std::pair<SingleResidue,double> EResidue; //energy
typedef std::vector<SingleResidue> IResidue; //island

typedef std::vector<SingleResidue> MultiResidue; //size==2
typedef std::pair<MultiResidue,double> EMultiResidue;
typedef std::vector<MultiResidue> IMultiResidue;

typedef std::pair<SingleResidue,std::string> SingleAtom;
typedef std::vector<SingleAtom> MultiAtom;
typedef std::pair<MultiAtom,double> EAtom; // 0->chainid, 1->residueid
typedef std::vector<MultiAtom> IAtom;

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


class EnergyFilter {
public:
	enum ENE {
		TOTAL, BOND, ANG, IMPDIH, STERIC, SCCONF, PHIPSI, LOCAL, S1, S2, SCPACKING, LOCALBBHB, SITEPAIR, ETERM
	};
	static std::vector<std::string> eneterm;
	static std::map<std::string,int> etseq;
	static std::set<std::string> ipdint;
	EnergyFilter(const std::vector<std::vector<Residue>> &chs) {chains_ = chs;}
	void setislandpar();
	void mcene();
	void scene(std::string abacusfile, bool readss);
	void sitetotalene();
	void enelist();
	void island();
	std::vector<std::string> enetermfilter(); // bigger than max island
	std::map<std::string,std::pair<int,double>> highene();
	std::map<std::string,double> totalene();

	std::string mainchainfilter(const std::vector<std::vector<Residue>> &chs);
	std::string allatomfilter(const std::vector<std::vector<Residue>> &chs, std::string abacusfile);
private:
	std::vector<std::vector<Residue>> chains_;
	//secondary structure
	std::vector<std::string> ss_;
	//energy term (contain repeat)
	std::vector<std::vector<std::vector<EAtom>>> ebond_, eang_, eimpdih_, esteric_;
	std::vector<std::vector<std::vector<EMultiResidue>>> esitepair_, escpacking_, elocalbbhb_, es2_;
	std::vector<std::vector<double>> ephipsi_, elocal_, escconf_, es1_, etot_;
	//energy list (no repeat)
	std::map<MultiAtom,double> lbond_, lang_, limpdih_, lsteric_;
	std::map<SingleResidue,double> lphipsi_, llocal_, lscconf_, ls1_, ltot_;
	std::map<MultiResidue,double> lsitepair_, lscpacking_, llocalbbhb_, ls2_;
	//island parameter, energy, distance (int double), maxislandsize
	std::vector<double> parene_;
	std::vector<int> pardisint_;
	std::vector<double> pardisdou_;
	std::vector<int> maxisland_;
	//island
	std::vector<IAtom> ibond_, iang_, iimpdih_, isteric_;
	std::vector<IMultiResidue> isitepair_, iscpacking_, ilocalbbhb_, is2_;
	std::vector<IResidue> iphipsi_, ilocal_, iscconf_, is1_, itot_;


	void assign_ene(std::vector<std::vector<std::vector<EAtom>>> &et,
			const std::vector<std::pair<std::vector<SingleAtom>,double>> &es);
	void assign_ene(std::vector<std::vector<std::vector<EMultiResidue>>> &et,
			const std::vector<std::pair<std::vector<std::vector<int>>,double>> &es);
	void assign_ene(std::vector<std::vector<double>> &et,
			const std::vector<std::pair<std::vector<int>,double>> &es);

	void enelist(std::map<SingleResidue,double>&el, const std::vector<std::vector<double>>&es);
	void enelist(std::map<MultiResidue,double>&el,
			const std::vector<std::vector<std::vector<EMultiResidue>>>&es);
	void enelist(std::map<MultiAtom,double>&el,const std::vector<std::vector<std::vector<EAtom>>>&es);

	NSPgeometry::XYZ getcrd(const SingleAtom &sa) {
		return chains_[sa.first[0]][sa.first[1]].rds().at(sa.second).crd;
	}
	std::vector<NSPgeometry::XYZ> getcrds(const SingleResidue &sr, bool onlymc) {
		std::set<std::string> mc=AA_Atom::mainchain();
		std::vector<NSPgeometry::XYZ> cs;
		const auto & rs = chains_[sr[0]][sr[1]].rds();
		for(auto &r:rs) {
			if(onlymc) {
				if(mc.find(r.first)==mc.end()) continue;
			}
			cs.push_back(r.second.crd);
		}
		return cs;
	}
	void island_atom(std::vector<IAtom>&il,const std::map<MultiAtom,double>&el,double ec,double ed);
	void island_singleresidue(std::vector<IResidue>&il,const std::map<SingleResidue,double>&el,double ec,int ed);
	void island_multiresidue_int(std::vector<IMultiResidue>&il,const std::map<MultiResidue,double>&el,double ec,int ed);
	void island_multiresidue_double(std::vector<IMultiResidue>&il,const std::map<MultiResidue,double>&el,double ec,double ed,bool onlymc);

};



typedef std::vector<std::vector<Residue>> DECOY;

struct EnergyTermRecord {
	DECOY decoy;
	std::map<MultiAtom,double> ebond, eang, eimpdih, esteric;
	std::map<SingleResidue,double> ephipsi, elocal, escconf, etot, es1, erot;
	std::map<MultiResidue,double> esitepair, escpacking, elocalbbhb, es2, epck;
	std::map<SingleResidue,char> sss;
	std::map<SingleResidue,std::string> seq;
	void readss(std::string ssfile);
	//std::map<char,std::map<std::string,std::vector<double>>> es1, es2;
	//std::vector<std::pair<std::pair<char,std::string>,double>> es1,es2;
	//void reads1s2(std::string efile);
	EnergyTermRecord() {
		;
	}
	void readenergy(std::string enefile, bool mconly=false);
	std::vector<std::vector<int>> v1_v2();
	void read_s1_rot(std::string enefile, std::map<SingleResidue,double> &es);
	void read_s2_pck(std::string enefile, std::map<MultiResidue,double> &es);

	void eneinit(NSPsd::EnergyTerm &et);
	EnergyTermRecord(const DECOY &dy, NSPsd::EnergyTerm &et):decoy(dy) {
		eneinit(et);
	}
	EnergyTermRecord(NSPsd::EnergyTerm &et) {
		eneinit(et);
	}
	void rank(MultiAtom &ma) {
		std::sort(ma.begin(),ma.end(),[](SingleAtom sa1, SingleAtom sa2)->bool{
			return (sa1.first[0]<sa2.first[0]?true:(sa1.first[0]>sa2.first[0]?false:(
					sa1.first[1]<sa2.first[1]?true:(sa1.first[1]>sa2.first[1]?false:(
							sa1.second<sa2.second)))));});
	}
	void rank(MultiResidue &mr) {
		std::sort(mr.begin(),mr.end(),[](std::vector<int>v1,
				std::vector<int>v2)->bool{return (v1[0]<v2[0]?true:
						(v1[0]>v2[0]?false:(v1[1]<v2[1])));});
	}
	void printene(std::string outene);
	void print(std::string outpdb, std::string outene) {
		residue2pdb(outpdb,decoy);
		printene(outene);
	}

	EnergyTermRecord(std::string pdbfile, std::string enefile) {
		PdbReader_xy pr;
		pr.init(pdbfile);
		decoy = pr.chains_new();
		for(int i=0;i<decoy.size();i++) {
			for(int j=0;j<decoy[i].size();j++) {
				seq.insert({{i,j},decoy[i][j].resname()});
			}
		}
		readenergy(enefile);
	}
	EnergyTermRecord(std::string enefile) {
		readenergy(enefile);
	}
	std::pair<MultiAtom,double> readmultiatom(std::stringstream &ss, int na);
	std::pair<MultiResidue,double> readmultiresidue(std::stringstream &ss);
	std::pair<SingleResidue,double> readsingleresidue(std::stringstream &ss);

	static bool allatomismc(const MultiAtom &ma);
	double totalene();
	EnergyTermRecord(const DECOY &dy);
	static std::map<SingleResidue,double> singleresidueenergy(const std::map<MultiResidue,double>&mes);
};

class TrajectoryEnergy {
public:
	TrajectoryEnergy(DECOY &dy,std::vector<std::vector<int>>&f):decoy_(dy), fixed_(f) {
		;
	}
	TrajectoryEnergy(const std::vector<std::string> &pdbfs, const std::vector<std::string> &enefs,
			const std::vector<std::vector<int>>&f):fixed_(f) {
		assert(pdbfs.size()==enefs.size());
		for(int i=0;i<pdbfs.size();i++) {
			eneterm_.push_back(EnergyTermRecord(pdbfs[i],enefs[i]));
		}
	}
	TrajectoryEnergy(const std::vector<std::string> &enefs) {
		for(int i=0;i<enefs.size();i++) {
			std::cout <<i <<std::endl;
			eneterm_.push_back(EnergyTermRecord(enefs[i]));
		}
	}
	void dosd(int nstep, int seed, std::map<std::string,double> &sdctrl, std::map<std::string,double> &ffctrl);
	enum ENE {
		TOTAL, BOND, ANG, IMPDIH, STERIC, SCCONF, PHIPSI, LOCAL, SCPACKING, LOCALBBHB, SITEPAIR, ETERM
	};
	void readcut(std::string parfile);
	struct HighEnergyTime {
		std::map<MultiAtom,std::vector<int>> tbond, tang, timpdih, tsteric;
		std::map<SingleResidue,std::vector<int>> tphipsi, tlocal, tscconf, ttot;
		std::map<MultiResidue,std::vector<int>> tsitepair, tscpacking, tlocalbbhb;
	};
	HighEnergyTime & highenergytime() {return het_;}
	struct EnergyMerge {
		std::map<MultiAtom,std::vector<double>> ebond, eang, eimpdih, esteric;
		std::map<SingleResidue,std::vector<double>> ephipsi, elocal, escconf, etot;
		std::map<MultiResidue,std::vector<double>> esitepair, escpacking, elocalbbhb;
		void standard(const std::vector<double> &cut11, const std::vector<double> &cut22);
	};
	EnergyMerge & energymerge() {return em_;}
	struct EnergyAverage {
		std::map<MultiAtom,double> ebond, eang, eimpdih, esteric;
		std::map<SingleResidue,double> ephipsi, elocal, escconf, etot;
		std::map<MultiResidue,double> esitepair, escpacking, elocalbbhb;
	};
	EnergyAverage & energyaverage() {return ea_;}
	struct HighEnergy {
		std::map<MultiAtom,double> ebond, eang, eimpdih, esteric;
		std::map<SingleResidue,double> ephipsi, elocal, escconf, etot;
		std::map<MultiResidue,double> esitepair, escpacking, elocalbbhb;
	};
	HighEnergy & highenergy() {return he_;}
	void energyanalysis(int mintime=75);
	std::map<std::string,double> ene2rank();
	std::vector<EnergyTermRecord> eneterm() {return eneterm_;}
	static double index01(double c1, double c2, double v); //same to below
	std::vector<double> c0() {return cut0;}
	std::vector<double> c1() {return cut1;}
	std::vector<double> c2() {return cut2;}
private:
	std::vector<EnergyTermRecord> eneterm_;
	DECOY 	decoy_;
	std::vector<std::vector<int>> fixed_; //chainid, st, en
	std::vector<double> cut0; // to delete
	std::vector<double> cut1; // high bound
	std::vector<double> cut2; // low bound
	HighEnergyTime het_;
	EnergyMerge em_;
	EnergyAverage ea_;
	HighEnergy he_;
	double index(double c1, double c2, double v); // 0.0 =< ? <= 1.0
	void changeweight(std::vector<std::string>&ctrl, std::string label, double wgt);
};











class Filter {
public:
	Filter(const std::vector<std::vector<Residue>> &chs) {chains_ = chs;}
	static std::vector<std::vector<int>> looplthfilter(
			std::vector<std::vector<int>>&acts, std::string pdbfile); //chainid, st, en
	struct HB {
		int chainid{-1};
		int resid{-1};
		std::string resname;
		std::map<std::string,std::set<SingleAtom>> hbs;
	};
	void mainchainHBfilter(const std::vector<std::vector<int>>&helixs, std::vector<std::vector<int>> &hbh,
			const std::vector<std::vector<int>>&strands, std::vector<std::vector<int>> &hbs); //start, end+1
	std::vector<std::vector<double>> mainchainHBfilter(
			const std::vector<std::vector<int>>&helixs, const std::vector<std::vector<int>>&strands);
	std::vector<std::map<std::string,std::vector<std::pair<SingleAtom,int>>>> mcpolar(); // 0->NOT in HB, 1->in HB
	//std::vector<std::map<std::string,std::vector<std::pair<SingleAtom,int>>>> mcpolar(); // 0->NOT in HB, 1->in HB
	//std::vector<std::map<std::vector<std::string>,std::vector<int>>> expose(bool considerHB=true); // 0->NOT in HB, 1->in HB
	void allatomexpose();
	void allatomexpose_psudosc(); // just for sidechain is psudo
	std::map<SingleAtom,int> exp() {
		//if(expose_.empty()) allatomexpose();
		return expose_;
	}
	void expose_remove_self();
	static DECOY psudosidechain(const DECOY &chs);
	std::set<std::vector<SingleAtom>> gethbpair(bool diffpro=false);
private:
	std::vector<std::vector<Residue>> chains_;
	//std::vector<std::vector<int>> acts_; //chainid, st, en
	double rad_add_{0.2}; //radius_a1 + radius_a2 + rad_add_ > distance_a1_a2  ==> HB
	std::vector<std::vector<HB>> hbps_;
	void gethbps(bool diffpro=true);
	std::vector<int> hbinhelix(int cid, int st, int en);
	std::vector<int> hbinstrand(int c0, int s0, int e0, int c1, int s1, int e1);
	std::map<SingleAtom,int> expose_;
	void mcpolarexpose(); // add psudo residue
};

bool mainchainHBfilter(const std::vector<std::vector<Residue>>&chs,
		const std::vector<std::vector<int>>&helixs, const std::vector<std::vector<int>>&strands);







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

/*
class Filter0 {
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
};*/
}


#endif /* ALLATOM_FILTER_H_ */
