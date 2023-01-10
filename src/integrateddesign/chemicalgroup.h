/*
 * chemicalgroup.h
 *
 *  Created on: 2018年4月19日
 *      Author: hyliu
 */

#ifndef CHEMICALGROUP_H_
#define CHEMICALGROUP_H_
#include "dataio/inputlines.h"
#include "geometry/xyz.h"
#include "fullsite/fullsite.h"
#include <string>
#include <sstream>
#include <vector>
#include <memory>
#include <map>
#include <set>
using namespace NSPgeometry;

namespace NSPproteinrep{
typedef std::pair<std::string,XYZ> AtomInfo;
typedef std::pair<std::string,std::vector<AtomInfo>> UnitInfo;//==one unit
typedef std::vector<UnitInfo> OneGrpInfo;//==one chemical group
typedef std::map<std::string,std::vector<OneGrpInfo>> GrpsInfo;//==assemble of chemical group
typedef std::map<std::pair<std::string,std::string>,std::vector<std::pair<std::vector<XYZ>,std::vector<XYZ>>>> AtomsPairCrd;
typedef std::map<std::pair<std::string,std::string>,std::vector<std::vector<double>>> UnitPairMatrix;
typedef std::map<std::pair<std::string,std::string>,std::vector<std::vector<XYZ>>> UnitPairCrd;

struct InteractionUnit{
	static int getunittype(const std::string &unitname,
			NSPdataio::InputLines *il=nullptr,int *lidx=nullptr);
	InteractionUnit(std::string un, int ut):name(un),unittype(ut) {}
	std::string name;
	int unittype; //+1:hb donor;-1:hb acceptor;0: both
};
struct ChemGrp {
	static const ChemGrp & getchemgrp(const std::string &name);
	static void readchemgrps(std::map<std::string,ChemGrp> & chemgrps,
			const std::string & filename="chemgrps.dat");
	struct UnitInGrp{
		std::string unitname;
		std::vector<int> atoms;
	};
	//std::multimap<std::string,std::vector<std::string>> modulesandatoms;
	std::map<std::string,std::vector<std::string>> modulesandatoms;
	std::vector<UnitInGrp> units;
	static std::vector<XYZ> buildcrd(std::string grpname);
	bool contain(std::string un) {
		for(UnitInGrp &uig:units) {
			if(uig.unitname==un) return true;
		}
		return false;
	}
};
typedef std::vector<NSPgeometry::XYZ> ChemGrpCrd;

struct GrpInProtein{
	int chainid{-1};
	int resid{-1};
	bool ismainchain{false};
	GrpInProtein(){;}
	GrpInProtein(int cid,int rid,const std::string &gpn,const std::string &mln,
			bool mc=false):chainid(cid),resid(rid),grpname(gpn),molname(mln),
							ismainchain(mc){grp=&(ChemGrp::getchemgrp(gpn));}
	OneGrpInfo unitsingroup();
	std::string grpname;
	std::string molname{""};
	const ChemGrp *grp;
	ChemGrpCrd crd;
	void getcrd(std::string &line) {
		int ix=grp->modulesandatoms.begin()->second.size();
		crd.clear();
		std::istringstream iss(line);
		for(int i=0;i<ix;i++) {
			double x,y,z;
			iss >>x >>y >>z;
			crd.push_back(NSPgeometry::XYZ(x,y,z));
		}
	}
};

struct UnitPairInGrpPair {
	static bool findunitpair(GrpInProtein &gip1, GrpInProtein &gip2, int &n1, int &n2);
	static bool findunitpair(GrpInProtein &gip1, GrpInProtein &gip2, std::string &str1, std::string &str2);
	static void getgrppaircrd(GrpInProtein &gip1, GrpInProtein &gip2, AtomsPairCrd &apc);
	static void getunitpaircrd(GrpInProtein &gip1, GrpInProtein &gip2, AtomsPairCrd &apc);
	static std::vector<double> upmatrix(std::vector<XYZ> &cs1, std::vector<XYZ> &cs2);
	static void getunitpairmatrix(GrpInProtein gip1, GrpInProtein gip2, UnitPairMatrix &upm);
	static std::set<std::pair<std::string,std::string>> unitpairingrppair(GrpInProtein &gip1, GrpInProtein &gip2);
	static void randomunitpair(GrpInProtein &gip1, GrpInProtein &gip2, std::vector<ChemGrpCrd> &cgcs, int ntimes);
	static void printgrppaircrd(AtomsPairCrd &apc, std::string dir);
	static void printunitpaircrd(AtomsPairCrd &apc, std::string dir);
	static void printunitpairmatrix(UnitPairMatrix &upm, std::string dir);
};

class ChemGrpsInProtein {
public:
	void init(std::vector<std::vector<FullSite>> & chains);
	ChemGrpsInProtein(std::vector<std::vector<FullSite>> & chains) {
		init(chains);
	}
	ChemGrpsInProtein(std::string pdbfile) {
		std::vector<std::vector<FullSite>> chains=readfullsitesfrompdb(pdbfile,true);
		init(chains);
	}
	void getgrppaircrd(AtomsPairCrd &apc);
	void getunitpaircrd(AtomsPairCrd &apc);
	void getunitpairmatrix(UnitPairMatrix &upms);
private:
	std::vector<GrpInProtein> gips;
	int mcunitinaa(FullSite &site,std::string chemgrpname,std::string molname);
	int mcunitlinkaa(FullSite &site,FullSite &siteb,std::vector<int> &pb,
			std::string chemgrpname,std::string molname);
	int grpssidechain(FullSite &site);
	int grpsmainchain(std::vector<FullSite> &chain,int posi);
};















void atom3inpdb(std::string pdbfile, GrpsInfo &gis);
void atom3inpdb(std::vector<std::vector<FullSite>> & chains, GrpsInfo &gis);
void atom3inpdb(std::vector<GrpInProtein> &grps, GrpsInfo &gis);
bool getupmx(UnitInfo &u1, UnitInfo &u2, std::vector<double> &ds);
void printupmatrix(UnitPairMatrix &upms, std::ostream &os);
int mcunitinaa(FullSite &site,std::vector<GrpInProtein> &grps,
		std::string chemgrpname,std::string molname);
int mcunitlinkaa(FullSite &site,FullSite &siteb,std::vector<int> &pb,
		std::vector<GrpInProtein> &grps,std::string chemgrpname,std::string molname);
struct UnitPair {
	struct UnitPairGeom{
		double rda;
		double thd,tha;
		double torsiond,torsiona,torsionda;
		UnitPairGeom(const ChemGrp::UnitInGrp *donorunit,
				const ChemGrpCrd &crd1, const ChemGrp::UnitInGrp *acceptorunit,
				const ChemGrpCrd &crd2,double r);
		double match(const UnitPairGeom & g2) const;
	};
	const ChemGrp::UnitInGrp *donorunit{nullptr};
	const ChemGrp::UnitInGrp *acceptorunit{nullptr};
	std::shared_ptr<UnitPairGeom> geometry;
	UnitPair(){;}
	UnitPair(const ChemGrp & cg1, const ChemGrpCrd &crd1,
			const ChemGrp &cg2,const ChemGrpCrd &crd2);
	bool valid() const { return (donorunit and acceptorunit);}
	std::string name() const {return donorunit->unitname
			+"_"+acceptorunit->unitname;}
	double match(const UnitPair &up2) const{
		if(name() != up2.name()) return 0.0;
		return geometry->match(*(up2.geometry));
	}
};
UnitPair makeunitpair(const GrpInProtein &gip1, const GrpInProtein &gip2);



std::vector<GrpInProtein> grpsinpdb(std::vector<std::vector<FullSite>> & chains);
std::vector<GrpInProtein> grpsinpdb(std::string pdbfile);
int mcunitinaa(FullSite &site,std::vector<GrpInProtein> &grps,
		std::string chemgrpname,std::string molname);
int mcunitlinkaa(FullSite &site,FullSite &siteb,std::vector<int> &pb,
		std::vector<GrpInProtein> &grps,std::string chemgrpname,std::string molname);
int grpssidechain(FullSite &site,std::vector<GrpInProtein> &grps);
int grpsmainchain(std::vector<FullSite> &chain,int posi,std::vector<GrpInProtein> &grps);
void getunitpaircrd(std::vector<GrpInProtein> &gips, UnitPairCrd &upms);
void getunitpairmatrix(std::vector<GrpInProtein> &gips, UnitPairMatrix &upms);
}



#endif /* CHEMICALGROUP_H_ */
