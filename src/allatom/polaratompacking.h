/*
 * polaratompacking.h
 *
 *  Created on: Apr 25, 2019
 *      Author: xuyang
 */

#ifndef ALLATOM_POLARATOMPACKING_H_
#define ALLATOM_POLARATOMPACKING_H_
#include "allatom/pdbreader_xy.h"
#include "allatom/basic_info.h"
#include "geometry/localframe.h"
#include <memory>
using namespace NSPgeometry;

namespace NSPallatom {

struct HBUnit {
	XYZ p0;
	XYZ c1;
	XYZ c2;
	int chainseq;
	int resseq;
	std::string resname;
	std::string np0, nc1, nc2;
	HBUnit(XYZ &crd0, XYZ &crd1, XYZ &crd2, std::string rn, std::string n0, std::string n1, std::string n2,
			int cs, int rs):p0(crd0),c1(crd1),c2(crd2),resname(rn),np0(n0),nc1(n1),nc2(n2),chainseq(cs),resseq(rs) {
		;
	}
};

struct HBUnitPair {
	std::shared_ptr<HBUnit> u1 { nullptr };
	std::shared_ptr<HBUnit> u2 { nullptr };
	HBUnitPair(std::shared_ptr<HBUnit> h1, std::shared_ptr<HBUnit> h2):u1(h1),u2(h2) {
		;
	}
	double dp0() {return sqrt((u1->p0-u2->p0).squarednorm());}
	double dc1() {return sqrt((u1->c1-u2->c1).squarednorm());}
	double dc2() {return sqrt((u1->c2-u2->c2).squarednorm());}
	double dp0c1() {return sqrt((u1->p0-u2->c1).squarednorm());}
	double dp0c2() {return sqrt((u1->p0-u2->c2).squarednorm());}
	double dc1p0() {return sqrt((u1->c1-u2->p0).squarednorm());}
	double dc1c2() {return sqrt((u1->c1-u2->c2).squarednorm());}
	double dc2p0() {return sqrt((u1->c2-u2->p0).squarednorm());}
	double dc2c1() {return sqrt((u1->c2-u2->c1).squarednorm());}
	double ac1p0p0() {return angle(u1->c1,u1->p0,u2->p0)*180.0/3.14159265;}
	double ap0p0c1() {return angle(u1->p0,u2->p0,u2->c1)*180.0/3.14159265;}
	double tc2c1p0p0() {return torsion(u1->c2,u1->c1,u1->p0,u2->p0)*180.0/3.14159265;}
	double tc1p0p0c1() {return torsion(u1->c1,u1->p0,u2->p0,u2->c1)*180.0/3.14159265;}
	double tp0p0c1c2() {return torsion(u1->p0,u2->p0,u2->c1,u2->c2)*180.0/3.14159265;}
};

struct PolarAtomGroup {
	PolarAtomGroup() {;}

	bool inMC{false};

	std::string resname;
	int chainseq{-1};
	int resseq{-1};
	std::map<std::string,XYZ> crds;
	PolarAtomGroup(const std::map<std::string,XYZ>&cs, std::string rn,
			int cq, int rq):crds(cs),resname(rn),chainseq(cq),resseq(rq) {
		;
	}

	std::set<std::vector<std::string>> units;
	std::set<std::string> polaratoms;
	void findunits();
	std::vector<HBUnit> extractHBunits();

	XYZ center();
	LocalFrame buildlocalframe();
	void global2local(const LocalFrame &lf);

	XYZ center_all();
	void translate(const XYZ &lth);
	void rotate(const Rotation &rot);
};

struct PolarAtomGroupPair {
	std::shared_ptr<PolarAtomGroup> g1 { nullptr };
	std::shared_ptr<PolarAtomGroup> g2 { nullptr };
	PolarAtomGroupPair(std::shared_ptr<PolarAtomGroup>p1, std::shared_ptr<PolarAtomGroup>p2):g1(p1),g2(p2) {
		;
	}
	void extracthbunits(std::vector<std::shared_ptr<HBUnit>> &u1s, std::vector<std::shared_ptr<HBUnit>> &u2s);
	void extracthbpairs(std::vector<HBUnitPair> &hbps, bool adjusthbp=false);
};

class PolarAtomPacking {
public:
	PolarAtomPacking(const std::vector<std::vector<Residue>>&ress, double dc=4.5):dcutoff_(dc), chains_(ress) {
		;
	}

	void extracthbunitsinMC(std::vector<std::vector<std::shared_ptr<HBUnit>>> &us);
	void extracthbunitsinSC(std::vector<std::vector<std::shared_ptr<HBUnit>>> &us);
	std::vector<std::vector<std::shared_ptr<HBUnit>>> findhbps(bool adjusthbp=false);
	std::vector<HBUnitPair>& hbps() {return hbpairs_;}

	void extractpagsinMC(std::vector<std::shared_ptr<PolarAtomGroup>> &gs);
	void extractpagsinSC(std::vector<std::shared_ptr<PolarAtomGroup>> &gs);
	void findgpps();
	std::vector<PolarAtomGroupPair>& gpps() {return gpairs_;}

private:
	std::vector<std::vector<Residue>> chains_;
	double dcutoff_;
	std::vector<PolarAtomGroupPair> gpairs_;
	std::vector<HBUnitPair> hbpairs_;
};

}


#endif /* ALLATOM_POLARATOMPACKING_H_ */
